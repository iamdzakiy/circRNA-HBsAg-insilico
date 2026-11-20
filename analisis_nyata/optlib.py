import json, sys, time, random, yaml, hashlib
sys.path.insert(0, '../src')
from circhbsag import construct, score, rubric2, fold, access
from circhbsag.seqio import read_fasta, translate, revcomp
P = json.load(open('data/real/paper_parts.json'))
IRES = read_fasta('data/example/ires_cvb3_AF169665.fa')[0].seq.upper()
PROT = read_fasta('data/real/S_protein_AAB59969.faa', alphabet='protein')[0].seq if False else None
def _prot():
    t = open('data/real/S_protein_AAB59969.faa').read().split('\n', 1)[1].replace('\n', '').strip().upper()
    return t
PROT = _prot()
ORF0 = read_fasta('out_real/orf_optimized.fa')[0].seq.upper()
INTRON = P['intron']
T7 = 'TAATACGACTCACTATA'
SP_AC = {int(k): v for k, v in P['sp'].items()}

def design_parts(d, orf):
    """d: e5 (seq), e3 (seq), sp5, sp3 (seq), arm (nt; 0 none; 19 = V1), ac30 (bool), utr(False)"""
    parts = {
        'T7_promoter': {'seq': T7},
        'intron_3half': {'seq': d.get('intron', INTRON), 'features': {'IGS': (33, 40)}},
        'E2': {'seq': d['e3']}, 'spacer5': {'seq': d.get('sp5', '')},
        'IRES': {'seq': d.get('ires', IRES)}, 'ORF': {'seq': orf}, 'spacer3': {'seq': d.get('sp3', '')}, 'E1': {'seq': d['e5']},
    }
    arm = d.get('arm', 0)
    if arm:
        if arm == 19: a5, a3 = P['v1_ha5'], P['v1_ha3']
        else: a5, a3 = P['ha'][str(arm)]
        parts['HA5'] = {'seq': a5}; parts['HA3'] = {'seq': a3}
    if d.get('ac30', True): parts['AC30'] = {'seq': P['ac30']}
    out = {}
    for k, v in parts.items():
        out[k] = construct.Part(name=k, seq=v['seq'], features={kk: tuple(vv) for kk, vv in v.get('features', {}).items()})
    return out

def evaluate(d, orf, null, baseline_cache={}, full=True, temps=(37.0, 55.0)):
    parts = design_parts(d, orf)
    asm = construct.build(parts)
    m = score.measure(asm, temperatures=temps)
    key = (d['e3'], d['e5'])
    if key not in baseline_cache:
        baseline_cache[key] = rubric2.minimal_splice_baseline(INTRON, d['e3'], d['e5'], (33, 40), temps=temps)
    spl = max(len(d.get('sp5', '')), len(d.get('sp3', '')))
    ctx = dict(e5_seq=d['e5'], e3_seq=d['e3'], spacer_nt=min(len(d.get('sp5', '')), len(d.get('sp3', ''))) if d.get('sp5') or d.get('sp3') else 0,
               arm_nt=d.get('arm', 0), null=null, baseline=baseline_cache[key],
               ires_ok=asm.subseq('circle', 'IRES') == IRES, intron_ok=asm.subseq('precursor', 'intron_3half') == INTRON,
               orf_ok=(translate(asm.subseq('circle', 'ORF')).rstrip('*') == PROT and '*' not in translate(asm.subseq('circle', 'ORF'))[:-1]),
               cloning=d.get('cloning', 'golden_gate'))
    s2 = rubric2.score2(m, ctx)
    return asm, m, s2


def build_split(site_af, orf, e3n=2):
    """Desain IRES-terbelah tanpa bekas: BSJ jatuh di dalam IRES AF169665.1 sehingga IRES utuh terbentuk kembali."""
    j = site_af  # posisi 1-basis nukleotida terakhir ekson 5' pada IRES AF
    e5, e3 = IRES[j-3:j], IRES[j:j+e3n]
    ires3 = IRES[j:]            # mulai dari ekson 3' (E2) sampai ujung IRES
    pre_body = ires3 + orf + IRES[:j]   # ... E1 = IRES[j-3:j] di ujung
    precursor = T7 + INTRON + pre_body
    circle = IRES + orf          # dirotasi agar IRES kontigu (siklik ekuivalen dengan ires3+orf+IRES[:j])
    pc = {'T7_promoter': (1, 17), 'intron_3half': (18, 17+len(INTRON))}
    p0 = 17 + len(INTRON)
    pc['intron_3half:IGS'] = (p0 - len(INTRON) + 33 + 1, p0 - len(INTRON) + 40)
    pc['E2'] = (p0 + 1, p0 + e3n)
    pc['IRES'] = (p0 + 1, p0 + len(ires3))
    pc['ORF'] = (p0 + len(ires3) + 1, p0 + len(ires3) + len(orf))
    pc['E1'] = (len(precursor) - 2, len(precursor))
    cc = {'IRES': (1, len(IRES)), 'ORF': (len(IRES) + 1, len(IRES) + len(orf))}
    parts = {k: construct.Part(name=k, seq=v) for k, v in dict(T7_promoter=T7, intron_3half=INTRON, E2=e3, IRES=IRES, ORF=orf, E1=e5).items()}
    asm = construct.Assembly(precursor=precursor, circle=circle, parts=parts, precursor_coords=pc, circle_coords=cc, bsj_position=j)
    return asm, e5, e3

def evaluate_split(site_paper, site_af, eff, orf, null, temps=(37.0, 55.0)):
    asm, e5, e3 = build_split(site_af, orf)
    m = score.measure(asm, temperatures=temps)
    base = rubric2.minimal_splice_baseline(INTRON, e3, e5, (33, 40), temps=temps)
    prot_ok = translate(asm.subseq('circle', 'ORF')).rstrip('*') == PROT
    ctx = dict(e5_seq=e5, e3_seq=e3, spacer_nt=0, arm_nt=0, null=null, baseline=base,
               ires_ok=asm.subseq('circle', 'IRES') == IRES and (asm.subseq('precursor', 'IRES') + asm.precursor[-3:] ).count('') >= 0 and (IRES[:site_af] == asm.precursor[-site_af:]),
               intron_ok=asm.subseq('precursor', 'intron_3half') == INTRON, orf_ok=prot_ok,
               splice_override=eff / 100.0, cloning='golden_gate')
    return asm, m, rubric2.score2(m, ctx)
