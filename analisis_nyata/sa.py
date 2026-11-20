"""Mutagenesis in silico: substitusi kodon sinonim pada ORF dengan simulated annealing (objektif cepat, MFE)."""
import sys, json, math, random, time
sys.path.insert(0, '.')
from optlib import *
from circhbsag import codon as cd, motifs
from circhbsag.seqio import CODON_TABLE, SYN
import RNA

W = cd.relative_adaptiveness(cd.load_codon_usage())
import os
CAI_MIN = float(os.environ.get('CAI_MIN', '0'))
NULL = json.load(open('data/real/null_shuffle.json'))
_ramp = rubric2._ramp

def mfe_circ(seq, T=37.0):
    md = RNA.md(); md.circ = 1; md.temperature = T
    fc = RNA.fold_compound(seq.replace('T', 'U'), md)
    st, e = fc.mfe(); return st, e

def fast_objective(d, orf, ires_alone_struct):
    """Skor cepat (MFE saja) yang meniru kriteria rubrik v2 yang peka terhadap ORF."""
    asm = build_split(d['split'], orf)[0] if 'split' in d else construct.build(design_parts(d, orf))
    circ = asm.circle
    st, e = mfe_circ(circ)
    ia, ib = asm.circle_coords['IRES']
    ident = fold.subsequence_structure_identity(st, ires_alone_struct, ia - 1)
    ds = access.longest_dsrna(st)['longest_bp']
    oa, _ = asm.circle_coords['ORF']
    win = (max(1, oa - 6), min(len(circ), oa + 14))
    a = access.accessibility(circ, {'s': win}, ulength=10, window=80, temperature_c=37.0)['s']['mean']
    unp = sum(1 for c in st[win[0] - 1:win[1]] if c == '.') / (win[1] - win[0] + 1)
    pre = motifs.summarize(asm.precursor)
    t2 = pre['restriction'].get('BsaI', 0) + pre['restriction'].get('BsmBI_Esp3I', 0)
    hard = (t2 == 0 and not pre['t7_terminator'] and not motifs.upstream_atg_in_frame(circ, oa - 1, within=60))
    z = ((e / len(circ)) - NULL['mfe_per_nt_mean']) / NULL['mfe_per_nt_sd']
    hp = max(len(x) for x in __import__('re').findall(r'(A+|C+|G+|T+)', asm.precursor))
    gcc = sum(c in 'GC' for c in circ) / len(circ)
    band = 1.0 if 0.35 <= gcc <= 0.70 else max(0.0, 1 - ((0.35 - gcc) if gcc < 0.35 else (gcc - 0.70)) / 0.10)
    clean = 0.4 + 0.2 * (1 - _ramp(len(pre['gc_violations']), 0, 10)) + 0.2 * (1 - _ramp(hp, 6, 10)) + 0.2 * band
    m6a = pre['m6a_density_per_100nt']
    mot = 0.6 * (1 - _ramp(m6a, 1.0, 4.0)) + 0.4
    J = (0.15 * _ramp(ident, 0.70, 0.90) + 0.10 * (1 - _ramp(ds, 18, 45)) + 0.10 * (0.5 * _ramp(a, 0.01, 0.35) + 0.5 * _ramp(unp, 0.05, 0.60))
         + 0.10 * 0.5 * _ramp(-z, 0, 6) + 0.10 * clean + 0.05 * mot)
    return (J if hard else J - 1.0), dict(J=round(J, 4), ident=round(ident, 3), ds=ds, start_plfold=round(a, 3), start_unp=round(unp, 2), z_mfe=round(z, 2), hard_ok=hard, t2=t2)

def mutate(orf, rng):
    cods = [orf[i:i + 3] for i in range(0, len(orf), 3)]
    for _ in range(50):
        k = rng.randrange(1, len(cods) - 1)      # tanpa Met awal dan stop
        aa = CODON_TABLE[cods[k]]
        alts = [c for c in SYN[aa] if c != cods[k] and W[c] >= 0.30]
        if alts:
            cods[k] = rng.choice(alts); return ''.join(cods)
    return orf

def anneal(d, orf0, steps, seed, T0=0.02, T1=0.002, tag=''):
    rng = random.Random(seed)
    asm = build_split(d['split'], orf0)[0] if 'split' in d else construct.build(design_parts(d, orf0))
    ia, ib = asm.circle_coords['IRES']
    ires_alone = fold.fold(asm.circle[ia - 1:ib], circular=False, temperature_c=37.0).mfe_structure
    cur, info = fast_objective(d, orf0, ires_alone); best = (cur, orf0, info); orf = orf0
    log = [(0, cur, info)]
    t0 = time.time()
    for s in range(1, steps + 1):
        T = T0 * (T1 / T0) ** (s / steps)
        cand = mutate(orf, rng)
        if cd.cai(cand, W) < CAI_MIN:
            continue
        sc, inf = fast_objective(d, cand, ires_alone)
        if sc >= cur or rng.random() < math.exp((sc - cur) / T):
            orf, cur, info = cand, sc, inf
            if cur > best[0]: best = (cur, orf, info)
        if s % 25 == 0:
            log.append((s, cur, info)); print(tag, s, round(cur, 4), info, 'best', round(best[0], 4), round(time.time() - t0), 's', flush=True)
    return best, log

if __name__ == '__main__':
    name = sys.argv[1]; steps = int(sys.argv[2]); seed = int(sys.argv[3])
    from scan_designs import design
    if name.startswith('split'):
        d = {'split': json.load(open('data/real/split_sites_AF.json'))[name.split(':')[1]]['af_pos']}
        best, log = anneal(d, ORF0, steps, seed, tag=name)
        out = dict(design=name, seed=seed, steps=steps, best_fast=best[0], best_info=best[2], orf=best[1], nt_changed=sum(a != b for a, b in zip(ORF0, best[1])), protein_ok=translate(best[1]).rstrip('*') == PROT, log=[(a, round(b, 4), c) for a, b, c in log])
        json.dump(out, open(f"sa/{name.replace(':','_')}_s{seed}.json", 'w'), indent=1, default=str); print('saved'); sys.exit()
    ex, L, arm = name.split('|'); d = design(ex, int(L), int(arm))
    best, log = anneal(d, ORF0, steps, seed, tag=name)
    nmut = sum(a != b for a, b in zip(ORF0, best[1]))
    out = dict(design=name, seed=seed, steps=steps, start=fast_objective(d, ORF0, fold.fold(construct.build(design_parts(d, ORF0)).circle[300:960] if False else IRES, circular=False).mfe_structure)[1] if False else None,
               best_fast=best[0], best_info=best[2], orf=best[1], nt_changed=nmut, protein_ok=translate(best[1]).rstrip('*') == PROT, log=[(a, round(b, 4), c) for a, b, c in log])
    fn = f"sa/{name.replace('/', '-').replace('|', '_')}_s{seed}.json"
    json.dump(out, open(fn, 'w'), indent=1, default=str)
    print('saved', fn, best[0], nmut)
