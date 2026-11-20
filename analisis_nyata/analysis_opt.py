import os, json, glob, random, statistics as st, sys
sys.path.insert(0, '../src')
from circhbsag import rubric2
W = rubric2.WEIGHTS2
K = list(W)
def L(f): return json.load(open(f))
ctrl = {d['name']: d for d in map(L, glob.glob('ctrl/*.json'))}
scan = {d['name']: d for d in map(L, glob.glob('scan/S_*.json'))}
fin = {d['name']: d for d in map(L, glob.glob('final/*.json'))}
out = {}
# --- kontrol
out['controls'] = [{'name': n, 'composite': d['composite'], 'verdict': d['verdict'], 'gates_failed': d['gates_failed'],
                    'e_total': d['criteria']['splicing_competence']['raw']['e_total']} for n, d in sorted(ctrl.items())]
# --- pindai
rows = []
for n, d in scan.items():
    ex, sp, arm = n[2:].split('_'); rows.append(dict(name=n, exon=ex, spacer=int(sp[2:]), arm=int(arm[3:]), composite=d['composite'], verdict=d['verdict'],
                 e_total=d['criteria']['splicing_competence']['raw']['e_total'], circle=d['circle_length'], sub={k: d['criteria'][k]['subscore'] for k in K}))
rows.sort(key=lambda r: -r['composite']); out['scan'] = rows
def marg(key):
    g = {}
    for r in rows: g.setdefault(r[key], []).append(r['composite'])
    return {str(k): [round(st.mean(v), 3), round(min(v), 3), round(max(v), 3), len(v)] for k, v in sorted(g.items(), key=lambda x: str(x[0]))}
out['scan_marginal'] = {k: marg(k) for k in ('exon', 'spacer', 'arm')}
out['scan_counts'] = {v: sum(r['verdict'] == v for r in rows) for v in ('GO', 'GO-BERSYARAT', 'NO-GO')}
# --- final
ID = {'SA_11-11_0_19_c85': 'V14', 'SA_16-16_0_19_c85': 'V15', 'SA_3-2GCT_20_19_c85': 'V16', 'SA_3-2GCT_20_0_c85': 'V17', 'SPLIT_403_sa_c85': 'V18'}
UNC = {'SA_11-11_0_19': 'V14', 'SA_16-16_0_19': 'V15', 'SA_3-2GCT_20_19': 'V16', 'SPLIT_403_sa': 'V18'}
base_of = {'V14': 'S_11-11_sp0_arm19', 'V15': 'S_16-16_sp0_arm19', 'V16': 'S_3-2GCT_sp20_arm19', 'V17': 'S_3-2GCT_sp20_arm0', 'V18': 'SPLIT_403'}
def jitter(sub, n=200, seed=7, j=0.5):
    rng = random.Random(seed); sc = []; vd = {}
    for _ in range(n):
        w = {k: rng.uniform(W[k] * (1 - j), W[k] * (1 + j)) for k in K}
        c = sum(w[k] * sub[k] for k in K) / sum(w.values()); sc.append(c)
        vd['GO' if c >= .70 else 'GO-BERSYARAT' if c >= .55 else 'NO-GO'] = vd.get('GO' if c >= .70 else 'GO-BERSYARAT' if c >= .55 else 'NO-GO', 0) + 1
    sc.sort(); return dict(p05=round(sc[int(.05 * n)], 3), median=round(sc[n // 2], 3), p95=round(sc[int(.95 * n)], 3), verdicts=vd)
def loo(sub):
    r = {}
    for k in K:
        ws = {x: W[x] for x in K if x != k}; r[k] = round(sum(ws[x] * sub[x] for x in ws) / sum(ws.values()), 3)
    return r
F = []
for fn, vid in ID.items():
    d = fin[fn]; sub = {k: d['criteria'][k]['subscore'] for k in K}
    b = (scan.get(base_of[vid]) or fin.get(base_of[vid]))
    F.append(dict(id=vid, name=fn, composite=d['composite'], verdict=d['verdict'], gates_failed=d['gates_failed'], circle=d['circle_length'],
                  precursor=d.get('precursor_length'), sub=sub, e_total=d['criteria']['splicing_competence']['raw']['e_total'],
                  composite_before_mutagenesis=b['composite'] if b else None, jitter=jitter(sub), leave_one_out=loo(sub),
                  raw={k: d['criteria'][k]['raw'] for k in K}))
F.sort(key=lambda x: x['id']); out['final'] = F
out['split'] = []
for n in ('403', '406', '436', '443', '468'):
    d = fin[f'SPLIT_{n}']; out['split'].append(dict(site=int(n), af_pos=d['site']['af_pos'], ctx=d['site']['af_ctx'], eff=round(d['site']['eff'], 1), composite=d['composite'],
        sub={k: d['criteria'][k]['subscore'] for k in K}, ires37=d['criteria']['ires_integrity']['raw']['identity_37C']))
for n, key in (('403', 'SPLIT_403_sa_c85'),):
    if key not in fin: continue
    d = fin[key]; out['split'].append(dict(site=f'{n} (ORF termutasi)', composite=d['composite'], sub={k: d['criteria'][k]['subscore'] for k in K}, ires37=d['criteria']['ires_integrity']['raw']['identity_37C']))
out['unconstrained'] = [dict(id=v, composite=fin[n]['composite']) for n, v in UNC.items() if n in fin]
import sys as _s
from circhbsag import codon as cd
Wc = cd.relative_adaptiveness(cd.load_codon_usage())
orf0 = ''.join(l.strip() for l in open('out_real/orf_optimized.fa') if not l.startswith('>')).upper()
FO = json.load(open('final_orfs.json'))
out['orf_stats'] = {'ORF_awal': dict(cai=round(cd.cai(orf0, Wc), 3), gc=round(sum(c in 'GC' for c in orf0) / len(orf0), 3), nt=0, codons=0)}
for v, o in FO.items():
    out['orf_stats'][v] = dict(cai=round(cd.cai(o, Wc), 3), gc=round(sum(c in 'GC' for c in o) / len(o), 3), nt=sum(a != b for a, b in zip(o, orf0)), codons=sum(o[i:i+3] != orf0[i:i+3] for i in range(0, 681, 3)))
unc_orf = {}
for f, v in (('11-11_0_19', 'V14'), ('16-16_0_19', 'V15'), ('3-2GCT_20_19', 'V16'), ('split_403', 'V18')):
    if not os.path.exists(f'sa/{f}_s1.json'): continue
    o = L(f'sa/{f}_s1.json')['orf']; unc_orf[v] = dict(cai=round(cd.cai(o, Wc), 3), nt=sum(a != b for a, b in zip(o, orf0)))
out['orf_stats_unconstrained'] = unc_orf
# --- SA ringkas
out['sa'] = []
for f in sorted(glob.glob('sa/*_s[12].json')):
    d = L(f); out['sa'].append(dict(file=f, design=d['design'], J0=d['log'][0][1], J=round(d['best_fast'], 4), info0=d['log'][0][2], info=d['best_info'], nt=d['nt_changed'], prot=d['protein_ok'], steps=d['steps']))
# --- ablasi v1->v2
v1w = {'junction_accessibility': .25, 'ires_integrity': .20, 'dsrna_burden': .15, 'start_codon_accessibility': .15, 'circle_stability': .10, 'sequence_cleanliness': .10, 'motif_burden': .05}
def v1sub(v):
    x = L(f'varR/{v}/tahap0.json')['results']['decision']; return {c['name']: c['subscore'] for c in x['criteria']}, x['composite_score']
abl = []
for v, cn in (('V1', 'PC1_V1_paper_CIRC-V1'), ('V11', 'V11_3-2_min'), ('V12', 'V12_3-2_spacer_arm19')):
    if not os.path.exists(f'varR/{v}/tahap0.json'): continue  # ablasi v1->v2 memerlukan putaran rubrik v1 (tidak dijalankan ulang)
    s1, c1 = v1sub(v); s2 = {k: ctrl[cn]['criteria'][k]['subscore'] for k in K}
    cur = dict(s1); steps = [('rubrik v1', round(sum(v1w[k] * cur[k] for k in v1w), 3))]
    for lab, key in (('IRES diukur pada 37 °C', 'ires_integrity'), ('aksesibilitas situs splice relatif', 'junction_accessibility'), ('stabilitas terhadap nol acak', 'circle_stability')):
        cur[key] = s2[key]; steps.append((lab, round(sum(v1w[k] * cur[k] for k in v1w), 3)))
    steps.append(('kriteria splicing ditambah, bobot v2', round(ctrl[cn]['composite'], 3)))
    abl.append(dict(id=v, steps=steps))
out['ablation'] = abl
json.dump(out, open('analisis_opt.json', 'w'), indent=1, default=str)
for f in F: print(f['id'], f['composite'], f['verdict'], f['circle'], f['jitter'], 'before', f['composite_before_mutagenesis'])
print(out['scan_counts'], out['scan_marginal'])
for a in abl: print(a)
