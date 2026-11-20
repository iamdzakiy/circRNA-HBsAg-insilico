"""Skor rubrik v2 untuk varian V1-V13 langsung dari config/parts.Vn.yaml (ORF diganti ORF hasil optimasi kodon)."""
import sys, json, yaml
sys.path.insert(0, '.')
from optlib import *
v = sys.argv[1]
cfg = yaml.safe_load(open(f'config/parts.{v}.yaml'))
P_ = cfg['parts']; parts = {}
for k, x in P_.items():
    if k == 'IRES': seq = IRES
    elif k == 'ORF': seq = ORF0
    else: seq = x['inline']
    feats = {kk: tuple(vv) for kk, vv in (x.get('features') or {}).items()}
    parts[k] = construct.Part(name=k, seq=seq, features=feats)
asm = construct.build(parts)
m = score.measure(asm, temperatures=(37.0, 55.0))
g = lambda k: P_[k]['inline'] if k in P_ else ''
e5, e3 = g('E1'), g('E2')
null = json.load(open('data/real/null_shuffle.json'))
base = rubric2.minimal_splice_baseline(INTRON, e3, e5, (33, 40), temps=(37.0, 55.0))
ctx = dict(e5_seq=e5, e3_seq=e3, spacer_nt=min(len(g('spacer5')), len(g('spacer3'))) if g('spacer5') or g('spacer3') else 0,
           arm_nt=len(g('HA5')), null=null, baseline=base, ires_ok=asm.subseq('circle', 'IRES') == IRES,
           intron_ok=asm.subseq('precursor', 'intron_3half') == INTRON,
           orf_ok=(translate(asm.subseq('circle', 'ORF')).rstrip('*') == PROT), cloning='golden_gate')
s2 = rubric2.score2(m, ctx); s2['name'] = v; s2['circle_length'] = len(asm.circle); s2['precursor_length'] = len(asm.precursor)
s2['note'] = cfg['design'].get('note', '')
json.dump(s2, open(f'var2/{v}.json', 'w'), indent=1, default=str); print(v, s2['composite'], s2['verdict'], s2['gates_failed'], flush=True)
