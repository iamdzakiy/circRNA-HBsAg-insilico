"""Primer, probe, dan simulasi BSJ untuk satu kandidat (V1-V17 dari config, V18 dari desain IRES-terbelah)."""
import sys, json
sys.path.insert(0, '.')
from optlib import *
from circhbsag import primers, bsj_stats
v = sys.argv[1]
if v == 'V18':
    orf = json.load(open('final_orfs.json'))['V18']; af = json.load(open('data/real/split_sites_AF.json'))['403']['af_pos']
    asm = build_split(af, orf)[0]
else:
    import yaml
    from circhbsag import parts as pm
    parts, meta = pm.load_parts(f'config/parts.{v}.yaml')
    orf = json.load(open('final_orfs.json')).get(v, ORF0)
    parts['ORF'] = construct.Part(name='ORF', seq=orf)
    asm = construct.build(parts)
res = primers.design_all(asm)
sim = bsj_stats.run_simulation(asm.circle, asm.precursor, n_circular=1500, n_linear=1500, read_len=150, anchor=20, seed=20260101)
json.dump(dict(id=v, circle=asm.circle, precursor=asm.precursor, primers=res, bsj_simulation=sim), open(f'prim/{v}.json', 'w'), default=str)
print(v, len(asm.circle), len(asm.precursor), sim['support']['n_candidates'], sim['support']['confusion_on_simulated_origin'], flush=True)
