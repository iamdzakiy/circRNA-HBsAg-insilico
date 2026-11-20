import sys, json
sys.path.insert(0, '.')
from optlib import *
from scan_designs import design
null = json.load(open('data/real/null_shuffle.json'))
kind, arg = sys.argv[1], sys.argv[2]
if kind == 'sa':
    SEED = __import__('os').environ.get('SEED', '1')
    sa = json.load(open(f'sa/{arg}_s{SEED}.json')); ex, L, ar = sa['design'].split('|')
    d = design(ex, int(L), int(ar)); orf = sa['orf']
    asm, m, s2 = evaluate(d, orf, null); s2['design'] = d; s2['orf'] = orf; name = 'SA_' + arg + ('' if SEED == '1' else '_c85')
else:
    sites = json.load(open('data/real/split_sites_AF.json'))
    s = sites[arg]; orf = ORF0
    if len(sys.argv) > 3:
        orf = json.load(open(f"sa/{sys.argv[3]}_s{__import__('os').environ.get('SEED','1')}.json"))['orf']
    asm, m, s2 = evaluate_split(arg, s['af_pos'], s['eff'], orf, null); s2['site'] = dict(paper=arg, **s); s2['orf'] = orf
    name = f'SPLIT_{arg}' + (('_sa' if __import__('os').environ.get('SEED','1') == '1' else '_sa_c85') if len(sys.argv) > 3 else '')
s2['name'] = name; s2['circle_length'] = len(asm.circle); s2['precursor_length'] = len(asm.precursor)
json.dump(s2, open(f'final/{name}.json', 'w'), indent=1, default=str)
print(name, s2['composite'], s2['verdict'], s2['gates_failed'], flush=True)
