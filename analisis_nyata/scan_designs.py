import sys, json
sys.path.insert(0, '.')
from optlib import *
E1_16 = 'GAGAGGCTACGGACTT'   # E1 16 nt, situs BsmBI dibuang (6 nt ujung 5' diubah; 6 nt P1 di ujung 3' tetap)
def exons(tag):
    if tag == '3/2CTT': return 'CTT', 'AA'
    if tag == '3/2GCT': return 'GCT', 'AA'
    n = int(tag.split('/')[0]); return E1_16[-n:], P['e2'][:n]
def spacers(L):
    if L == 0: return '', ''
    if L == 20: return P['sp5'], P['sp3']
    return P['sp'][str(L)], P['sp'][str(L)]
def design(ex, L, arm):
    e5, e3 = exons(ex); s5, s3 = spacers(L)
    return dict(e5=e5, e3=e3, sp5=s5, sp3=s3, arm=arm, ac30=True)
GRID = [(ex, L, arm) for ex in ['3/2CTT', '3/2GCT', '6/6', '11/11', '16/16'] for L in [0, 20, 50] for arm in [0, 19, 34]]
def name(ex, L, arm): return f"S_{ex.replace('/', '-')}_sp{L}_arm{arm}"
if __name__ == '__main__':
    if sys.argv[1] == 'list':
        # prediksi analitik, keluarkan daftar yang layak dilipat penuh
        for ex, L, arm in GRID:
            d = design(ex, L, arm)
            sp = rubric2.predict_splicing(d['e5'], d['e3'], min(len(d['sp5']), len(d['sp3'])) if d['sp5'] else 0, arm)
            print(name(ex, L, arm), sp['e_total'])
    else:
        ex, L, arm = GRID[int(sys.argv[1])]
        null = json.load(open('data/real/null_shuffle.json'))
        asm, m, s2 = evaluate(design(ex, L, arm), ORF0, null)
        s2['name'] = name(ex, L, arm); s2['design'] = design(ex, L, arm); s2['circle_length'] = len(asm.circle)
        json.dump(s2, open(f'scan/{s2["name"]}.json', 'w'), indent=1, default=str)
        print(s2['name'], s2['composite'], s2['verdict'], s2['gates_failed'], flush=True)
