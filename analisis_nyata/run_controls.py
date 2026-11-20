import sys, json, random
sys.path.insert(0, '.')
from optlib import *
null = json.load(open('data/real/null_shuffle.json'))
name = sys.argv[1]
rng = random.Random(5)
V1 = dict(e5=P['e1'], e3=P['e2'], sp5=P['sp5'], sp3=P['sp3'], arm=19, ac30=True)
D = {
 'PC1_V1_paper_CIRC-V1': V1,
 'PC2_V1_BsmBI_free_exon': dict(V1, e5='GAGAGGCTACGGACTT'),
 'PC3_16-16_spacer20_noarm': dict(e5='GAGAGGCTACGGACTT', e3=P['e2'], sp5=P['sp5'], sp3=P['sp3'], arm=0, ac30=True),
 'V10_8-8': dict(e5=P['e1'][-8:], e3=P['e2'][:8], sp5=P['sp'].get('10', '') , sp3=P['sp']['10'], arm=0, ac30=False),
 'V11_3-2_min': dict(e5='CTT', e3='AA', sp5='', sp3='', arm=0, ac30=True),
 'V12_3-2_spacer_arm19': dict(e5='CTT', e3='AA', sp5=P['sp5'], sp3=P['sp3'], arm=19, ac30=True),
 'NC1_exon_1-1': dict(e5='T', e3='A', sp5='', sp3='', arm=0, ac30=True),
 'NC2_IRES_scrambled': dict(e5='CTT', e3='AA', sp5='', sp3='', arm=0, ac30=True, ires=rubric2.dinuc_shuffle(IRES, rng)),
 'NC3_ORF_nonsense_codon50': dict(e5='CTT', e3='AA', sp5='', sp3='', arm=0, ac30=True),
 'NC4_intron_scrambled': dict(e5='CTT', e3='AA', sp5='', sp3='', arm=0, ac30=True, intron=rubric2.dinuc_shuffle(INTRON, rng)),
 'NC5_arm_100nt': dict(e5='GAGAGGCTACGGACTT', e3=P['e2'], sp5=P['sp5'], sp3=P['sp3'], arm=100, ac30=True),
}
d = D[name]; orf = ORF0
if name.startswith('NC3'):
    orf = ORF0[:147] + 'TAA' + ORF0[150:]
asm, m, s2 = evaluate(d, orf, null)
s2['name'] = name; s2['circle_length'] = len(asm.circle)
json.dump(s2, open(f'ctrl/{name}.json', 'w'), indent=1, default=str)
print(name, s2['composite'], s2['verdict'], s2['gates_failed'])
