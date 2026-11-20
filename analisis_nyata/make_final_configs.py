import json, yaml, sys
sys.path.insert(0, '.')
from optlib import P, INTRON
from scan_designs import design
PAPER = "Shen dkk. 2025, Nat Commun 16:7376, Supplementary Data 1 (MOESM3)"
FIN = {  # id: (sa file, design string, keterangan)
 'V14': ('11-11_0_19', '11/11|0|19', 'ekson 11/11 nt (situs BsmBI dibuang), lengan homologi 19 nt, tanpa spacer, AC30, ORF hasil mutagenesis sinonim'),
 'V15': ('16-16_0_19', '16/16|0|19', 'ekson 16/16 nt (situs BsmBI dibuang), lengan homologi 19 nt, tanpa spacer, AC30, ORF hasil mutagenesis sinonim'),
 'V16': ('3-2GCT_20_19', '3/2GCT|20|19', 'ekson minimal GCU dan AA, spacer 20/18, lengan 19 nt, AC30, ORF hasil mutagenesis sinonim'),
 'V17': ('3-2GCT_20_0', '3/2GCT|20|0', 'ekson minimal GCU dan AA, spacer 20/18, tanpa lengan, AC30, ORF hasil mutagenesis sinonim'),
}
for v, (saf, dstr, note) in FIN.items():
    ex, L, arm = dstr.split('|'); d = design(ex, int(L), int(arm)); orf = json.load(open('final_orfs.json'))[v]
    def part(seq, role, **kw):
        x = {'inline': seq, 'source': PAPER, 'role': role}; x.update(kw); return x
    parts = {'T7_promoter': {'inline': 'TAATACGACTCACTATA', 'source': 'Promotor T7 konsensus (literatur)'},
      'intron_3half': part(INTRON, 'intron Ana utuh 251 nt (grup I)', features={'IGS': [33, 40]}, note="[PERKIRAAN] letak IGS dari komplementaritas; bukan anotasi paper"),
      'E2': part(d['e3'], "ekson 3'"), 'spacer5': part(d['sp5'], "spacer 5'")}
    if int(arm):
        a5, a3 = (P['v1_ha5'], P['v1_ha3']) if int(arm) == 19 else P['ha'][arm]
        parts['HA5'] = part(a5, "lengan homologi 5'"); parts['HA3'] = part(a3, "lengan homologi 3'")
    parts['AC30'] = part(P['ac30'], 'AC30')
    parts['IRES'] = {'file': 'data/example/ires_cvb3_AF169665.fa', 'accession': 'AF169665.1', 'role': 'IRES CVB3 (GenBank AF169665.1)'}
    parts['ORF'] = {'inline': orf, 'note': 'ORF S U95551.1 CDS 157..837 (AAB59969.1); kodon sinonim dimutasi in silico (simulated annealing)'}
    parts['spacer3'] = part(d['sp3'], "spacer 3'"); parts['E1'] = part(d['e5'], "ekson 5'")
    cfg = {'provenance_policy': 'strict', 'design': {'ires_choice': 'CVB3 AF169665.1', 'note': note},
           'protein': {'file': 'data/real/S_protein_AAB59969.faa', 'accession': 'AAB59969.1 (U95551.1 CDS 157..837)', 'genotype': 'D (komputasi); subtipe ayw (catatan GenBank)'}, 'parts': parts}
    yaml.safe_dump(cfg, open(f'config/parts.{v}.yaml', 'w'), sort_keys=False, allow_unicode=True)
    print(v, 'ok')
