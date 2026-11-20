"""Ekspor GenBank (.gb) beranotasi dan FASTA per kandidat untuk SnapGene/Benchling.
Jalankan di folder kerja pipeline: PYTHONPATH=src python3 export_gb.py"""
import json, os, sys, datetime
sys.path.insert(0, '../src')
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from circhbsag.parts import load_parts
from circhbsag.construct import build
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
ORF = ''.join(l.strip() for l in open('out_real/orf_optimized.fa') if not l.startswith('>'))
S_PROT = ''.join(l.strip() for l in open('data/real/S_protein_AAB59969.faa') if not l.startswith('>'))
LABEL = {'T7_promoter': ('promoter', "promotor T7"), 'intron_3half': ('misc_feature', "intron Ana (grup I) utuh"), 'E2': ('misc_feature', "ekson 3' (E2)"),
         'spacer5': ('misc_feature', "spacer 5'"), 'HA5': ('misc_feature', "lengan homologi 5'"), 'AC30': ('misc_feature', "AC30 sebelum IRES"),
         'IRES': ('misc_feature', "IRES CVB3 AF169665.1"), 'ORF': ('CDS', "ORF protein S HBsAg (U95551.1, kodon dioptimasi)"), 'UTR3': ('3UTR', "3' UTR HBA1 (paper)"),
         'HA3': ('misc_feature', "lengan homologi 3'"), 'spacer3': ('misc_feature', "spacer 3'"), 'E1': ('misc_feature', "ekson 5' (E1)")}
SRC = {'intron_3half': 'Shen dkk. 2025 Supp Data 1', 'E2': 'Shen dkk. 2025', 'E1': 'Shen dkk. 2025', 'spacer5': 'Shen dkk. 2025', 'spacer3': 'Shen dkk. 2025',
       'HA5': 'Shen dkk. 2025', 'HA3': 'Shen dkk. 2025', 'AC30': 'Shen dkk. 2025', 'UTR3': 'Shen dkk. 2025', 'IRES': 'GenBank AF169665.1', 'ORF': 'GenBank U95551.1 (protein S)'}
DESC = {'V14': 'ekson 11/11 nt, lengan 19 nt, tanpa spacer, AC30, ORF dimutasi sinonim', 'V15': 'ekson 16/16 nt (BsmBI dibuang), lengan 19 nt, tanpa spacer, AC30, ORF dimutasi sinonim',
        'V16': 'ekson minimal GCU/AA, spacer 20/18, lengan 19 nt, AC30, ORF dimutasi sinonim', 'V17': 'ekson minimal GCU/AA, spacer 20/18, tanpa lengan, AC30, ORF dimutasi sinonim',
        'V18': 'IRES-terbelah tanpa bekas: sambungan balik di dalam IRES AF169665.1 (posisi 322|323), ORF dimutasi sinonim'}
def find(seq, q, circular=False):
    s = seq + (seq[:len(q)] if circular else ''); i = s.find(q)
    return i if i >= 0 else None
def feat(a, b, kind, label, note=None, strand=1, extra=None):
    q = {'label': [label]}
    if note: q['note'] = [note]
    if extra: q.update(extra)
    return SeqFeature(FeatureLocation(a, b, strand=strand), type=kind, qualifiers=q)
def split_asm():
    sys.path.insert(0, '.')
    from optlib import build_split
    orf = json.load(open('final_orfs.json'))['V18']; af = json.load(open('data/real/split_sites_AF.json'))['403']['af_pos']
    return build_split(af, orf)[0]
def build_rec(v, circular):
    if v != 'V18':
        parts, prov = load_parts(f'config/parts.{v}.yaml', root='.')
    else:
        parts, prov = {}, {}
    if v != 'V18' and parts['ORF'].seq.upper() == 'ATGTAA': parts['ORF'].seq = ORF
    if v == 'V18':
        asm = split_asm(); coords = asm.circle_coords if circular else asm.precursor_coords
    else:
        asm = build(parts); coords = asm.circle_coords if circular else asm.precursor_coords
    seq = asm.circle if circular else asm.precursor
    desc = DESC.get(v) or (__import__('yaml').safe_load(open(f'config/parts.{v}.yaml'))['design'].get('note', '') if os.path.exists(f'config/parts.{v}.yaml') else '')
    rec = SeqRecord(Seq(seq), id=f'circHBsAg_{v}', name=f'circHBsAg_{v}'[:16],
        description=f"Kandidat {v} circRNA HBsAg (S, U95551.1) dengan IRES CVB3 AF169665.1, desain CIRC (Shen dkk. 2025): " + ('lingkaran hasil splicing (DNA setara RNA)' if circular else 'templat IVT: promotor T7 + prekursor') + f". {desc}")
    rec.annotations['molecule_type'] = 'DNA'; rec.annotations['topology'] = 'circular' if circular else 'linear'
    rec.annotations['data_file_division'] = 'SYN'; rec.annotations['date'] = '20-NOV-2026'; rec.annotations['organism'] = 'synthetic construct'; rec.annotations['source'] = 'synthetic construct'
    rec.features.append(SeqFeature(FeatureLocation(0, len(seq), strand=1), type='source', qualifiers={'organism': ['synthetic construct'], 'mol_type': ['other DNA'], 'note': ['in silico design; ' + ('circular (RNA product written as DNA)' if circular else 'linear DNA template')]}))
    rec.annotations['comment'] = ("In silico. Sekuens: Shen dkk. 2025 Nat Commun 16:7376 (Supplementary Data 1), GenBank AF169665.1 dan U95551.1. "
        "Ini rancangan, belum diuji di laboratorium. Letak IGS pada intron adalah perkiraan.")
    if v == 'V18':
        for name, (a, b) in coords.items():
            if ':' in name: continue
            kind, label = LABEL[name]
            if name == 'IRES': label = 'IRES CVB3 AF169665.1' + ('' if circular else " (bagian 3', posisi 323-660)")
            if name == 'E2': label = "ekson 3' (E2) = IRES 323-324"
            if name == 'E1': label = "ekson 5' (E1) = IRES 320-322"
            extra = {}
            if name == 'ORF': extra = {'codon_start': [1], 'transl_table': [1], 'product': ['HBsAg S protein (226 aa)'], 'translation': [S_PROT]}
            rec.features.append(feat(a - 1, b, kind, label, SRC.get(name, 'konsensus'), extra=extra))
        if circular:
            rec.features.append(feat(321, 322, 'misc_feature', "BSJ: ujung ekson 5' (IRES 322)|ekson 3' (IRES 323)", 'sambungan balik di dalam IRES; IRES utuh terbentuk kembali'))
        else:
            i = coords['intron_3half'][0]; rec.features.append(feat(i - 1, i, 'misc_feature', 'awal transkripsi (+1 G)', 'GG pertama intron'))
            rec.features.append(feat(len(seq) - 322, len(seq), 'misc_feature', "IRES posisi 1-322 (bagian 5')", 'diakhiri ekson 5\' (E1)'))
    for name, (a, b) in (coords.items() if v != 'V18' else []):
        if ':' in name: continue
        kind, label = LABEL[name]
        extra = {}
        if name == 'ORF': extra = {'codon_start': [1], 'transl_table': [1], 'product': ['HBsAg S protein (226 aa)'], 'translation': [S_PROT]}
        rec.features.append(feat(a - 1, b, kind, label, SRC.get(name, 'konsensus'), extra=extra))
    if v == 'V18': pass
    elif circular:
        n = len(seq); rec.features.append(feat(n - 1, n, 'misc_feature', 'BSJ: ujung ekson 5′|ekson 3′', 'sambungan balik; primer pengapit BSJ menargetkan lokasi ini'))
    else:
        i = coords['intron_3half'][0]; rec.features.append(feat(i - 1, i, 'misc_feature', "awal transkripsi (+1 G)", 'GG pertama intron'))
    # primer
    pf = f'prim/{v}.json'
    if os.path.exists(pf):
        p = json.load(open(pf))['primers']; circ = asm.circle; add = []
        d = p['divergent']['pairs'][0]; add += [('Primer_konvergen_F (kontrol total RNA)', d['left_seq'], 1), ('Primer_konvergen_R (kontrol total RNA)', d['right_seq'], -1)]
        b = p['bsj']['pairs'][0] if p['bsj']['pairs'] else None
        if b: add += [('Primer_BSJ_F (pengapit sambungan)', b['left_seq'], 1), ('Primer_BSJ_R (pengapit sambungan)', b['right_seq'], -1)]
        if p['rnase_h']['probes']: add += [('Probe_RNaseH', p['rnase_h']['probes'][0]['probe_dna_5to3'], -1)]
        for name, q, strand in add:
            tgt = q if strand == 1 else str(Seq(q).reverse_complement()); ref = seq
            i = find(ref, tgt, circular)
            if i is None and name == 'Primer_BSJ_F': continue
            if i is None: continue
            if i + len(tgt) <= len(ref):
                loc = FeatureLocation(i, i + len(tgt), strand=strand)
                rec.features.append(SeqFeature(loc, type='primer_bind', qualifiers={'label': [name], 'note': [f'{q} (5′→3′)']}))
            else:
                # melintasi origin: dipecah menjadi dua fitur agar valid di semua pembaca
                k = len(ref) - i
                rec.features.append(SeqFeature(FeatureLocation(i, len(ref), strand=strand), type='primer_bind', qualifiers={'label': [name + ' (bagian 1/2, melintasi BSJ)'], 'note': [f'{q} (5′→3′)']}))
                rec.features.append(SeqFeature(FeatureLocation(0, len(tgt) - k, strand=strand), type='primer_bind', qualifiers={'label': [name + ' (bagian 2/2)'], 'note': [f'{q} (5′→3′)']}))
                continue
    return rec
from Bio import SeqIO
allrec = []
for v in ['V%d' % i for i in range(1, 19)]:
    for circ in (False, True):
        rec = build_rec(v, circ); suffix = 'lingkaran' if circ else 'templat_T7'
        SeqIO.write(rec, f'{OUT}/circHBsAg_{v}_{suffix}.gb', 'genbank'); SeqIO.write(rec, f'{OUT}/circHBsAg_{v}_{suffix}.fasta', 'fasta')
        allrec.append(rec)
SeqIO.write(allrec, f'{OUT}/semua_kandidat.fasta', 'fasta')
print(len(allrec), 'rekaman')
