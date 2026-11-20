#!/bin/bash
# Alur otomatis lengkap: protein -> ORF -> semua varian, kontrol, pindai, mutagenesis, kandidat akhir, primer -> JSON.
# Pemakaian: bash run_all.sh   (butuh ViennaRNA, primer3-py, biopython, pyyaml; jam-an pada 4 inti)
set -e
export PYTHONPATH=../src
P=${PARALLEL:-4}
mkdir -p out_real scan sa final ctrl var2 prim
printf "%s\n" PC1_V1_paper_CIRC-V1 PC2_V1_BsmBI_free_exon PC3_16-16_spacer20_noarm V10_8-8 V11_3-2_min V12_3-2_spacer_arm19 NC1_exon_1-1 NC2_IRES_scrambled NC3_ORF_nonsense_codon50 NC4_intron_scrambled NC5_arm_100nt | xargs -P$P -I{} sh -c '[ -f ctrl/{}.json ] || python3 run_controls.py {} > ctrl/{}.log 2>&1'
echo "[3/7 kontrol] selesai $(date +%T)"
seq 0 44 | xargs -P$P -I{} sh -c 'python3 scan_designs.py {} > scan/{}.log 2>&1'
echo "[4/7 pindai 45 desain] selesai $(date +%T)"
printf "%s\n" "11/11|0|19" "16/16|0|19" "3/2GCT|20|19" "3/2GCT|20|0" "split:403" | xargs -P$P -I{} sh -c 'CAI_MIN=0.85 python3 sa.py "{}" 300 2 > sa/$(echo "{}" | tr "/|:" "-__").log 2>&1'
echo "[5/7 mutagenesis SA] selesai $(date +%T)"
for a in 11-11_0_19 16-16_0_19 3-2GCT_20_19 3-2GCT_20_0; do SEED=2 python3 run_final.py sa $a & done; wait
for s in 403 406 436 443 468; do python3 run_final.py split $s & done; wait
SEED=2 python3 run_final.py split 403 split_403
python3 collect_final_orfs.py
python3 make_final_configs.py
echo "[6/7 kandidat akhir] selesai $(date +%T)"
printf "%s\n" V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 | xargs -P$P -I{} sh -c 'python3 prim_all.py {} > prim/{}.log 2>&1'
echo "[7/7 primer + BSJ] selesai $(date +%T)"
echo ALLDONE > ALLDONE
