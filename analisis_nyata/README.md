# analisis_nyata: varian V1–V18 dengan sekuens nyata

Isi folder ini adalah salinan alur kerja yang menghasilkan skor kandidat V1–V18, kontrol positif/negatif, pemindaian 45 desain,
mutagenesis kodon sinonim, dan primer, untuk protein **AAB59969.1** (genom **U95551.1**, CDS 157..837). Rubrik yang dipakai adalah
**rubrik v2** (`src/circhbsag/rubric2.py`, 8 kriteria). Hasil dibuat pada 20 Nov 2026 dan disimpan apa adanya.

> Proyeksi berbasis literatur dan perhitungan in silico, bukan data eksperimen. Skor dioptimasi terhadap rubrik yang sama dengan penilainya.

| Folder/berkas | Isi |
|---|---|
| `config/parts.V*.yaml` | konfigurasi bagian untuk tiap varian |
| `data/real/` | protein, parameter model splicing dari Shen dkk. 2025 (`[PERKIRAAN]` lintas konteks), situs split, null pengacakan |
| `var2/`, `ctrl/`, `scan/`, `sa/`, `final/`, `prim/` | hasil antara: varian V1–V13, kontrol, 45 desain, mutagenesis, kandidat akhir, primer dan simulasi BSJ |
| `final_orfs.json`, `analisis_opt.json` | ORF teroptimasi dan ringkasan analisis |
| `hasil_insilico_U95551.xlsx` | tabel ringkas hasil (20 sheet) |
| `genbank_U95551/` | file GenBank untuk SnapGene/Benchling |

## Menjalankan ulang

```bash
cd analisis_nyata
bash run_all.sh        # jam-an pada 4 inti; memakai ../src (paket circhbsag di repo)
```

## Catatan reprodusibilitas (jujur)

Hasil di folder ini dihitung dengan **definisi motif lama** (`polyU_6`). Paket `circhbsag` di repo sekarang memakai definisi proyek
BLOK 4.2 (`polyU_5` dan `polyA_7`, batas GC jendela 0,35–0,70 yang satu untuk optimasi kodon dan Gate 0). Menilai ulang **V1** dengan paket repo
memberi komposit **0,7611** (hasil tersimpan: 0,7839); verdict sama (NO-GO, gerbang Tipe IIS). Selisihnya ada pada `sequence_cleanliness`
dan `motif_burden`; kriteria lain identik. Karena itu angka tersimpan belum identik dengan paket terbaru. Pilihannya: (a) menjalankan ulang
`run_all.sh` dengan paket repo supaya semua angka konsisten, atau (b) menetapkan definisi lama sebagai definisi proyek. Belum diputuskan.
