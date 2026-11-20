# Alur skrip, hasil antara, dan penilaian in silico

> Proyeksi berbasis literatur dan perhitungan in silico, bukan data eksperimen.

## 1. Satu skrip utuh?

Ada **dua pintu masuk ke alur yang sama**, keduanya memakai pustaka `src/circhbsag/`:

| Pintu | Perintah | Cocok untuk |
|---|---|---|
| Snakemake (berantai otomatis) | `snakemake --cores 4 --config parts=<parts.yaml> outdir=<folder>` | satu perintah, semua langkah, hanya menghitung ulang yang berubah |
| Skrip bernomor `scripts/01…13` | `python scripts/NN_nama.py --parts <parts.yaml> --results <folder>` | menjalankan satu langkah sendiri; tiap keluaran punya `*.params.json` (argumen, seed, versi alat, waktu) |

`rule all` menghasilkan: ORF teroptimasi, Tahap 0, primer, simulasi BSJ, off-target, lisensi, `SUMMARY.md`, dan **tabel xlsx lengkap**.
Di luar `rule all` (karena lama atau bergantung keputusan): `gate0` (sembilan uji + kunci desain), `optimize` (Optuna).

## 2. Di mana hasil antara tiap langkah

`<out>` = folder keluaran. Peta lengkapnya juga ada di sheet **Peta_alur** pada xlsx.

| Langkah | Berkas |
|---|---|
| Kodon | `<out>/orf_optimized.fa`, `<out>/codon.json` |
| Rakit | `<out>/tables/construct.fa`, `construct_features.tsv`, `precursor.gb` |
| Gate 0 | `<out>/reports/gate0.json`, `gate0.md`, `design_lock.json` |
| Lipatan, IRES, m6A, dsRNA | `<out>/tables/fold_analysis.tsv`, `ires_integrity.json`, `m6a_rrach_sites.tsv`, `dsrna_stems.json` |
| Primer, BSJ, off-target | `<out>/primers.json`, `bsj_simulation.json`, `offtarget_report.json` |
| Optuna | `<out>/optimization.json` |
| Redesain Gate 0 | `<out>/gate0_iterations.json`, `.md` |
| **Tabel lengkap** | `<out>/tables/circhbsag_hasil_lengkap.xlsx` dan `<out>/tables/csv/*.csv` |

Setiap JSON memuat label data, seed, versi alat, dan waktu UTC. Xlsx tidak menghitung angka baru, hanya meratakan berkas di atas.

## 3. Mana hasil sungguhan, mana demo

- **`results/demo/`** memakai sekuens acak (`config/parts.example.yaml`). Semua angkanya placeholder dan **tidak boleh dikutip**. Xlsx dari folder ini diberi label `[DEMO — SEKUENS PLACEHOLDER]` di sheet Ringkasan dan Sumber_dan_status.
- **Angka di TA v3/v4** (V1–V18, kontrol positif dan negatif) dihitung dengan **rubrik v2** di salinan kerja thread TA (`ta-docx/insilico`), bukan oleh repo ini. Bobot v2 disalin sebagai snapshot di sheet **Rubrik_v2_TA**. Setelah rerun dengan protein AAB59969.1 (genom U95551.1) selesai, kode v2 dimasukkan ke repo agar satu skrip menghasilkan semua angka TA.
- Konfigurasi sungguhan `config/parts.yaml` (strict) menolak berjalan sampai IRES, intron, ekson, lengan, dan spacer disediakan di `data/reference/`.

## 4. Penilaian (rubrik v1 di repo) dan bobot

Skor komposit = jumlah (bobot × subskor tiap kriteria, 0–1). **GO ≥ 0,70**, **GO-BERSYARAT ≥ 0,55**, di bawah itu NO-GO.
Satu gerbang keras yang gagal memaksa NO-GO tanpa memandang skor.

| Kriteria | Bobot | Alasan mekanistik |
|---|---:|---|
| junction_accessibility | 0,25 | situs sambungan dan P1/P10 menentukan langsung sirkularisasi |
| ires_integrity | 0,20 | IRES harus tetap terlipat agar translasi berjalan |
| dsrna_burden | 0,15 | dsRNA panjang memicu respons imun bawaan (PKR/RIG-I) |
| start_codon_accessibility | 0,15 | AUG harus terbuka agar ribosom dapat mulai |
| circle_stability | 0,10 | lipatan lingkaran kokoh |
| sequence_cleanliness | 0,10 | situs enzim, GC, homopolimer, terminator T7 |
| motif_burden | 0,05 | motif m6A/imunostimulan |

**Status bobot: `[PERKIRAAN]`.** Bobot dan ambang adalah keputusan desain pipeline, **bukan nilai dari literatur**; rujukan hanya dipakai untuk ambang tertentu (mis. dsRNA ≈ 30 bp, `[LITERATUR]` dengan kutipan di `content/references.yaml`). Konsekuensinya diuji dengan **uji sensitivitas**: 200 undian bobot dengan variasi ±50%, hasilnya di sheet **Sensitivitas_bobot**. Pada run demo, keputusan berpindah antara GO dan GO-BERSYARAT, artinya keputusan sensitif terhadap pilihan bobot dan harus dilaporkan demikian.
Rumus normalisasi tiap kriteria dan ambangnya ada di sheet **Bobot_dan_justifikasi** dan **Ambang_berlabel**.

**Kontrol.** Rubrik v1 di repo tidak punya kontrol positif/negatif. Kontrol (PC1–PC3, NC1–NC5) hanya ada pada rubrik v2 TA. Temuan v1: kontrol negatif sempat skor lebih tinggi daripada konstruk dasar, sehingga rubrik dinyatakan belum valid dan diganti v2 (lihat Rubrik_v2_TA dan logbook 19–20 Nov).

## 5. Menjalankan di macOS

```bash
# 0. Sekali saja: alat dasar (Homebrew: https://brew.sh)
brew install python@3.12 git
# 1. Ambil kode
git clone https://github.com/iamdzakiy/circRNA-HBsAg-insilico.git
cd circRNA-HBsAg-insilico
# 2. Lingkungan
python3.12 -m venv .venv && source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt && pip install -e .
# 3. Jalankan seluruh alur (demo)
snakemake --cores 4 --config parts=config/parts.example.yaml outdir=results/demo
# 4. Buka tabel lengkap
open results/demo/tables/circhbsag_hasil_lengkap.xlsx
```

Dependensi non-pip: **tidak ada** untuk jalur demo. ViennaRNA, primer3-py, openpyxl, dan Snakemake semua dipasang oleh `pip`.
Opsional: `blastn` (BLAST+, `brew install blast`) untuk off-target; NetMHCpan (unduh lokal, lisensi akademik) untuk skrip 09.

Masalah umum:

| Gejala | Penyebab dan solusi |
|---|---|
| `pip install ViennaRNA` gagal membangun | Pakai conda: `brew install --cask miniforge` lalu `mamba env create -f environment.yml && conda activate circhbsag` |
| Apple Silicon: wheel tidak ada | Pakai jalur conda di atas (osx-arm64 tersedia di conda-forge/bioconda) |
| `snakemake: command not found` | Aktifkan venv: `source .venv/bin/activate` |
| `ModuleNotFoundError: circhbsag` | Jalankan `pip install -e .` dari folder repo |
| `Missing input: ...` | Berkas acuan belum ada di `data/reference/`; pesan galat menyebut berkas dan cara mengisinya |
| `snakemake` tidak menghitung ulang | Hapus folder keluaran atau tambahkan `--forcerun all` |
| Gate 0 keluar dengan kode 2 | Itu keputusan FAIL-PINDAH-JALUR, bukan galat program; baca `reports/gate0.md` |
