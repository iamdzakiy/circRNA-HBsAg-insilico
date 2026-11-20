# circhbsag — pipeline in silico circRNA-HBsAg (PIE, IRES CVB3)

Pipeline reprodusibel untuk **Tahap 0 (keputusan GO/NO-GO sebelum masuk
laboratorium)** dan **Tahap 1 (desain serta optimasi in silico)** konstruk RNA
sirkular pembawa ORF HBsAg:

```
prekursor:  T7 | GGG | HA5 | intron 3' PIE | E2 | spacer | IRES | ORF HBsAg | spacer | E1 | intron 5' PIE | HA3
lingkaran:                                   E2 | spacer | IRES | ORF HBsAg | spacer | E1
                                             ^-------------------------------------------^ sambungan balik (BSJ)
```

Pipeline menjawab satu pertanyaan praktis: **apakah konstruk ini layak
disintesis, atau sebaiknya didesain ulang dulu?** Jawabannya diberikan sebagai
keputusan GO / GO-BERSYARAT / NO-GO dengan gerbang keras dan skor komposit
berbobot yang seluruhnya dapat diaudit.

Lisensi kode: **MIT**. Lisensi setiap dependensi: lihat [LICENSE](LICENSE) dan
jalankan `circhbsag licenses`.

---

## Tiga aturan yang dipegang pipeline ini

1. **Tidak ada sekuens yang dikarang.** Sekuens referensi (IRES, intron, gen
   HBsAg) tidak diikutkan ke repositori. Anda menyediakannya sendiri, dengan
   nomor aksesi tercatat di `config/parts.yaml`. Bila berkasnya belum ada,
   pipeline **berhenti dengan galat** dan menyebut berkas mana yang kurang.
2. **Tidak ada angka hasil yang dihardcode.** Setiap nilai pada laporan
   dihitung saat dijalankan oleh ViennaRNA, primer3, atau kode di repositori
   ini. Versi tiap alat dicatat di setiap berkas JSON.
3. **Setiap laporan membawa label provenans.** Bila satu saja bagian masukan
   bertanda sintetis, seluruh laporan dilabeli `[DEMO — SEKUENS PLACEHOLDER]` dan
   menyatakan bahwa angkanya tidak boleh dikutip sebagai hasil ilmiah. Dengan
   masukan Anda sendiri, labelnya `[HASIL IN SILICO]`.

---

## Mulai cepat

```bash
# 1. Lingkungan (pilih salah satu)
mamba env create -f environment.yml && conda activate circhbsag
#   atau: python -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt
pip install -e .

# 2. Jalankan DEMO dengan sekuens placeholder (tanpa unduhan apa pun)
snakemake --cores 4 --config parts=config/parts.example.yaml outdir=results/demo

# 3. Lihat hasilnya
cat results/demo/SUMMARY.md
open results/demo/tahap0.html         # laporan lengkap, HTML mandiri

# 4. Satu berkas tabel lengkap (semua varian, uji, skor, primer, motif, kodon, sekuens)
#    results/demo/tables/circhbsag_hasil_lengkap.xlsx  (+ csv/ per sheet untuk Google Sheets)
python scripts/13_export_workbook.py --results results/demo   # otomatis ikut `snakemake`

# 5. Tes
pytest -q                              # 109 tes
```

Untuk **analisis sungguhan**: isi `data/reference/` (lihat
[data/reference/README.md](data/reference/README.md)), lengkapi nomor aksesi di
`config/parts.yaml`, lalu:

```bash
snakemake --cores 4                    # memakai config/parts.yaml (provenance_policy: strict)
```

### Apa arti angka pada jalankan DEMO

Jalankan demo di atas memakai `config/parts.example.yaml`, yang seluruh
bagiannya ditandai `synthetic: true` dan dibangkitkan secara deterministik oleh
`data/example/make_example.py`. **Tidak satu nukleotida pun berasal dari CVB3,
dari intron grup I Anabaena, maupun dari genom HBV.** Panjang tiap bagian
disamakan dengan konstruk sungguhan (protein demo 226 aa, sepanjang SHBs;
lingkaran 1407 nt) semata agar dimensi alurnya sebanding.

Karena itu angka yang keluar — skor komposit, MFE, aksesibilitas, CAI, jumlah
situs m6A, presisi filter BSJ — **mengilustrasikan bahwa alur kerjanya berjalan,
bukan sifat konstruk circRNA-HBsAg yang sebenarnya.** Angka itu tidak boleh
dikutip di dokumen proyek sebagai temuan. Setiap laporan yang dihasilkan
dari konfigurasi ini membawa label `[DEMO — SEKUENS PLACEHOLDER]` di kepala
berkas, dan langkah `workflow` pada CI memastikan label itu benar-benar muncul.

Untuk rujukan, jalankan demo pada commit ini menghasilkan keputusan
**GO-BERSYARAT dengan skor komposit 0.6200** (tidak ada gerbang keras yang
gagal); subskor terendah adalah `dsrna_burden` 0.185 dan
`start_codon_accessibility` 0.336. Pola itu wajar untuk sekuens acak, yang tidak
punya alasan apa pun untuk menjaga kodon awal tetap terbuka — dan justru itulah
yang ditunjukkan: mekanisme penilaiannya bekerja dan menghukum konstruk yang
belum dioptimasi. Angka pada sekuens Anda sendiri akan berbeda.

Verdikt GO-BERSYARAT itu sendiri memerintahkan siklus optimasi, dan siklus itu
dijalankan pada konstruk demo yang sama:

```bash
circhbsag optimize --parts config/parts.example.yaml --out results/demo \
                   --n-trials 24 --codon-iterations 3000 --top 24
```

40 trial Optuna (NSGA-II, 6 objektif) menghasilkan 24 titik di front Pareto.
Setelah front itu diperingkat ulang dengan skor komposit Tahap 0, konstruk
terbaik mencapai **0.7437 (GO)** — naik dari 0.6922 pada konstruk hasil redesain
Gate 0 — dengan spacer 50/30 nt dari pustaka `polyAC` dan lengan homologi 9 nt. Jadi loop
NO-GO/GO-BERSYARAT → optimasi → penilaian ulang memang tertutup dan terbukti
mampu memindahkan konstruk melewati ambang GO. Arah parameternya (spacer tak
berstruktur memperbaiki aksesibilitas sambungan) wajar secara termodinamika,
tetapi nilai optimumnya **tetap khusus untuk sekuens placeholder ini** dan harus
dihitung ulang pada sekuens sungguhan.

### Antarmuka web

```bash
streamlit run app/streamlit_app.py
```

Tempel sekuens tiap bagian, jalankan skor GO/NO-GO, lihat struktur lingkaran,
MFE, aksesibilitas, motif m6A, dan primer, lalu unduh laporan HTML/JSON.
Situs statis penjelasan metode ada di `docs/` (siap dipakai GitHub Pages).

---

## Apa yang dihitung

### Tahap 0 — kelayakan sirkularisasi

| Pengukuran | Alat | Keluaran |
|---|---|---|
| Lipatan sirkular pada 37 dan 55 °C | ViennaRNA `RNAfold --circ` (via API Python) | MFE, MFE/nt, energi bebas ensemble, keragaman ensemble, frekuensi MFE |
| Probabilitas pasangan basa di sambungan balik | partition function ViennaRNA | P(tak-berpasangan) jendela ±10 nt di BSJ, melingkar melewati batas |
| Aksesibilitas situs splice 5'/3' dan IGS (P1/P10) | `RNAplfold` | rerata dan minimum P(u nt tak-berpasangan) per wilayah |
| Integritas domain IRES | lipatan IRES sendiri vs di dalam konstruk | fraksi posisi yang mempertahankan status berpasangan/tidak |
| dsRNA beruntun | penguraian heliks dari struktur | panjang heliks terpanjang, daftar pelanggaran ≥ ambang |
| Keterbukaan start codon | `RNAplfold` + lipatan global | aksesibilitas jendela −6..+14 relatif AUG |
| Situs Tipe IIS internal (BsaI, BsmBI/Esp3I, SapI, BbsI, AarI) | pemindai regex dua untai | posisi tiap situs |
| Terminator T7 prematur | inti Tphi + pola intrinsik (batang GC + runtun U) | posisi tiap kecocokan |
| m6A RRACH | regex konsensus | jumlah, kerapatan per 100 nt, posisi adenosin sasaran |
| GC jendela, homopolimer | jendela 50 nt langkah 10 nt | pelanggaran rentang 30–70%, homopolimer terpanjang |

### Keputusan: gerbang keras lalu skor berbobot

**Lima gerbang keras.** Satu pelanggaran memaksa NO-GO, karena kegagalannya
tidak bisa dikompensasi kriteria lain:

| Gerbang | Alasan |
|---|---|
| `no_internal_typeIIS` | situs BsaI/BsmBI internal akan dipotong saat perakitan Golden Gate |
| `no_premature_T7_terminator` | transkrip prekursor terpotong sebelum intron 5' |
| `orf_intact_and_unique_start` | ATG di awal, panjang kelipatan 3, tanpa ATG sekerangka dalam 60 nt hulu |
| `dsrna_below_hard_limit` | dsRNA ≥ 45 bp memicu PKR/RIG-I dan menekan translasi |
| `circle_length_within_limit` | batas praktis sintesis gen (≤ 3,0 kb) |

**Tujuh kriteria berbobot** (jumlah bobot = 1,00), masing-masing dinormalkan
ke [0, 1] dengan fungsi yang didokumentasikan di `src/circhbsag/score.py`:

| Kriteria | Bobot | Mengapa sebesar itu |
|---|---:|---|
| `junction_accessibility` | 0,25 | penentu paling langsung sirkularisasi: situs sambungan dan P1/P10 harus dapat dicapai |
| `ires_integrity` | 0,20 | IRES yang lipatannya rusak membuat lingkaran terbentuk tetapi tidak tertranslasi |
| `dsrna_burden` | 0,15 | imunogenisitas intrinsik menentukan kegunaan sebagai vaksin, bukan hanya hasil sintesis |
| `start_codon_accessibility` | 0,15 | start codon tertutup struktur menekan inisiasi |
| `circle_stability` | 0,10 | lipatan harus kompak dan kokoh, tetapi tidak sampai rapat berlebihan |
| `sequence_cleanliness` | 0,10 | sebagian besar isinya sudah ditangani gerbang keras, jadi bobotnya kecil |
| `motif_burden` | 0,05 | m6A dan motif imunostimulan penting tetapi paling mudah diperbaiki lewat optimasi kodon |

**Ambang:** skor ≥ 0,70 dan tanpa gerbang keras gagal → **GO**;
0,55 ≤ skor < 0,70 → **GO-BERSYARAT** (jalankan siklus optimasi dulu);
selain itu **NO-GO**.

Semua bobot dan ambang ditandai `[PERKIRAAN]` — ditetapkan pipeline ini, bukan
diambil dari pustaka — dan dapat diubah lewat config. Karena itu setiap laporan
menyertakan **uji sensitivitas**: 200 set bobot dengan jitter ±50% dijalankan
ulang, dan distribusi skor serta keputusannya dilaporkan, sehingga pembaca tahu
apakah keputusan bergantung pada pilihan bobot.

### Tahap 1 — desain dan optimasi

- **Optimasi kodon** ORF HBsAg: hill-climbing stokastik dengan seed tetap atas
  objektif gabungan CAI, GC total, GC jendela, situs Tipe IIS, motif
  imunostimulan/destabilisasi, dan kerapatan m6A RRACH.
  *LinearDesign tidak dipakai* karena lisensinya membatasi pemakaian komersial
  dan redistribusi; optimizer di sini seluruhnya MIT.
- **Templat DNA terpermutasi** + pemeriksaan promotor T7 dan GGG awal transkrip.
- **Primer reverse IVT dengan 2'-O-metil** pada dua nukleotida 5' untuk menekan
  penambahan nukleotida non-templat oleh T7.
- **Primer divergen** (dirancang pada titik antipodal BSJ, sehingga tidak
  memberi produk pada templat linear) dan **primer yang menumpang BSJ**,
  keduanya lewat primer3.
- **Probe RNase H** dipilih dengan filter Tm/GC/homopolimer + aksesibilitas
  RNAplfold pada RNA lingkaran, berjarak ≥ 40 nt dari BSJ.
- **Filter statistik BSJ ala Szabo**: simulasi bacaan (seed tetap) → pemanggilan
  jangkar dua sisi → **uji keunikan** (bacaan yang dapat dijelaskan referensi
  linear dikeluarkan) → **posterior Bayes** per bacaan dari panjang overlap,
  jumlah mismatch, dan keunikan, lalu digabung ke tingkat sambungan.
  Perintah STAR chimeric yang setara dicetak untuk validasi silang bila STAR
  terpasang.
- **Homologi off-target**: BLAST+ bila tersedia, atau filter k-mer 21 nt
  terhadap FASTA lokal; bila keduanya tidak ada, status **SKIPPED** eksplisit.
- **Optimasi multi-objektif** (Optuna NSGA-II, 6 objektif) atas panjang dan
  komposisi spacer, panjang lengan homologi, serta seed/bobot optimasi kodon.
  Hasilnya front Pareto, lalu diperingkat ulang dengan skor komposit Tahap 0.

---

## Struktur repositori

```
src/circhbsag/        14 modul: seqio, motifs, fold, access, codon, construct,
                      parts, score, primers, bsj_stats, offtarget, optimize,
                      report, cli
config/parts.yaml     konfigurasi sungguhan (provenance_policy: strict)
config/parts.example.yaml   konfigurasi uji dengan data sintetis
config/config.yaml    parameter alur kerja (seed, suhu, jumlah trial, dll)
Snakefile             8 aturan, dari optimasi kodon sampai ringkasan
data/reference/       KOSONG — Anda isi sendiri; lihat README di dalamnya
data/example/         data sintetis + generatornya (deterministik, seed tetap)
tests/                102 tes pytest
notebooks/            notebook Tahap 0 dengan penjelasan langkah demi langkah
app/streamlit_app.py  antarmuka web
docs/                 situs statis (GitHub Pages) + DESIGN.md
.github/workflows/    CI: tes, alur ujung-ke-ujung, pengaman provenans
```

## Perintah CLI

```bash
circhbsag licenses                  # versi + lisensi semua alat terpasang
circhbsag tahap0      --parts config/parts.yaml --out results
circhbsag codon       --parts config/parts.yaml --iterations 20000
circhbsag primers     --parts config/parts.yaml --orf-fasta results/orf_optimized.fa
circhbsag bsj-sim     --parts config/parts.yaml --n-circular 1500 --anchor 20
circhbsag offtarget   --parts config/parts.yaml --db refseq_rna
circhbsag optimize    --parts config/parts.yaml --n-trials 40
circhbsag codon-table --fasta cds_inang.fa      # ganti tabel kodon [PERKIRAAN]
```

`tahap0` keluar dengan kode **2** bila keputusannya NO-GO, sehingga alur
Snakemake berhenti dan keputusan NO-GO tidak pernah terlewat diam-diam.

## Reprodusibilitas

Seed tetap di seluruh langkah stokastik (`config/config.yaml: seed`), versi alat
dicatat di setiap laporan, dan CI memverifikasi bahwa data contoh dapat
dibangkitkan ulang bit-per-bit. Jalankan `snakemake --dag | dot -Tsvg` untuk
diagram alurnya.

## Keterbatasan yang perlu dinyatakan di dokumen proyek

- Semua keluaran adalah **prediksi** termodinamika dan heuristik sekuens, bukan
  pengukuran. Tidak satu pun menggantikan uji RNase R, RT-PCR divergen dengan
  sekuensing Sanger, atau uji RNase H.
- Model energi ViennaRNA **tidak memodelkan pseudoknot**, sehingga domain IRES
  dan struktur tersier intron grup I tidak terwakili sepenuhnya. Metrik
  integritas IRES di sini adalah **proksi struktural**, bukan prediksi fungsi
  translasi.
- Skor komposit adalah alat **triase** untuk mengurutkan kandidat sebelum
  sintesis. Ia belum dikalibrasi terhadap efisiensi sirkularisasi hasil
  laboratorium; korelasi skor dengan hasil lab harus dilaporkan sebagai
  eksploratif.
- Tabel penggunaan kodon bawaan ditandai `[PERKIRAAN]` dan harus diganti dengan
  unduhan terverifikasi atau dihitung dari CDS inang (`circhbsag codon-table`)
  sebelum nilai CAI dilaporkan sebagai angka absolut.
- Prediksi epitop B/T (BepiPred, NetMHCpan/NetMHCIIpan, VaxiJen) dan struktur
  (ColabFold) **tidak** termasuk dalam pipeline ini: alat-alat tersebut berjalan
  sebagai layanan daring atau butuh unduhan model besar dengan pendaftaran
  akademik. Lihat [docs/DESIGN.md](docs/DESIGN.md) untuk prosedur manualnya
  beserta status lisensi masing-masing.
