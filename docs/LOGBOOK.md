# Logbook

Catatan kerja kronologis. Satu entri per langkah yang menghasilkan sesuatu yang
bisa diperiksa: berkas, angka, atau keputusan. Entri ditulis saat langkahnya
selesai, bukan direkonstruksi dari ingatan di akhir.

Angka pada logbook ini dikutip apa adanya dari jalankan yang disebut. Kalau
sebuah angka berubah karena pipeline berubah, entri lamanya **tidak diperbaiki**
— itulah gunanya logbook: menunjukkan bahwa angkanya bergerak, dan kapan.

---

## 2026-11-16 — Kerangka pipeline Tahap 0 dan Tahap 1

Dibangun paket `circhbsag` dengan modul: `seqio`, `motifs`, `fold`, `access`,
`construct`, `codon`, `primers`, `bsj_stats`, `offtarget`, `optimize`, `score`,
`report`, `parts`, `cli`. Keputusan dibuat dua lapis: lima gerbang keras yang
memveto, dan tujuh kriteria berbobot yang jumlah bobotnya 1,00.

Dua bug nyata ditemukan oleh tes, bukan oleh pembacaan kode:

1. Gerbang keras `orf_intact_and_unique_start` memaksa NO-GO karena mencari ATG
   sekerangka di **seluruh** hulu ORF dan menemukan tujuh di dalam IRES. IRES
   pikornavirus memang memuat banyak AUG, jadi gerbang itu akan menggagalkan
   konstruk yang sah. Diperbaiki menjadi jendela 60 nt, dan ATG sekerangka yang
   jauh dilaporkan sebagai informasi, bukan pelanggaran.
2. Filter statistik BSJ meloloskan 50 bacaan linear sebagai positif palsu karena
   penalti "dapat dijelaskan referensi linear" terlalu lemah untuk sekuens
   repetitif. Diubah menjadi diskualifikasi keras, yaitu kriteria keunikan gaya
   Szabo, dan jumlah bacaan terdiskualifikasi kini dilaporkan sebagai peringatan
   bahwa konstruknya terlalu repetitif.

Satu asumsi saya sendiri juga patah di tes: saya menuliskan bahwa MFE sirkular
selalu lebih rendah daripada linear. ViennaRNA mengembalikan +3,83 (sirkular)
berbanding −0,03 (linear) untuk satu sekuens uji, karena model sirkular
membebankan energi loop eksterior pada lingkaran tak berstruktur. Tesnya diganti
dan alasannya ditulis di docstring.

Keluaran: 109 tes lulus, ruff bersih, Snakemake 8 aturan, CI tiga job.

## 2026-11-16 — Jalankan demo dengan sekuens placeholder

Atas permintaan "fake dulu aja", seluruh alur dijalankan pada
`config/parts.example.yaml`, yang setiap bagiannya ditandai `synthetic: true`
dan dibangkitkan deterministik oleh `data/example/make_example.py`. Tidak satu
nukleotida pun berasal dari CVB3, intron Anabaena, atau genom HBV.

Hasil: prekursor 1867 nt, lingkaran 1407 nt, BSJ pada 1407, ORF CAI 0,9087 dan
GC 0,5492, skor komposit **0,6200 (GO-BERSYARAT)**, tanpa gerbang keras yang
gagal. Simulasi BSJ 102/102 lolos dengan presisi 1,0 dan nol positif palsu dari
bacaan linear.

Label laporan diperkuat menjadi `[DEMO — SEKUENS PLACEHOLDER]` yang menyatakan
secara eksplisit bahwa angkanya tidak boleh dikutip sebagai temuan, dan CI
diberi asersi yang memastikan label itu benar-benar muncul.

## 2026-11-17 — Siklus optimasi menutup loop

`circhbsag optimize` dijalankan pada konstruk demo: 24 trial Optuna (NSGA-II, 6
objektif) menghasilkan 16 titik di front Pareto. Setelah front diperingkat ulang
dengan skor komposit Tahap 0, konstruk terbaik mencapai **0,7228 (GO)**, naik
dari 0,6200, dengan spacer 25/50 nt pustaka polyAC tanpa lengan homologi.

Jadi rantai nilai → optimasi → peringkat ulang → nilai lagi terbukti tertutup.
Nilai optimum spacer itu khusus untuk sekuens placeholder dan harus dihitung
ulang pada sekuens sungguhan.

## 2026-11-17 — Rujukan diverifikasi lewat PubMed

Sepuluh sitasi ditarik dari PubMed beserta PMID dan DOI dari metadatanya, bukan
dari ingatan, dan disimpan di `content/references.yaml`. Satu di antaranya
mengubah argumen relevansi lokal proyek: Thedja dkk. 2011 (Arch Virol
156:855–868) memetakan 419 individu dari 40 populasi Indonesia dan menemukan
HBV/B dominan pada populasi penutur Austronesia dengan **subgenotipe B3 dominan
di paruh barat kepulauan**, sementara HBV/C dominan di Papua dan Maluku yang
terpengaruh Papua. Thedja dkk. 2015 menambahkan bahwa HBV/C terpecah menjadi
tipe Asia dan Papua-Pasifik yang berbeda justru di epitop core dan surface.

Konsekuensinya ditulis apa adanya di `content/decisions.yaml`: satu konstruk
tunggal tidak menutup seluruh kepulauan, dan itu harus dinyatakan sebagai batas
cakupan, bukan dibiarkan terbaca sebagai klaim nasional.

## 2026-11-17 — Gate 0 sembilan uji

Modul `gate0.py` menjalankan sembilan uji berurutan, masing-masing dengan angka,
ambang, dasar ambang, dan verdict. Pada konstruk demo: **6 PASS, 2 REVISE, 1
belum bisa dijalankan**, sehingga keputusannya ditulis `BELUM-LENGKAP`, bukan
PASS.

- Uji 2 (integritas ribozim grup I terhadap Rfam RF00028 dengan `cmsearch`)
  berstatus `UNAVAILABLE`: Infernal belum terpasang dan Rfam tidak dapat
  dijangkau dari sesi cloud ini. Uji ini mengembalikan perintah persisnya, bukan
  angka tebakan, dan tidak dihitung sebagai lulus.
- Uji 3 memberi REVISE: lengan homologi 40+40 nt pada data demo menghasilkan ΔG
  pengikatan −45,08 kcal/mol pada 55 °C (RNAcofold), lebih kuat daripada ambang
  −8 sampai −25 kcal/mol yang saya tetapkan sendiri. Ambang itu berlabel
  `[PERKIRAAN]` dan `[PERLU VERIFIKASI]`; lengan demo memang komplemen sempurna
  40 nt, jadi hasil ini wajar untuk data buatan.
- Uji 9 memberi REVISE: tujuh jendela GC 50 nt di luar rentang 35–70% dan dua
  poly-U 5 nt.

Lima dari sembilan ambang masih berlabel `[PERKIRAAN]`. Itu dinyatakan di
laporan, bukan disamarkan.

Keluaran: `results/demo/gate0.json`, `gate0.md`, dan `design_lock.json`.

## 2026-11-17 — Sumber tunggal dan generator

`content/` menjadi satu-satunya tempat identitas proyek, daftar pustaka,
keputusan desain, tabel software, Assumption Register, dan Open Questions
disimpan. `tools/build.py` menghasilkan situs multi-halaman, dokumen `.docx`,
dan naskah dek dari sumber itu ditambah `design_lock.json`.

Validator menolak konten yang menyitasi kunci referensi yang tidak ada, memakai
label provenans yang tidak dideklarasikan, atau mencantumkan alat tanpa lisensi.
`python -m tools.build check` membangun ulang situs ke direktori sementara dan
gagal bila yang di repositori berbeda, sehingga ketidakcocokan antar-deliverable
menjadi berkas kedaluwarsa yang terdeteksi CI, bukan angka yang tidak bisa
dilacak.

## 2026-11-18 — Optimasi diulang setelah Gate 0 diperbaiki

Angka 0,7228 pada entri 2026-11-17 dihitung sebelum batas GC/poly-U Gate 0
disatukan dengan optimasi kodon dan sebelum arah pemangkasan lengan diperbaiki,
jadi tidak lagi berlaku. Setelah dua perbaikan itu dan loop redesain (lengan 20 nt,
GC spacer 0,55; `results/demo/gate0_iterations.md`), skor komposit konstruk dasar
0,6922 (GO-BERSYARAT). `circhbsag optimize` diulang: 40 trial, 24 titik front
Pareto; setelah peringkat ulang konstruk terbaik 0,7437 (GO), spacer 50/30 nt
polyAC, lengan homologi 9 nt. Angka ini khusus sekuens placeholder (`[HASIL IN
SILICO]`, bukan data eksperimen). Gate 0 pada data demo tetap FAIL-PINDAH-JALUR
karena Uji 2 (intron placeholder), terlepas dari skor komposit.
