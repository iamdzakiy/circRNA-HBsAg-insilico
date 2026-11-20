# Rencana optimasi in silico vaksin circRNA-HBsAg

Dokumen ini adalah rencana kerja, bukan laporan hasil. Ia menjawab satu
pertanyaan: **bagaimana setiap komponen konstruk circRNA-HBsAg dioptimasi
secara in silico, di macOS, tanpa biaya, berjalan otomatis, dengan memanfaatkan
data berskala besar — dan bagaimana kita tahu optimasinya berhasil?**

Yang sudah berjalan hari ini ada di [README](../README.md) dan
[LOGBOOK](LOGBOOK.md). Yang belum, tertulis di sini sebagai rencana dengan
kriteria keberhasilan yang bisa diuji.

## Konvensi label

Setiap pernyataan di dokumen ini membawa satu label:

| Label | Arti |
|---|---|
| `[LITERATUR]` | diambil dari publikasi yang sitasinya ada di [Daftar pustaka](#daftar-pustaka) |
| `[PERKIRAAN]` | penalaran, asumsi, atau angka kasar saya sendiri; belum diverifikasi |
| `[HASIL IN SILICO]` | dihitung oleh pipeline ini pada sekuens sungguhan |
| `[DEMO]` | dihitung oleh pipeline ini pada sekuens placeholder; hanya ilustrasi alur |
| `[RENCANA]` | belum dikerjakan |

Tidak ada angka di dokumen ini yang dikarang. Bila sebuah angka tidak punya
sumber dan belum dihitung, ia ditulis sebagai pertanyaan terbuka, bukan sebagai
nilai.

---

## 1. Lima prinsip yang mengikat seluruh rencana

1. **Satu keputusan, satu fungsi tujuan.** Optimasi tanpa fungsi tujuan yang
   ditulis eksplisit tidak bisa diaudit. Setiap komponen di Bab 3 punya tabel
   variabel keputusan, fungsi tujuan, dan kriteria sukses.
2. **Keputusan dua lapis.** Gerbang keras (hard gate) memveto; kriteria
   berbobot memberi peringkat. Konstruk yang menabrak satu gerbang keras tidak
   diselamatkan oleh skor tinggi di kriteria lain.
3. **Prediksi bukan pengukuran.** Semua yang dihasilkan tahap ini adalah
   prediksi termodinamika dan statistika. Nilainya terletak pada **urutan**
   antar-varian desain, bukan pada nilai absolutnya.
4. **Provenans ikut ke hilir.** Setiap laporan membawa label sumber masukannya.
   Pipeline menolak jalan bila sekuens referensi belum ada, bukan melanjutkan
   dengan pengganti.
5. **Gratis dan reprodusibel di laptop.** Tidak ada komponen yang mewajibkan
   lisensi berbayar, GPU, atau klaster. Bila sebuah alat memberi hasil lebih
   baik tetapi lisensinya membatasi, alat itu dicatat dan **tidak** dipakai.

---

## 2. Lingkungan macOS, gratis, dan otomatis

### 2.1 Basis lingkungan

`[PERKIRAAN]` Rekomendasi: **Miniforge (conda-forge)** dengan `mamba`, bukan
Anaconda Distribution, karena saluran `defaults` milik Anaconda memuat syarat
lisensi komersial sedangkan conda-forge dan bioconda tidak. Di Apple Silicon,
jalankan dua kanal arsitektur:

```bash
# Apple Silicon (M1-M4): lingkungan utama native arm64
CONDA_SUBDIR=osx-arm64 mamba env create -f environment.yml
conda activate circhbsag
conda config --env --set subdir osx-arm64

# Bila satu paket belum punya build osx-arm64, buat lingkungan kedua di
# bawah Rosetta 2 dan panggil hanya paket itu dari Snakemake:
CONDA_SUBDIR=osx-64 mamba create -n circhbsag-x64 -c bioconda -c conda-forge <paket>
```

`[PERKIRAAN]` Paket yang ketersediaan `osx-arm64`-nya perlu diperiksa saat
pemasangan, bukan diasumsikan: STAR, CIRI2, dan beberapa alat NGS lama lebih
sering hanya tersedia `osx-64`. Cara memeriksa, bukan menduga:

```bash
mamba repoquery search 'star' --channel bioconda --platform osx-arm64
```

Snakemake memegang satu `conda:` per-rule, sehingga aturan yang butuh `osx-64`
bisa hidup berdampingan dengan aturan native tanpa mengubah lingkungan utama.

### 2.2 Tabel alat: lisensi dan ketersediaan macOS

`[LITERATUR]` untuk kolom lisensi (teks lisensi masing-masing proyek);
`[PERKIRAAN]` untuk kolom macOS, yang harus dikonfirmasi saat pemasangan.

| Alat | Fungsi | Lisensi | macOS |
|---|---|---|---|
| ViennaRNA | lipatan sirkular, MFE, aksesibilitas | ViennaRNA License (gratis untuk akademik, non-komersial) | arm64 + x64 |
| LinearFold | pembanding lipatan | akademik gratis, komersial perlu izin | sumber, kompilasi sendiri |
| MXfold2 | pembanding lipatan berbasis ML | MIT | pip, arm64 |
| primer3 / primer3-py | desain primer | GPLv2 / GPLv2 | arm64 |
| Biopython | I/O sekuens, utilitas | Biopython License (mirip MIT) | arm64 |
| NumPy, pandas, SciPy | numerik | BSD-3-Clause | arm64 |
| Optuna | optimasi multi-objektif (NSGA-II) | MIT | arm64 |
| Snakemake | orkestrasi alur | MIT | arm64 |
| BLAST+ | homologi off-target | public domain (karya pemerintah AS) | x64 (Rosetta) |
| ART / wgsim | simulasi bacaan | GPL / MIT | x64 |
| STAR | penyelarasan kimerik | MIT | periksa; mungkin x64 |
| CIRI2 | pemanggil circRNA | GPL | Perl, lintas platform |
| Streamlit | antarmuka web | Apache-2.0 | arm64 |
| DeepCIP | prediksi IRES circRNA | periksa repositori sebelum dipakai | PyTorch, arm64 |
| IEDB tools (BepiPred, NetMHCpan, NetMHCIIpan) | epitop | gratis untuk akademik, **wajib unduh dengan akun** dan dilarang redistribusi | x64/arm64 |
| **Tidak dipakai: LinearDesign** | desain mRNA optimal | membatasi penggunaan komersial | — |
| **Tidak dipakai: NUPACK** | termodinamika asam nukleat | akademik, melarang redistribusi | — |

Dua baris terakhir penting untuk dokumen proyek: `[LITERATUR]` LinearDesign
melaporkan kenaikan titer antibodi hingga 128 kali dibanding pembanding
optimasi kodon pada vaksin mRNA (Zhang et al. 2023), jadi menolaknya adalah
pengorbanan nyata, bukan kelalaian. Alasannya lisensi, dan harus ditulis
sebagai keterbatasan, bukan disembunyikan.

### 2.3 Otomatisasi penuh

`[HASIL IN SILICO]` Sudah berjalan: Snakemake 8 aturan dari parts → keputusan
Tahap 0 → primer → simulasi BSJ → off-target → optimasi → ringkasan, dengan
seed tetap, dan GitHub Actions tiga job (tes, alur ujung-ke-ujung, penjaga
provenans).

`[RENCANA]` Yang ditambahkan oleh rencana ini:

1. `make all` satu perintah di macOS yang membuat lingkungan, menjalankan
   alur, lalu membangun situs.
2. Job GitHub Actions `pages` yang membangun dan menerbitkan situs progres
   setiap kali `main` berubah.
3. Jadwal `cron` mingguan yang menjalankan ulang seluruh alur pada sekuens
   yang sama dan **gagal bila angkanya bergeser** di luar toleransi — uji
   reprodusibilitas yang berjalan sendiri.
4. Satu berkas `logbook/entries/*.md` per jalankan, ditulis oleh pipeline,
   sehingga logbook tumbuh otomatis dan bukan hasil ingatan.

`[PERKIRAAN]` Biaya: nol untuk repositori publik, karena GitHub Actions dan
GitHub Pages gratis tanpa batas menit untuk repositori publik, sementara
repositori privat punya kuota menit bulanan dan Pages privat butuh paket
berbayar. Konsekuensinya langsung ke permintaan "website yang bisa diakses
semua": **situs publik praktis mengharuskan repositorinya publik.** Batas kuota
yang berlaku harus dibaca di halaman penagihan GitHub saat repositori dibuat,
bukan dikutip dari dokumen ini.

---

## 3. Optimasi per komponen

Konstruk, dari hulu ke hilir:

```
T7 | GGG | HA5 | intron 3' PIE | E2 | spacer5 | IRES | ORF HBsAg | spacer3 | E1 | intron 5' PIE | HA3
                                 └──────────────── lingkaran setelah splicing ────────────────┘
```

### 3.1 Antigen: genotipe dan bentuk HBsAg

| | |
|---|---|
| Variabel keputusan | genotipe/subgenotipe; S saja vs L/M/S; mutasi penghindar vaksin yang perlu dihindari |
| Data skala besar | korpus genom HBV lengkap (HBVdb, NCBI Virus), puluhan ribu sekuens lintas genotipe `[LITERATUR]` |
| Metode | penyelarasan berganda per subgenotipe → sekuens konsensus → peta konservasi per residu → pemeriksaan determinan "a" (residu ~124–147) utuh |
| Kriteria sukses | `[RENCANA]` konsensus menutup ≥95% posisi dengan dukungan ≥90% sekuens; determinan "a" identik dengan konsensus subgenotipe sasaran |

`[LITERATUR]` Dasar pemilihan genotipe untuk Indonesia kini punya sumber yang
bisa disitasi: Thedja et al. (2011) memetakan 419 individu dari 40 populasi dan
menemukan HBV/B dominan pada populasi penutur Austronesia, dengan **subgenotipe
B3 dominan di paruh barat kepulauan** (penutur Malayo-Polinesia Barat),
sementara **HBV/C dominan di Papua dan populasi Maluku yang terpengaruh Papua**;
subgenotipe B7, B8, dan B9 spesifik Nusa Tenggara. Thedja et al. (2015)
menambahkan bahwa HBV/C di kepulauan ini terpecah menjadi dua tipe, Asia dan
Papua-Pasifik, yang berbeda pada beberapa substitusi asam amino **di epitop imun
protein core dan surface** — jadi pilihan genotipe bukan soal label, tetapi
langsung mengubah epitop yang divaksinkan.

`[PERKIRAAN]` Implikasi desain: untuk sasaran "Indonesia" secara umum, satu
konstruk tunggal tidak menutup seluruh kepulauan. Pilihan yang jujur ada tiga,
dan ini keputusan penulis proyek, bukan keputusan pipeline:

1. **B3** sebagai sasaran utama, dengan cakupan dinyatakan terbatas pada paruh
   barat yang berpenduduk paling banyak;
2. **konstruk konsensus B3+C** yang dioptimasi untuk epitop bersama, dengan
   konsekuensi kecocokan yang lebih rendah pada keduanya;
3. **dua konstruk** sebagai kandidat paralel, yang menggandakan beban in vitro.

### 3.2 ORF: optimasi kodon

| | |
|---|---|
| Variabel keputusan | 1 kodon per residu; untuk 226 aa ruangnya `[PERKIRAAN]` berorde 10^80 |
| Fungsi tujuan | maksimalkan CAI dan stabilitas struktur lokal; minimalkan GC ekstrem, homopolimer, situs Tipe IIS, motif terminator T7, situs m6A tak diinginkan |
| Data skala besar | tabel penggunaan kodon inang dari korpus CDS manusia; `[RENCANA]` korpus efisiensi translasi untuk validasi silang |
| Metode | `[HASIL IN SILICO]` sudah: pendakian bukit berseed tetap. `[RENCANA]` tambah: pencarian berkas (beam search) atas lintasan kodon, pembanding wajib |
| Kriteria sukses | CAI varian terpilih ≥ CAI pendakian bukit, tanpa menambah pelanggaran motif; MFE jendela di sekitar kodon awal tidak memburuk |

`[LITERATUR]` Zhang et al. (2023) menunjukkan struktur sekunder yang lebih
panjang memperpanjang waktu paruh mRNA dan, bersama kodon optimal, memperbaiki
ekspresi protein — jadi fungsi tujuan yang hanya memaksimalkan CAI memang
kurang. Itu yang membuat MFE masuk ke fungsi tujuan di atas.

### 3.3 IRES

| | |
|---|---|
| Variabel keputusan | sudah ditetapkan satu IRES (CVB3); yang tersisa: titik potong 5'/3' dan apakah domain terjaga utuh dalam lipatan lingkaran |
| Data skala besar | `[RENCANA]` korpus IRES untuk penyaringan in silico |
| Metode | `[HASIL IN SILICO]` sudah: kriteria `ires_integrity`, memeriksa domain IRES tidak dirusak pasangan basa dari luar. `[RENCANA]` tambah: skor DeepCIP sebagai pembanding independen, dengan label eksperimental |
| Kriteria sukses | `ires_integrity` = 1.0 pada suhu 37 °C dan 55 °C; peringkat varian dari DeepCIP tidak berlawanan arah dengan peringkat termodinamika |

`[LITERATUR]` DeepCIP (Zhou et al. 2023) adalah prediktor IRES yang dilatih
khusus untuk circRNA, bukan untuk RNA linear, dan menggabungkan informasi
sekuens dengan struktur. `[LITERATUR]` LAMAR (Zhou et al. 2025), model bahasa
RNA yang dilatih pada ~15 juta sekuens dari 225 mamalia dan 1569 virus,
melaporkan penyaringan in silico IRES baru yang menghasilkan IRES aktif tinggi
pendorong translasi circRNA — inilah contoh paling dekat dengan "big data untuk
mengoptimasi komponen vaksin circRNA" yang menjadi tujuan proyek. `[PERKIRAAN]` Keduanya
masuk sebagai **pembanding berlabel eksperimental**, bukan pengganti kriteria
termodinamika, karena keduanya belum pernah divalidasi pada konstruk PIE-HBsAg.

### 3.4 Spacer dan lengan homologi

| | |
|---|---|
| Variabel keputusan | panjang spacer 5' dan 3', pustaka komposisi, GC sasaran, panjang lengan homologi |
| Fungsi tujuan | maksimalkan aksesibilitas sambungan balik; minimalkan dsRNA panjang; jaga aksesibilitas kodon awal |
| Metode | `[HASIL IN SILICO]` NSGA-II Optuna 6 objektif, front Pareto diperingkat ulang dengan skor komposit Tahap 0 |
| Kriteria sukses | skor komposit varian terbaik melewati ambang GO 0.70 |

`[LITERATUR]` Wesselhoeft et al. (2018) menunjukkan lengan homologi dan spacer
bukan hiasan: efisiensi sirkularisasi naik dari 0% tanpa lengan menjadi 16%
dengan lengan 9 nt dan 48% dengan 19 nt, spacer permisif menaikkan 46% → 87%,
dan satu spacer yang mengganggu menghapus splicing sama sekali; intron Anabaena
unggul dibanding alternatif yang diuji. `[DEMO]` Pada konstruk placeholder,
siklus optimasi memindahkan skor dari 0.6922 (GO-BERSYARAT) ke 0.7437 (GO)
dengan spacer 50/30 nt pustaka polyAC dan lengan homologi 9 nt — angka ini khusus
untuk sekuens placeholder dan tidak boleh digeneralkan.

### 3.5 Intron PIE dan IGS

| | |
|---|---|
| Variabel keputusan | pasangan P1/P10 dan sekuens IGS; keduanya harus komplementer dengan ujung ekson |
| Metode | `[HASIL IN SILICO]` RNAplfold pada jendela P1 dan P10 untuk memastikan keduanya tidak terkubur struktur |
| Kriteria sukses | probabilitas tak-berpasangan rata-rata jendela P1 dan P10 di atas ambang yang ditetapkan di `config/` |

### 3.6 m6A

`[HASIL IN SILICO]` Pemindaian RRACH sudah berjalan dan masuk sebagai beban
motif. `[RENCANA]` Keputusan apakah menyisipkan m6A tetap **tidak menyisipkan
secara baku**, karena arah efeknya pada konstruk ini belum ditentukan oleh
literatur yang saya verifikasi. `[LITERATUR]` Bahwa translasi circRNA dapat
berjalan lewat mekanisme bergantung m6A **maupun** lewat IRES dicatat oleh Zhu
et al. (2025); karena konstruk ini sudah memakai IRES, menambah m6A berarti
menambah variabel tanpa hipotesis yang jelas.

### 3.7 Epitop

| | |
|---|---|
| Variabel keputusan | tidak ada; ORF ditentukan oleh antigen. Yang dinilai: apakah epitop B dan T tetap utuh setelah optimasi kodon dan pemilihan konsensus |
| Data skala besar | `[LITERATUR]` IEDB sebagai repositori epitop terkurasi; `[RENCANA]` frekuensi alel HLA populasi Indonesia untuk pembobotan cakupan |
| Metode | `[RENCANA]` prediksi pengikatan dengan NetMHCpan-4.1 dan NetMHCIIpan-4.0, lalu cakupan populasi |
| Kriteria sukses | `[RENCANA]` cakupan alel kelas I dan II di atas ambang yang ditetapkan sebelum perhitungan, bukan setelah melihat hasil |

`[LITERATUR]` NetMHCpan-4.1 dan NetMHCIIpan-4.0 (Reynisson et al. 2020)
menggabungkan data afinitas dan ligan tereluasi spektrometri massa. `[PERKIRAAN]`
Keduanya gratis untuk akademik tetapi **tidak boleh diredistribusi**, sehingga
tidak bisa dipasang oleh CI publik; langkah ini harus berjalan lokal di macOS
dan hanya hasilnya yang masuk repositori, dengan versi alat dicatat.

### 3.8 Kemurnian: RNA linear sisa

`[LITERATUR]` Wesselhoeft et al. (2019) menemukan imunogenisitas dan stabilitas
ekspresi circRNA bergantung pada kemurnian: sedikit saja RNA linear yang
mengontaminasi memicu respons imun seluler yang kuat. `[HASIL IN SILICO]`
Konsekuensi in silico-nya sudah dipasang: filter statistik BSJ gaya Szabo
mendiskualifikasi bacaan yang dapat dijelaskan oleh referensi linear, dan
jumlah bacaan terdiskualifikasi dilaporkan sebagai peringatan bahwa konstruk
terlalu repetitif. `[RENCANA]` Yang perlu ditambah: estimasi batas deteksi RNA
linear dari simulasi, sebagai angka target untuk tahap in vitro.

---

## 4. Lapis data skala besar

Empat korpus, dengan peran yang jelas. `[PERKIRAAN]` untuk ukuran; semuanya
harus dibaca ulang dari sumber saat diunduh, karena nomornya berubah.

| Korpus | Dipakai untuk | Catatan akses |
|---|---|---|
| Genom HBV lintas genotipe (HBVdb, NCBI Virus) | konsensus antigen, peta konservasi, deteksi mutasi penghindar vaksin | `[LITERATUR]` HBVdb memuat sekuens klinis seluruh genotipe A–H; satu studi melaporkan analisis hampir 70.000 sekuens darinya (Meier-Stephenson et al. 2018) |
| CDS manusia | tabel kodon inang yang dihitung sendiri, bukan diwarisi | pipeline sudah punya `circhbsag codon-table` untuk ini |
| Korpus IRES / model bahasa RNA | penyaringan dan pembanding IRES | `[LITERATUR]` LAMAR dilatih pada ~15 juta sekuens (Zhou et al. 2025) |
| IEDB + frekuensi alel HLA | cakupan epitop berbobot populasi | unduhan besar; jalankan lokal |

`[PERKIRAAN]` Semua korpus di atas berukuran ratusan MB sampai beberapa GB,
yang masih nyaman di laptop. Yang tidak nyaman di laptop adalah pelatihan ulang
model bahasa RNA; karena itu rencana ini hanya memakai model terlatih untuk
inferensi, tidak melatih apa pun.

**Catatan:** unduhan korpus di bab ini belum dijalankan; semuanya berstatus `[RENCANA]` untuk dijalankan di macOS.
Pipeline sengaja dibuat berhenti dengan galat bila berkasnya belum ada.

---

## 5. Mesin optimasi

Lima tahap, berurutan, semuanya sudah ada kodenya kecuali yang bertanda
`[RENCANA]`:

1. **Enumerasi terbatas** varian spacer dan lengan homologi dari pustaka
   literatur. Murah, dan menjadi garis dasar yang harus dikalahkan.
2. **NSGA-II** (Optuna, MIT) atas 6 objektif: aksesibilitas sambungan,
   keutuhan IRES, beban dsRNA, aksesibilitas kodon awal, stabilitas lingkaran,
   kebersihan sekuens.
3. **Peringkat ulang** front Pareto dengan skor komposit Tahap 0, supaya
   keputusan akhir memakai satu fungsi yang bobotnya tertulis.
4. **Uji ketahanan**: dua suhu (37 °C dan 55 °C, suhu reaksi sirkularisasi),
   dan jitter bobot gaya Dirichlet untuk melihat apakah peringkat bertahan bila
   bobotnya digeser. Keputusan yang berubah arah karena bobot digeser sedikit
   bukan keputusan.
5. **Validasi** `[RENCANA]`: homologi off-target terhadap genom manusia dengan
   BLAST+, dan simulasi bacaan → pemanggilan BSJ untuk memastikan konstruk
   terdeteksi dan tidak memalsukan deteksi.

---

## 6. Tonggak dan gerbang

| Tonggak | Keluaran | Gerbang |
|---|---|---|
| M1 | korpus terunduh, `parts.yaml` lengkap dengan aksesi | pipeline jalan tanpa galat provenans |
| M2 | konsensus antigen + peta konservasi | determinan "a" utuh |
| M3 | Tahap 0 pada sekuens sungguhan | tidak ada gerbang keras gagal |
| M4 | siklus optimasi | skor komposit ≥ 0.70 |
| M5 | epitop + cakupan populasi | ambang cakupan yang ditetapkan di muka |
| M6 | berkas desain final + protokol sintesis | semua angka berlabel `[HASIL IN SILICO]` |

`[PERKIRAAN]` M3 sampai M6 masing-masing berorde jam di laptop, bukan hari,
karena konstruknya pendek (~1,4 kb). Yang memakan waktu adalah M1 dan M2, yaitu
kerja kurasi data, bukan komputasi.

---

## 7. Apa yang tidak boleh diklaim dari rencana ini

1. Tidak ada klaim imunogenisitas. Prediksi epitop bukan pengukuran respons imun.
2. Tidak ada klaim efisiensi sirkularisasi. Aksesibilitas sambungan berkorelasi
   dengan harapan sirkularisasi, dan itu bukan nilai efisiensi.
3. Model termodinamika yang dipakai tidak menangani pseudoknot, dan IRES
   pikornavirus memuat struktur tersier yang tidak diwakili model itu.
4. Prediktor berbasis pembelajaran mesin dipakai sebagai pembanding berlabel
   eksperimental. Bila ia bertentangan dengan termodinamika, yang dilaporkan
   adalah pertentangannya, bukan yang hasilnya lebih enak dibaca.
5. Seluruh angka pada data placeholder berlabel `[DEMO]` dan tidak boleh masuk
   dokumen proyek sebagai temuan.

---

## Daftar pustaka

Sitasi di bawah diambil dari PubMed pada 17 November 2026. Setiap DOI sudah
diverifikasi lewat metadata PubMed, bukan ditulis dari ingatan.

1. Wesselhoeft RA, Kowalski PS, Anderson DG. Engineering circular RNA for
   potent and stable translation in eukaryotic cells. *Nat Commun* 2018;
   9(1):2629. <https://doi.org/10.1038/s41467-018-05096-6>
2. Wesselhoeft RA, Kowalski PS, Parker-Hale FC, Huang Y, Bisaria N, Anderson DG.
   RNA circularization diminishes immunogenicity and can extend translation
   duration in vivo. *Mol Cell* 2019; 74(3):508–520.e4.
   <https://doi.org/10.1016/j.molcel.2019.02.015>
3. Zhang H, Zhang L, Lin A, et al. Algorithm for optimized mRNA design improves
   stability and immunogenicity. *Nature* 2023; 621(7978):396–403.
   <https://doi.org/10.1038/s41586-023-06127-z>
4. Thedja MD, Muljono DH, Nurainy N, Sukowati CHC, Verhoef J, Marzuki S.
   Ethnogeographical structure of hepatitis B virus genotype distribution in
   Indonesia and discovery of a new subgenotype, B9. *Arch Virol* 2011;
   156(5):855–868. <https://doi.org/10.1007/s00705-011-0926-y>
5. Thedja MD, Muljono DH, Ie SI, Sidarta E, Turyadi, Verhoef J, Marzuki S.
   Genogeography and immune epitope characteristics of hepatitis B virus
   genotype C reveals two distinct types: Asian and Papua-Pacific.
   *PLoS One* 2015; 10(7):e0132533.
   <https://doi.org/10.1371/journal.pone.0132533>
6. Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M. NetMHCpan-4.1 and
   NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by
   concurrent motif deconvolution and integration of MS MHC eluted ligand data.
   *Nucleic Acids Res* 2020; 48(W1):W449–W454.
   <https://doi.org/10.1093/nar/gkaa379>
7. Zhou Y, Wu J, Yao S, Xu Y, Zhao W, Tong Y, Zhou Z. DeepCIP: a multimodal deep
   learning method for the prediction of internal ribosome entry sites of
   circRNAs. *Comput Biol Med* 2023; 164:107288.
   <https://doi.org/10.1016/j.compbiomed.2023.107288>
8. Zhou H, Hu Y, Zheng Y, et al. A foundation language model to decipher diverse
   regulation of RNAs. *Genome Biol* 2025; 26(1):301.
   <https://doi.org/10.1186/s13059-025-03752-x>
9. Zhu Z, Liang X, Ma R, Yin S, Xu M, Li G. CircCode3: integrating deep learning
   to mine and evaluate translatable circular RNAs from ribosome profiling
   sequencing and mass spectrometry data. *Brief Bioinform* 2025; 26(5):bbaf458.
   <https://doi.org/10.1093/bib/bbaf458>
10. Meier-Stephenson V, Bremner WTR, Dalton CS, van Marle G, Coffin CS,
    Patel TR. Comprehensive analysis of hepatitis B virus promoter region
    mutations. *Viruses* 2018; 10(11):603.
    <https://doi.org/10.3390/v10110603>

Sumber metadata: PubMed (NCBI).
