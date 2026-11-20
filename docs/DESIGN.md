# Keputusan desain dan justifikasinya

Dokumen ini mencatat **mengapa** konstruk dirancang seperti ini, dengan
penandaan sumber yang konsisten:

- **[LITERATUR]**: dinyatakan pustaka yang dapat diperiksa (sitasi diberikan)
- **[PERKIRAAN]**: ditetapkan pipeline ini sebagai nilai kerja; dapat diubah
- **[HASIL IN SILICO]**: dihitung oleh pipeline ini saat dijalankan

---

## 1. Mengapa satu IRES, dan mengapa CVB3

Lima varian IRES sempat dipertimbangkan (CVB3, EMCV, poliovirus, HCV, dan satu
kandidat hasil prioritisasi model). Pipeline ini mempersempit pilihan ke
**satu IRES, yaitu CVB3**, dengan alasan berikut.

**[LITERATUR]** Wesselhoeft, Kowalski & Anderson (2018), *Engineering circular
RNA for potent and stable translation in eukaryotic cells*, Nature
Communications 9:2629:

- IRES CVB3 **1,5 kali lebih efektif** daripada IRES EMCV yang lazim dipakai, di
  sel HEK293.
- Ketika diuji pada konteks ORF yang berbeda (Firefly luciferase), CVB3 **tetap
  superior terhadap semua IRES lain yang diuji**, sementara efikasi beberapa
  IRES lain, terutama poliovirus, **berubah drastis**.
- Pada beberapa jenis sel (HeLa, A549, Min6), efikasi IRES berbeda-beda, tetapi
  CVB3 tetap yang terbaik.

Poin kedua adalah alasan metodologis terkuat: efikasi IRES selain CVB3
**bergantung pada ORF yang mengikutinya**. Pada proyek ini payload-nya baru
(ORF HBsAg), sehingga membandingkan lima IRES berarti membandingkan lima sistem
yang masing-masing butuh optimasi spacer sendiri agar adil, pekerjaan yang jauh
melampaui lingkup satu proyek kecil. Dengan memilih IRES yang paling kokoh
lintas konteks, variabel IRES dikeluarkan dari percobaan dan sumber daya
dialihkan ke pertanyaan yang benar-benar baru, yaitu apakah payload HBsAg
mempertahankan sirkularisasi dan translasi.

**Konsekuensinya:** proyek ini tidak membandingkan varian IRES secara sistematis.
Kebaruannya ada di payload dan di alur keputusan GO/NO-GO in silico.
Perbandingan IRES tetap terbuka sebagai pekerjaan lanjutan.

### Catatan penting tentang konteks prokariot

**[LITERATUR]** IRES tidak berfungsi sebagai pemicu translasi pada ribosom
prokariot. Bila validasi awal dilakukan di *E. coli*,
IRES harus diperlakukan sebagai **konteks sekuens yang memengaruhi splicing dan
stabilitas**, bukan sebagai elemen yang fungsi translasinya diukur. Pipeline ini
konsisten dengan itu: metrik `ires_integrity` adalah **proksi struktural**
(apakah lipatan IRES terjaga), bukan prediksi aktivitas translasi.

---

## 2. Mengapa ORF S-HBsAg kecil, bukan L/M/S atau antigen multi-epitop

Tiga pilihan payload dipertimbangkan.

| Pilihan | Panjang ORF | Pertimbangan |
|---|---|---|
| **S kecil (SHBs, ~226 aa)** | ~681 nt | **dipilih** |
| M dan L (pre-S2/pre-S1) | ~846 / ~1.200 nt | menambah 0,2–0,5 kb; epitop pre-S bernilai, tetapi menyulitkan perakitan partikel |
| MEV multi-epitop (≤ 200 aa) | ~600 nt | bergantung penuh pada prediksi epitop in silico yang belum tervalidasi |

Alasan memilih S kecil:

1. **Antigen yang sama dengan vaksin hepatitis B berlisensi.** HBsAg adalah
   antigen pada vaksin HB rekombinan yang dipakai sejak 1980-an, sehingga
   korelasi antara antibodi anti-HBs dan perlindungan sudah mapan dan hasil uji
   dapat dibandingkan dengan standar yang ada. MEV multi-epitop tidak punya
   pembanding semacam itu.
2. **Protein S berukuran kecil mampu merakit diri menjadi partikel subviral**
   tanpa protein M/L, sehingga satu ORF saja sudah cukup untuk antigen
   berbentuk partikel, bentuk yang lebih imunogenik daripada protein terlarut.
3. **Anggaran panjang konstruk.** Pipeline membatasi konstruk ≤ 3,0 kb. Dengan
   IRES CVB3 (~0,6 kb), dua separuh intron + eksón (~0,5 kb), spacer dan lengan
   homologi (~0,1 kb), ORF 0,68 kb menyisakan margin; ORF L/M/S akan menghabiskan
   margin itu. Pipeline menegakkan batas ini sebagai gerbang keras
   `circle_length_within_limit`.
4. **Satu variabel baru saja.** Payload sudah merupakan variabel baru; menambah
   ketidakpastian prediksi epitop akan membuat kegagalan sulit ditafsirkan.

**Genotipe: B** (dengan C sebagai alternatif yang didokumentasikan).
Genotipe B dan C disebut dominan di Indonesia, sehingga sekuens HBsAg diambil dari
genotipe tersebut dan genotipe B dipilih sebagai desain utama. Klaim dominansi
ini belum diverifikasi pipeline, jadi periksa dengan sitasi epidemiologi
Indonesia yang dapat diperiksa, dan catat nomor aksesi sekuens yang dipakai di
`config/parts.yaml`.

**Catatan antigenik:** determinan "a" HBsAg (sekitar residu 124–147) adalah
epitop penetral utama dan terkonservasi antar-genotipe, sehingga pilihan
genotipe lebih memengaruhi kecocokan epitop sekunder daripada perlindungan
dasar. Nyatakan ini sebagai pertimbangan, bukan sebagai hasil.

---

## 3. Mengapa intron grup I Anabaena

**[LITERATUR]** Wesselhoeft dkk. (2018) menguji intron grup I alternatif dari
pre-tRNA *Anabaena* dan menyimpulkan intron tersebut **superior untuk
sirkularisasi** dibandingkan intron fag T4 yang dipermutasi: separuh intron
Anabaena panjangnya hampir sama dan **lebih kecil kemungkinannya tetap
berasosiasi setelah splicing** meski lengan homologinya lebih kuat. Pemisahan
intron setelah splicing penting karena intron yang masih menempel menjadi
kontaminan yang harus disingkirkan pada pemurnian dan menaikkan imunogenisitas
intrinsik.

---

## 4. Mengapa spacer dan lengan homologi dioptimasi, bukan ditetapkan

Angka-angka berikut adalah dasar mengapa dua elemen ini jadi variabel utama
pada optimasi multi-objektif.

**[LITERATUR]** Wesselhoeft dkk. (2018):

| Intervensi | Efisiensi splicing |
|---|---|
| tanpa lengan homologi | 0% |
| lengan homologi 9 nt ("weak") | 16% |
| lengan homologi 19 nt ("strong") | 48% |
| + spacer permisif (P1, P2) | 46% → **87%** |
| + spacer yang merusak struktur intron | **splicing hilang sama sekali** |

Dua hal yang dibaca dari tabel ini:

1. **Spacer bukan pengisi netral.** Spacer yang sama panjangnya bisa menaikkan
   splicing ke 87% atau menghilangkannya, tergantung apakah ia membiarkan
   struktur kedua separuh intron melipat independen. Inilah sebabnya pipeline
   menilai **aksesibilitas situs splice dan IGS (P1/P10)** dan bukan hanya
   panjang spacer.
2. **Hipotesis mekanistik yang dipakai:** situs splice PIE 3' berdekatan dengan
   IRES, dan karena keduanya sangat berstruktur, sekuens IRES dapat mengganggu
   pelipatan ribozim splicing, baik secara proksimal di situs splice 3' maupun
   secara jauh di situs splice 5' lewat kontak jarak-jauh. Karena itu kriteria
   `ires_integrity` dan `junction_accessibility` diukur pada **konstruk penuh**,
   bukan pada bagian-bagian terpisah.

**[LITERATUR]** Studi yang sama menemukan spacer **polyA** maupun **polyAC**
tak-berstruktur sama-sama menaikkan ekspresi protein pada semua konstruk, diduga
karena memberi pemisahan tak-berstruktur yang lebih lega di sekitar IRES; pada
IRES EMCV dan poliovirus, polyA memberi perbaikan melebihi polyAC, yang mungkin
menandakan peran pengikatan protein pengikat poliadenilat. Karena itu
`SPACER_LIBRARIES` pada pipeline berisi tiga varian tak-berstruktur
(`unstructured`, `polyA`, `polyAC`) dan `HOMOLOGY_ARM_LENGTHS` memuat 9 dan 19
nt dari studi tersebut, ditambah 0 (kontrol), 30, dan 40 nt sebagai ekstrapolasi
**[PERKIRAAN]**.

---

## 5. Mengapa dua suhu: 37 dan 55 °C

**[LITERATUR]** Reaksi splicing PIE in vitro lazim dijalankan pada sekitar
55 °C dengan GTP berlebih dan Mg²⁺, sedangkan circRNA hasilnya harus berfungsi
pada 37 °C di dalam sel. Keduanya menuntut sifat lipatan yang berbeda:

- pada **55 °C** yang penting adalah situs splice dan P1/P10 **terbuka** agar
  ribozim dapat bekerja;
- pada **37 °C** yang penting adalah IRES **terlipat benar** dan start codon
  tetap dapat dicapai, sementara dsRNA panjang tidak terbentuk.

Pipeline mengambil **nilai terburuk dari kedua suhu** untuk setiap kriteria,
sehingga konstruk yang hanya baik pada salah satu suhu tidak lolos.

---

## 6. Mengapa LinearDesign tidak dipakai

Pendekatan ala LinearDesign menarik untuk optimasi kodon, tetapi lisensi
LinearDesign membatasi pemakaian komersial dan redistribusi serta
memerlukan persetujuan terpisah, sehingga ia **tidak dijadikan dependensi**
agar repositori ini tetap dapat dirilis seluruhnya di bawah MIT. Penggantinya
adalah optimizer kodon sendiri di `src/circhbsag/codon.py`: hill-climbing
stokastik dengan seed tetap atas objektif gabungan CAI, GC total, GC jendela,
situs Tipe IIS, motif yang dihindari, dan kerapatan m6A RRACH.

Perbedaan yang perlu dinyatakan: LinearDesign mengoptimasi **MFE mRNA dan CAI
secara bersama lewat pemrograman dinamis pada kisi lipatan**, sehingga ia dapat
mencari struktur mRNA yang paling stabil. Optimizer di sini **tidak** mengoptimasi
MFE ORF secara langsung; struktur ORF ditangani di tingkat konstruk melalui
kriteria `start_codon_accessibility` dan `dsrna_burden` pada Tahap 0. Bila MFE
ORF ingin dioptimasi langsung, tambahkan objektif tersebut pada
`CodonObjective`, perhitungan ViennaRNA-nya sudah tersedia.

Hal serupa berlaku untuk **NUPACK** (gratis hanya untuk akademik dengan
pendaftaran, tidak boleh diredistribusi): tidak dijadikan dependensi wajib.
Pembanding lipatan alternatif yang lisensinya aman untuk dijadikan dependensi
opsional adalah **LinearFold** (MIT) dan **MXfold2** (MIT); keduanya belum
diintegrasikan dan ditandai sebagai pekerjaan terbuka di bawah.

---

## 7. Keputusan tentang m6A

Soal penyisipan m6A (misalnya 1%), pipeline **memindai** semua situs
RRACH (konsensus m6A) dan melaporkan kerapatannya, lalu secara bawaan
**menurunkan** jumlahnya lewat optimasi kodon.

**Dasar keputusan bawaan "jangan sisipkan m6A":** penyisipan m6A pada circRNA
sintetis adalah intervensi yang efeknya bergantung konteks dan belum ada
konsensus bahwa ia menguntungkan untuk vaksin circRNA; menambah satu variabel
kimia yang tidak terkendali ke konstruk yang payload-nya sudah baru akan membuat
kegagalan sulit ditafsirkan. Karena itu desain utama memakai
`m6a_mode: minimize`, dan penyisipan m6A diletakkan sebagai **lengan percobaan
terpisah**, bukan bagian desain dasar.

Bila Anda memutuskan sebaliknya, ubah `--m6a-mode ignore` dan **bobot ulang**
kriteria `motif_burden`, jangan biarkan pipeline meminimalkan sesuatu yang
sengaja Anda inginkan. Catatan ini juga tercetak pada penalaran kriteria di
setiap laporan.

---

## 8. Yang TIDAK dikerjakan pipeline ini, dan cara mengerjakannya

Prediksi epitop dan struktur tidak diintegrasikan karena alat-alatnya berjalan
sebagai layanan daring atau butuh unduhan model besar dengan pendaftaran
akademik, keduanya tidak dapat dijalankan reprodusibel di dalam CI. Jalankan
manual, dan catat versi serta tanggal akses pada notebook:

| Alat | Fungsi | Lisensi | Cara pakai |
|---|---|---|---|
| IEDB Analysis Resource | prediksi dan cakupan populasi epitop | gratis, non-komersial | layanan web |
| BepiPred-3.0 | epitop sel B linear | gratis akademik (DTU) | layanan web / unduhan akademik |
| NetMHCpan-4.1 | epitop CTL (MHC kelas I) | gratis akademik (DTU) | layanan web / unduhan akademik |
| NetMHCIIpan-4.0 | epitop HTL (MHC kelas II) | gratis akademik (DTU) | layanan web / unduhan akademik |
| VaxiJen 2.0 | antigenisitas (ambang ≥ 0,4) | gratis | layanan web |
| ColabFold / AlphaFold2 | struktur protein | MIT (ColabFold), model DeepMind CC BY 4.0 | Google Colab gratis |
| SRAMP | prediksi m6A berbasis model | gratis | layanan web; pipeline ini memakai regex RRACH sebagai pengganti yang dapat diaudit |

**Pekerjaan terbuka** yang paling bernilai bila ada waktu tambahan:

1. Integrasi **LinearFold** dan **MXfold2** (keduanya MIT) sebagai pembanding
   lipatan terhadap ViennaRNA, dilaporkan berdampingan bukan dirata-ratakan.
2. Kalibrasi skor komposit terhadap efisiensi sirkularisasi hasil laboratorium
   begitu data pilot tersedia; sampai itu terjadi, bobotnya tetap `[PERKIRAAN]`.
3. Fitur **model bahasa RNA** (mis. RNA-FM) sebagai prediktor tambahan, wajib
   dilabeli eksperimental dan **dibandingkan** dengan baseline termodinamika,
   bukan menggantikannya.

---

## 9. Daftar perangkat lunak gratis yang benar-benar dipakai

| Perangkat lunak | Fungsi di pipeline ini | Lisensi |
|---|---|---|
| ViennaRNA 2.x | `RNAfold --circ`, `RNAplfold`, partition function, probabilitas pasangan basa | ViennaRNA License (gratis akademik/non-komersial) |
| primer3 / primer3-py | primer divergen dan BSJ, Tm, hairpin, homodimer | GPLv2 |
| Biopython | utilitas sekuens | Biopython License (mirip MIT) |
| Optuna | optimasi multi-objektif NSGA-II | MIT |
| Snakemake | alur kerja reprodusibel | MIT |
| NumPy, pandas | numerik dan tabel | BSD-3-Clause |
| NCBI BLAST+ (opsional) | homologi off-target | domain publik |
| STAR (opsional) | validasi silang BSJ | MIT |
| CIRI2 (opsional) | pemanggil circRNA pembanding | GPL |
| Streamlit | antarmuka web | Apache-2.0 |
| circhbsag (repositori ini) | perakitan konstruk, skor GO/NO-GO, optimasi kodon, probe RNase H, filter Bayes BSJ, laporan | MIT |

Jalankan `circhbsag licenses` untuk melihat versi yang benar-benar terpasang di
mesin Anda, tabel ini menyebut lisensi, sedangkan versinya diambil saat runtime.

---

## Sumber yang disitasi di dokumen ini

Sitasi di atas diperiksa langsung dari naskah aslinya:

- Wesselhoeft RA, Kowalski PS, Anderson DG (2018). Engineering circular RNA for
  potent and stable translation in eukaryotic cells. *Nature Communications*
  9:2629. doi:10.1038/s41467-018-05096-6

Pustaka lain yang relevan, tetapi tidak dipakai sebagai dasar keputusan desain di atas, antara lain: Chen dkk.
(2025) tentang sistem sirkularisasi tanpa scar (*Molecular Therapy Nucleic
Acids*); Wang dkk. (2024) tentang sistem PIE kimerik yang ditingkatkan
(*Theranostics* 14:5869); Meganck dkk. (2021) tentang rekayasa backsplicing
sintetis; Lee dkk. (2023) tentang ribozim intron grup I *Tetrahymena*; Singh,
Dujsikova, Mueller & Chen (2026) tentang pembuatan circRNA dengan ligasi
enzimatik (*Nucleic Acids Research* 54:gkag405); Hwang & Kim (2024) dan
Margvelani dkk. (2025) sebagai tinjauan mekanisme translasi circRNA; serta
Zhang dkk. (2025) tentang pengurangan kontaminan pada circRNA sintetis.
Periksa sendiri sebelum menyitasi, pipeline ini tidak memverifikasi klaim
isi pustaka selain yang dikutip eksplisit di atas.
