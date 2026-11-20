# data/reference/ — sekuens yang WAJIB Anda sediakan sendiri

Direktori ini sengaja dibiarkan kosong (kecuali tabel kodon). Pipeline tidak
mengikutkan sekuens virus, intron, maupun gen HBsAg, dengan dua alasan:

1. **Provenans.** Nomor aksesi, versi sekuens, koordinat potongan, dan tanggal
   unduh adalah bagian dari metode yang harus Anda catat dan pertanggungjawabkan
   di dokumen tugas akhir. Sekuens yang "sudah ada di repo" menghilangkan jejak itu.
2. **Pipeline tidak boleh mengarang.** Dengan `provenance_policy: strict` pada
   `config/parts.yaml`, pipeline **berhenti dengan galat** bila berkas di bawah
   belum ada, bukan melanjutkan dengan pengganti.

Untuk sekadar menguji pemasangan tanpa mengunduh apa pun, pakai
`config/parts.example.yaml` yang memakai data sintetis di `data/example/`.

---

## Berkas yang dibutuhkan

| Berkas | Isi | Catatan |
|---|---|---|
| `hbsag_S_protein.fa` | protein S-HBsAg (SHBs, ~226 aa), genotipe pilihan | masukan optimasi kodon; ORF dirancang ulang oleh pipeline, jadi yang dibutuhkan adalah **proteinnya**, bukan CDS-nya |
| `ires_cvb3.fa` | 5'UTR Coxsackievirus B3 (IRES tipe I) | catat koordinat awal/akhir potongan IRES |
| `intron_5half.fa` | separuh 5' intron grup I permuted | berisi IGS; sesuaikan koordinat `features: IGS` di `parts.yaml` |
| `intron_3half.fa` | separuh 3' intron grup I permuted | berisi IGS untuk P10 |
| `E1.fa`, `E2.fa` | eksón 1 dan eksón 2 | pasangan ini membentuk sambungan balik pada lingkaran |
| `homology_arm_5.fa`, `homology_arm_3.fa` | lengan homologi | panjangnya dioptimasi oleh `circhbsag optimize` |
| `codon_usage_human.tsv` | tabel penggunaan kodon inang | sudah ada, tetapi ditandai `[PERKIRAAN]` — ganti (lihat bawah) |

Sumber yang lazim: GenBank/ENA untuk sekuens virus; peta plasmid PIE
terpublikasi (mis. di Addgene) untuk intron, eksón, dan lengan homologi, karena
koordinat intron grup I yang tepat lebih mudah diambil dari konstruk yang sudah
terbukti bekerja daripada dirakit ulang dari genom.

## Cara mengunduh

Dengan `ncbi-datasets` atau E-utilities (ganti `<AKSESI>`):

```bash
# satu rekaman FASTA dari nuccore
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=<AKSESI>&rettype=fasta&retmode=text" \
  -o data/reference/unduhan.fa

# memotong sub-wilayah (koordinat 1-based, inklusif) dengan seqkit
seqkit subseq -r <AWAL>:<AKHIR> data/reference/unduhan.fa > data/reference/ires_cvb3.fa
```

Setelah mengunduh, **catat di `config/parts.yaml`**: ganti setiap
`accession: "ISI_SENDIRI"` dengan aksesi dan versi yang benar-benar dipakai.
CI punya langkah `provenance-guard` yang gagal bila ada FASTA di direktori ini
sementara `parts.yaml` masih memuat placeholder `ISI_SENDIRI`.

## Tabel penggunaan kodon

`codon_usage_human.tsv` yang disertakan ditandai `[PERKIRAAN]`: nilainya adalah
tabel Kazusa *Homo sapiens* yang lazim direproduksi, tetapi belum diunduh ulang
dan diverifikasi oleh pipeline. Nilai **CAI relatif** antar-varian desain stabil
terhadap pilihan tabel; nilai **CAI absolut** tidak. Sebelum melaporkan angka
CAI absolut, lakukan salah satu:

```bash
# (a) unduh tabel resmi dari https://www.kazusa.or.jp/codon/ (Homo sapiens),
#     ubah ke format: codon <tab> aa <tab> frequency_per_thousand
# (b) atau hitung dari kumpulan CDS inang sendiri:
circhbsag codon-table --fasta cds_inang.fa --out data/reference/codon_usage_human.tsv
```

Opsi (b) menghasilkan berkas berlabel `[HASIL IN SILICO]` karena dihitung dari
data yang Anda sediakan.

## Setelah semuanya terisi

```bash
circhbsag tahap0 --parts config/parts.yaml --out results   # harus jalan tanpa galat
snakemake --cores 4                                        # alur lengkap
```

Laporan akan berlabel `[HASIL IN SILICO]`, bukan `[DEMO — SEKUENS PLACEHOLDER]`.


---

## Sekuens HBsAg nyata (ditambahkan 2026-11-20)

| Berkas | Isi | Asal |
|---|---|---|
| `hbsag_S_protein.fa` | protein S-HBsAg 226 aa | protein_id AAB59969.1 `[LITERATUR]` |
| `hbsag_S_CDS_AAB59969.1.fa` | CDS 681 nt (226 kodon + stop), posisi `157..837` pada catatan induk | disalin jek dari catatan GenBank; terjemahannya sama persis dengan terjemahan di catatan itu |

**Aksesi genom induk `[PERLU VERIFIKASI]`:** jek menyebut D00329.1, tetapi protein ini berbeda di 20 posisi dari protein D00329.1,
jadi label D00329.1 tidak dipakai sampai jek mengonfirmasi aksesi yang benar.
**Koordinat:** CDS di catatan itu `157..837`, bukan 155–835 (bergeser 2 nt, salah kerangka).
**Genotipe:** B menurut pohon filogeni lokal pada D00329.1, jadi genotipe untuk protein ini juga `[PERLU VERIFIKASI]`.
Bagian lain (IRES CVB3, intron, ekson) masih harus disediakan, jadi `config/parts.yaml` (strict) belum bisa berjalan penuh.
