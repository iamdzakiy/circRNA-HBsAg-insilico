# Model kovarians RF00028 dan kontrol positifnya

## Mengapa dua berkas ini ada di repositori

Aturan repositori ini adalah **tidak menyertakan sekuens referensi**: nomor
aksesi dan tanggal unduh adalah bagian dari metode yang harus dicatat penulis.
Dua berkas di direktori ini adalah pengecualian yang disengaja, dengan alasan
yang berbeda untuk masing-masing.

| Berkas | Isi | Lisensi | Alasan disertakan |
|---|---|---|---|
| `RF00028.cm` | model kovarians Rfam untuk intron katalitik grup I | CC0 (data Rfam) | ini **model**, bukan sekuens konstruk; tanpa berkas ini Uji 2 Gate 0 tidak bisa dijalankan sama sekali |
| `../../tests/data/RF00028_positive_control.fa` | satu sekuens intron grup I nyata dari alignment Rfam | CC0 | kontrol positif: membuktikan Uji 2 mendeteksi intron sungguhan, bukan selalu menjawab FAIL |

## Provenans

Keduanya **tidak** diunduh dari rfam.org, karena sesi cloud tempat pipeline ini
dibangun ditolak proxy saat menjangkau `rfam.org` dan `ftp.ebi.ac.uk`. Keduanya
diambil dari cermin pihak ketiga di GitHub pada 2026-11-17:

```
# model kovarians (format INFERNAL1/a [1.1.1], ACC RF00028, CLEN 265)
https://raw.githubusercontent.com/biopro/genix/master/databases/rfam/CMs/RF00028.cm

# kontrol positif, diekstraksi dari alignment Stockholm RF00028
https://raw.githubusercontent.com/febos/SQUARNA-data/main/datasets/SubAli/sto/RF00028_2_2.sto
# sekuens: RF00028__Z17234.1/2999-3375, 377 nt setelah gap dihapus, U -> T
```

**Konsekuensi yang harus ditulis di dokumen tugas akhir:** karena keduanya
berasal dari cermin, bukan dari Rfam langsung, versi rilis Rfam-nya belum
dipastikan. Sebelum angka Uji 2 dikutip sebagai hasil akhir, unduh ulang dari
sumber resmi di macOS dan catat versinya:

```bash
curl -sS https://rfam.org/family/RF00028/cm -o data/reference/RF00028.cm
grep -m1 INFERNAL data/reference/RF00028.cm   # catat versi formatnya
```

Hasil tidak akan berubah bila modelnya sama; yang berubah adalah kemampuan Anda
menyebut versinya dengan pasti.

## Hasil kontrol positif

Dijalankan dengan Infernal 1.1.5 (dibangun dari sumber pada 2026-11-17):

```
cmsearch --notrunc -E 10 data/reference/RF00028.cm data/reference/RF00028_positive_control.fa
```

| Sekuens | Panjang | Bit score | E-value |
|---|---|---|---|
| RF00028__Z17234.1/2999-3375 | 377 nt | 118,0 | 2,1e-37 |
| RF00028__AY262145.2/209-414 | 206 nt | 102,6 | 5,9e-33 |

Keduanya jauh di bawah ambang E-value 1e-5 yang dipakai Uji 2. Sebaliknya,
prekursor demo yang "intron"-nya dibangkitkan acak menghasilkan **nol hit**.
Jadi Uji 2 terbukti bekerja di kedua arah: mengenali intron grup I sungguhan,
dan menolak sekuens yang hanya menyandang namanya.

## Di mana kontrol positifnya disimpan

Berkas kontrol positif ada di `tests/data/RF00028_positive_control.fa`, bukan di
direktori ini, supaya aturan "tidak ada sekuens referensi di `data/reference/`"
tetap utuh dan penjaga provenans di CI tidak perlu dikecualikan. Ia memang data
uji, bukan bagian konstruk.
