#!/usr/bin/env python3
"""Bangkitkan DATA CONTOH SINTETIS untuk menguji perangkat lunak.

PENTING — baca sebelum memakai:

  Berkas yang dihasilkan skrip ini BUKAN sekuens biologis. Tidak ada satu pun
  nukleotida di sini yang berasal dari Coxsackievirus B3, dari intron grup I
  Anabaena, atau dari virus hepatitis B. Semuanya dibangkitkan oleh generator
  acak dengan seed tetap, dirancang hanya agar PANJANG dan SIFAT KASAR-nya
  (rentang GC, adanya batang-lengkung, adanya IGS pendek) cukup mirip sehingga
  pipeline, tes, CI, dan antarmuka web dapat dijalankan tanpa jaringan.

  Angka apa pun yang dihasilkan dari data ini HANYA sah sebagai uji perangkat
  lunak, dan dilabeli [DEMO — SEKUENS PLACEHOLDER] oleh modul report. Angka tersebut
  TIDAK BOLEH masuk ke dokumen tugas akhir sebagai hasil ilmiah.

  Untuk hasil sungguhan, ganti berkas di data/reference/ dengan sekuens
  terverifikasi (lihat data/reference/README.md) dan pakai config/parts.yaml
  dengan provenance_policy: strict.

Pemakaian:
    python data/example/make_example.py            # tulis ke data/example/
"""
from __future__ import annotations

import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEED = 20260101


def rand_seq(n: int, gc: float, rng: random.Random) -> str:
    """Sekuens acak dengan GC kira-kira `gc`, tanpa homopolimer >= 4."""
    out: list[str] = []
    while len(out) < n:
        b = rng.choice("GC") if rng.random() < gc else rng.choice("AT")
        if len(out) >= 3 and out[-1] == out[-2] == out[-3] == b:
            continue
        out.append(b)
    return "".join(out)


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def hairpin(stem_len: int, loop_len: int, rng: random.Random, gc: float = 0.6) -> str:
    """Batang-lengkung sederhana: batang berpasangan sempurna + lengkung acak."""
    stem = rand_seq(stem_len, gc, rng)
    loop = rand_seq(loop_len, 0.35, rng)
    return stem + loop + revcomp(stem)


def structured_block(total: int, rng: random.Random) -> str:
    """Blok berstruktur: beberapa batang-lengkung dipisahkan penghubung pendek.

    Dipakai sebagai pengganti sintetis untuk wilayah IRES, agar pipeline punya
    sesuatu yang berstruktur untuk diukur. Bukan IRES sungguhan.
    """
    parts: list[str] = []
    while sum(len(p) for p in parts) < total:
        parts.append(hairpin(rng.randint(7, 12), rng.randint(4, 9), rng))
        parts.append(rand_seq(rng.randint(5, 14), 0.4, rng))
    return "".join(parts)[:total]


def write_fa(path: Path, name: str, desc: str, seq: str) -> None:
    body = "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60))
    path.write_text(f">{name} {desc}\n{body}\n")
    print(f"  {path.name}: {len(seq)} nt")


def main() -> None:
    rng = random.Random(SEED)
    tag = "SINTETIS-UJI BUKAN-SEKUENS-NYATA"

    print("Menulis data contoh sintetis (seed", SEED, "):")

    # Separuh intron: IGS 6 nt di awal separuh 5' (anotasi features pada parts.yaml),
    # sisanya blok berstruktur agar mirip intron grup I secara kasar.
    igs5 = "GGTGGT"
    intron5 = igs5 + structured_block(174, rng)
    igs3 = "ACCACC"
    intron3 = structured_block(174, rng) + igs3
    write_fa(HERE / "intron_5half.fa", "DEMO_intron_5half", f"{tag} 180 nt, IGS pada 0..6", intron5)
    write_fa(HERE / "intron_3half.fa", "DEMO_intron_3half", f"{tag} 180 nt, IGS pada 174..180", intron3)

    # Eksón pendek yang membentuk sambungan balik.
    write_fa(HERE / "E1.fa", "DEMO_E1", f"{tag} eksón 1", rand_seq(30, 0.45, rng))
    write_fa(HERE / "E2.fa", "DEMO_E2", f"{tag} eksón 2", rand_seq(30, 0.45, rng))

    # Lengan homologi yang saling komplementer (meniru fungsi penyatuan ujung).
    ha5 = rand_seq(40, 0.45, rng)
    write_fa(HERE / "HA5.fa", "DEMO_HA5", f"{tag} lengan homologi 5'", ha5)
    write_fa(HERE / "HA3.fa", "DEMO_HA3", f"{tag} lengan homologi 3' (komplemen HA5)", revcomp(ha5))

    # Blok berstruktur sebagai pengganti IRES (panjang ~600 nt seperti IRES pikorna).
    write_fa(HERE / "ires_demo.fa", "DEMO_IRES", f"{tag} blok berstruktur 600 nt", structured_block(600, rng))

    # Protein demo: BUKAN HBsAg. Dibangkitkan acak dari asam amino umum.
    #
    # Panjangnya ditetapkan 226 residu, yaitu panjang protein S-HBsAg kecil
    # (SHBs), SEMATA agar dimensi konstruk demo (panjang ORF ~681 nt dan
    # panjang lingkaran total) sebanding dengan konstruk sungguhan, sehingga
    # angka MFE/nt, beban dsRNA, dan ukuran amplikon yang dihasilkan demo
    # berada pada orde yang sama. Komposisi residunya tetap acak: tidak satu
    # pun berasal dari HBsAg, dan determinan "a" tidak terwakili.
    aa_pool = "ACDEFGHIKLMNPQRSTVWY"
    protein = "M" + "".join(rng.choice(aa_pool) for _ in range(225))
    (HERE / "demo_antigen_protein.fa").write_text(
        f">DEMO_ANTIGEN {tag} protein acak 226 aa — BUKAN HBsAg, "
        f"panjang disamakan dengan SHBs hanya agar dimensi konstruk sebanding\n"
        + "\n".join(protein[i : i + 60] for i in range(0, len(protein), 60))
        + "\n"
    )
    print(f"  demo_antigen_protein.fa: {len(protein)} aa")
    print("Selesai. Semua berkas berlabel SINTETIS-UJI.")


if __name__ == "__main__":
    main()
