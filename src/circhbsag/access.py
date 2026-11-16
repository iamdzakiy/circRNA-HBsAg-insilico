"""Aksesibilitas lokal dengan RNAplfold (ViennaRNA) dan deteksi dsRNA panjang.

RNAplfold menghitung probabilitas suatu jendela beruntun tetap tak-berpasangan
dalam ansambel lipatan lokal; ini adalah ukuran aksesibilitas yang relevan untuk
(i) situs sambungan E1-E2, (ii) helai P1/P10 intron grup I, dan (iii) situs
hibridisasi primer/probe RNase H.

Lisensi ViennaRNA: ViennaRNA License (gratis untuk pemakaian akademik dan non-komersial).
"""
from __future__ import annotations

import RNA

from .seqio import to_rna


def plfold_unpaired(
    seq: str,
    *,
    ulength: int = 10,
    window: int = 80,
    max_bp_span: int = 40,
    temperature_c: float = 37.0,
) -> list[list[float]]:
    """Matriks aksesibilitas RNAplfold.

    Mengembalikan `up[i][u]` = probabilitas bahwa u nukleotida yang berakhir di
    posisi i (1-based) semuanya tak-berpasangan. Indeks 0 tidak dipakai.
    """
    prev = RNA.cvar.temperature
    try:
        RNA.cvar.temperature = float(temperature_c)
        up = RNA.pfl_fold_up(to_rna(seq), int(ulength), int(window), int(max_bp_span))
    finally:
        RNA.cvar.temperature = prev
    return up


def accessibility(
    seq: str,
    regions: dict[str, tuple[int, int]],
    *,
    ulength: int = 10,
    window: int = 80,
    max_bp_span: int = 40,
    temperature_c: float = 37.0,
) -> dict[str, dict]:
    """Aksesibilitas rata-rata untuk setiap wilayah bernama.

    `regions` memakai koordinat 1-based inklusif: {"nama": (awal, akhir)}.
    Nilai yang dilaporkan adalah rerata dan minimum dari P(u nt tak-berpasangan)
    pada seluruh posisi di wilayah tersebut, dengan u = `ulength`.
    """
    up = plfold_unpaired(
        seq, ulength=ulength, window=window, max_bp_span=max_bp_span, temperature_c=temperature_c
    )
    n = len(seq)
    out: dict[str, dict] = {}
    for name, (start, end) in regions.items():
        start, end = max(1, start), min(n, end)
        vals = []
        for i in range(start, end + 1):
            row = up[i] if i < len(up) else None
            if not row:
                continue
            u = min(ulength, i)
            v = row[u]
            if v is not None:
                vals.append(float(v))
        if not vals:
            out[name] = {"start": start, "end": end, "n": 0, "mean": None, "min": None}
            continue
        out[name] = {
            "start": start,
            "end": end,
            "n": len(vals),
            "ulength": ulength,
            "mean": round(sum(vals) / len(vals), 6),
            "min": round(min(vals), 6),
            "max": round(max(vals), 6),
        }
    return out


def helices_from_structure(structure: str) -> list[dict]:
    """Uraikan struktur titik-kurung menjadi daftar heliks (batang) beruntun.

    Mengembalikan daftan {"i5": awal_5', "j3": akhir_3', "length": bp} dengan
    koordinat 1-based. Heliks dipotong bila ada bulge/mismatch (pasangan tidak
    bersebelahan).
    """
    stack: list[int] = []
    pairs: list[tuple[int, int]] = []
    for k, ch in enumerate(structure, start=1):
        if ch == "(":
            stack.append(k)
        elif ch == ")":
            if stack:
                pairs.append((stack.pop(), k))
    pairs.sort()
    helices: list[dict] = []
    pairset = dict(pairs)
    used: set[int] = set()
    for i, j in pairs:
        if i in used:
            continue
        length = 0
        a, b = i, j
        while pairset.get(a) == b:
            used.add(a)
            length += 1
            a, b = a + 1, b - 1
        helices.append({"i5": i, "j3": j, "length": length})
    return sorted(helices, key=lambda h: -h["length"])


def longest_dsrna(structure: str) -> dict:
    """Heliks beruntun terpanjang (proksi dsRNA pemicu RIG-I/PKR).

    [LITERATUR] dsRNA >~30-40 bp beruntun dikaitkan dengan aktivasi PKR dan
    pengenalan oleh RIG-I/MDA5; pipeline ini memakai ambang yang dapat diubah
    pada config (default 30 bp).
    """
    hel = helices_from_structure(structure)
    if not hel:
        return {"longest_bp": 0, "helices_ge_threshold": [], "n_helices": 0}
    return {
        "longest_bp": hel[0]["length"],
        "longest_helix": hel[0],
        "n_helices": len(hel),
        "top_helices": hel[:10],
    }


def dsrna_violations(structure: str, threshold_bp: int = 30) -> list[dict]:
    """Semua heliks beruntun >= threshold_bp."""
    return [h for h in helices_from_structure(structure) if h["length"] >= threshold_bp]
