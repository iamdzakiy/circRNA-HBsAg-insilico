"""Desain primer dan probe dengan primer3-py (lisensi: GPLv2, gratis).

Tiga jenis oligo yang dirancang:

  1. Primer divergen (outward-facing) pada lingkaran. Pada templat linear kedua
     primer menjauh satu dari yang lain sehingga tidak ada produk; pada lingkaran
     keduanya menghasilkan amplikon yang melintasi sambungan balik (BSJ).
  2. Primer BSJ. Satu primer diletakkan menumpang tepat di atas sambungan
     E1-E2 sehingga hanya menempel pada molekul sirkular.
  3. Probe RNase H. Oligo DNA yang berhibridisasi dengan RNA di satu posisi
     yang diketahui; RNase H memotong untai RNA pada heteroduplex, sehingga
     lingkaran menjadi satu fragmen linear berukuran yang dapat diprediksi
     (uji nicking, pembeda lingkaran vs concatemer linear).

Primer reverse untuk transkripsi in vitro dengan 2'-O-metil pada dua nukleotida
terakhir dirancang oleh `ivt_reverse_primer`: modifikasi tersebut menekan
penambahan nukleotida non-templat oleh T7 di ujung 3' templat.
"""
from __future__ import annotations

import primer3

from . import access
from .seqio import gc, revcomp, to_dna

PRIMER3_LICENSE = "GPLv2 (primer3 / primer3-py) — gratis, open-source"

DEFAULT_GLOBALS: dict = {
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_MAX_SIZE": 25,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MIN_TM": 57.0,
    "PRIMER_MAX_TM": 63.0,
    "PRIMER_MIN_GC": 40.0,
    "PRIMER_MAX_GC": 60.0,
    "PRIMER_MAX_POLY_X": 4,
    "PRIMER_MAX_SELF_ANY": 8.0,
    "PRIMER_MAX_SELF_END": 3.0,
    "PRIMER_PAIR_MAX_COMPL_ANY": 8.0,
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_DNA_CONC": 50.0,
    "PRIMER_MAX_NS_ACCEPTED": 0,
    "PRIMER_GC_CLAMP": 1,
}


def _design(seq: str, *, product_range: list[list[int]], target: tuple[int, int] | None,
            n: int = 5, extra: dict | None = None) -> dict:
    seq_args = {"SEQUENCE_ID": "amplicon", "SEQUENCE_TEMPLATE": to_dna(seq)}
    if target:
        seq_args["SEQUENCE_TARGET"] = [int(target[0]), int(target[1])]
    g = {**DEFAULT_GLOBALS, "PRIMER_NUM_RETURN": n, "PRIMER_PRODUCT_SIZE_RANGE": product_range}
    if extra:
        g.update(extra)
    return primer3.bindings.design_primers(seq_args, g)


def _rc(x):
    return x[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))


def _pairs(res: dict, n: int) -> list[dict]:
    out = []
    for i in range(n):
        if f"PRIMER_LEFT_{i}_SEQUENCE" not in res:
            break
        out.append(
            {
                "rank": i,
                "left_seq": res[f"PRIMER_LEFT_{i}_SEQUENCE"],
                "left_start": res[f"PRIMER_LEFT_{i}"][0],
                "left_len": res[f"PRIMER_LEFT_{i}"][1],
                "left_tm": round(res[f"PRIMER_LEFT_{i}_TM"], 2),
                "left_gc": round(res[f"PRIMER_LEFT_{i}_GC_PERCENT"], 2),
                "right_seq": res[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                "right_start": res[f"PRIMER_RIGHT_{i}"][0],
                "right_len": res[f"PRIMER_RIGHT_{i}"][1],
                "right_tm": round(res[f"PRIMER_RIGHT_{i}_TM"], 2),
                "right_gc": round(res[f"PRIMER_RIGHT_{i}_GC_PERCENT"], 2),
                "product_size": res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                "pair_penalty": round(res[f"PRIMER_PAIR_{i}_PENALTY"], 4),
            }
        )
    return out


def divergent_primers(asm, *, n: int = 5, flank: int = 150, min_product: int = 100) -> dict:
    """Primer divergen: dirancang pada lingkaran yang dirotasi 180 derajat dari BSJ.

    Caranya: rotasi lingkaran sehingga titik ANTIPODAL BSJ berada di tengah
    amplikon. Primer3 merancang pasangan konvergen pada sekuens terotasi itu;
    pada lingkaran asli (atau pada prekursor linear) pasangan yang sama
    menghadap ke arah berlawanan, sehingga hanya templat sirkular memberi produk.
    """
    circle = asm.circle
    nlen = len(circle)
    anti = (asm.bsj_position + nlen // 2) % nlen  # titik antipodal BSJ
    flank = min(flank, nlen // 2 - 1)
    window = "".join(circle[(anti - flank + k) % nlen] for k in range(2 * flank))
    res = _design(
        window,
        product_range=[[min_product, len(window)]],
        target=(flank - 10, 20),
        n=n,
    )
    pairs = _pairs(res, n)
    for p in pairs:
        # Posisi pada lingkaran asli (1-based)
        p["left_start_in_circle"] = ((anti - flank + p["left_start"]) % nlen) + 1
        p["right_start_in_circle"] = ((anti - flank + p["right_start"]) % nlen) + 1
        p["expected_product_on_circle_nt"] = nlen - p["product_size"] + p["left_len"] + p["right_len"]
        # left_seq/right_seq berorientasi KONVERGEN: produk pendek (product_size) ada pada lingkaran DAN pada
        # prekursor linear (kontrol total RNA). Pasangan divergen sejati memakai komplemen terbalik keduanya
        # (maju = rc(right_seq), mundur = rc(left_seq)) dengan produk panjang pada lingkaran saja.
        p["convergent_product_nt"] = p["product_size"]
        p["divergent_fwd_5to3"] = _rc(p["right_seq"])
        p["divergent_rev_5to3"] = _rc(p["left_seq"])
        p["expected_product_on_linear"] = (
            f"{p['product_size']} nt (konvergen: ada pada linear); pasangan divergen sejati: tidak ada")
    return {
        "strategy": "konvergen (total RNA, kontrol) + versi divergen sejati bila urutan dibalik",
        "antipodal_position": anti,
        "window_length": len(window),
        "primer3_license": PRIMER3_LICENSE,
        "n_returned": len(pairs),
        "explain": res.get("PRIMER_LEFT_EXPLAIN", ""),
        "pairs": pairs,
    }


def bsj_primers(asm, *, n: int = 5, flank: int = 150, overlap: int = 8) -> dict:
    """Primer yang menumpang sambungan balik (BSJ-spanning).

    Primer3 dipakai untuk pasangan di sekitar BSJ pada amplikon yang telah
    dirotasi (BSJ di tengah). Selain itu dibuat satu primer forward manual yang
    berpusat pada sambungan dengan `overlap` nt di tiap sisi; primer semacam ini
    hanya dapat menempel bila E1 dan E2 benar-benar tersambung.
    """
    from .construct import bsj_amplicon

    amp = bsj_amplicon(asm, flank=flank)
    seq = amp["amplicon"]
    j = amp["junction_offset_in_amplicon"]
    res = _design(seq, product_range=[[80, len(seq)]], target=(j - 5, 10), n=n)
    pairs = _pairs(res, n)
    for p in pairs:
        p["spans_junction"] = p["left_start"] < j < (p["right_start"] + 1)

    manual = []
    for ov in (overlap, overlap + 4, overlap + 8):
        fwd = seq[j - ov : j + ov]
        if len(fwd) < 2 * ov:
            continue
        manual.append(
            {
                "name": f"BSJ_fwd_overlap{ov}",
                "seq": fwd,
                "length": len(fwd),
                "nt_on_E1_side": ov,
                "nt_on_E2_side": ov,
                "tm_primer3": round(primer3.calc_tm(fwd), 2),
                "gc": round(gc(fwd), 3),
                "hairpin_tm": round(primer3.calc_hairpin_tm(fwd), 2),
                "homodimer_tm": round(primer3.calc_homodimer_tm(fwd), 2),
            }
        )
    return {
        "strategy": "primer menumpang BSJ (junction-spanning)",
        "amplicon_length": len(seq),
        "junction_offset": j,
        "primer3_license": PRIMER3_LICENSE,
        "pairs": pairs,
        "junction_spanning_primers": manual,
        "note": "Verifikasi sambungan dengan sekuensing Sanger amplikon, bukan hanya ukuran pita.",
    }


def ivt_reverse_primer(template_dna: str, *, length: int = 22, two_o_methyl_nt: int = 2) -> dict:
    """Primer reverse untuk membuat templat IVT, dengan 2'-O-metil di ujung 5'.

    Catatan kimia: primer reverse mendefinisikan ujung 3' templat untai atas.
    Modifikasi 2'-O-metil pada dua nukleotida 5' primer reverse menghasilkan
    ujung templat yang menekan penambahan nukleotida non-templat ("untemplated
    addition") oleh T7 RNA polimerase, sehingga ujung 3' transkrip lebih homogen.
    """
    rev = revcomp(to_dna(template_dna))[:length]
    mod = "".join(f"m{b}" for b in rev[:two_o_methyl_nt]) + rev[two_o_methyl_nt:]
    return {
        "primer_reverse_5to3": rev,
        "with_modifications": mod,
        "modification": f"2'-O-metil pada {two_o_methyl_nt} nukleotida pertama (ujung 5')",
        "length": len(rev),
        "tm": round(primer3.calc_tm(rev), 2),
        "gc": round(gc(rev), 3),
        "hairpin_tm": round(primer3.calc_hairpin_tm(rev), 2),
        "purpose": "menekan penambahan nukleotida non-templat T7 di ujung 3' transkrip",
    }


def rnase_h_probes(
    asm, *, n: int = 3, probe_len: int = 20, temperature_c: float = 37.0, min_sep: int = 60
) -> dict:
    """Pilih probe DNA untuk uji nicking RNase H pada lingkaran.

    Kriteria pemilihan (dapat diaudit):
      - Tm 55-70 C (primer3 calc_tm)
      - GC 40-60%
      - tanpa homopolimer >= 5
      - aksesibilitas RNAplfold tinggi pada RNA lingkaran (hibridisasi lebih mudah)
      - jarak antar-probe >= `min_sep` nt, dan >= 40 nt dari BSJ agar produk
        linearisasi tetap mudah dibedakan dari fragmen BSJ
    """
    circle = asm.circle
    nlen = len(circle)
    up = access.plfold_unpaired(circle, ulength=probe_len, window=max(80, probe_len * 4),
                                temperature_c=temperature_c)
    cands = []
    for start in range(0, nlen - probe_len):
        sub = circle[start : start + probe_len]
        g = gc(sub)
        if not 0.40 <= g <= 0.60:
            continue
        if max(len(x) for x in __import__("re").findall(r"(.)\1*", sub)) >= 5:
            continue
        tm = primer3.calc_tm(revcomp(sub))
        if not 55.0 <= tm <= 70.0:
            continue
        end1 = start + probe_len  # 1-based akhir
        row = up[end1] if end1 < len(up) else None
        acc = float(row[probe_len]) if row and row[probe_len] is not None else 0.0
        dist_bsj = min(abs(start + probe_len // 2 - asm.bsj_position),
                       nlen - abs(start + probe_len // 2 - asm.bsj_position))
        if dist_bsj < 40:
            continue
        cands.append(
            {
                "start_1based": start + 1,
                "end_1based": start + probe_len,
                "rna_target": sub,
                "probe_dna_5to3": revcomp(sub),
                "tm": round(tm, 2),
                "gc": round(g, 3),
                "accessibility_plfold": round(acc, 6),
                "distance_to_bsj_nt": dist_bsj,
                "expected_linear_product_nt": nlen,
            }
        )
    cands.sort(key=lambda d: -d["accessibility_plfold"])
    chosen: list[dict] = []
    for c in cands:
        if all(abs(c["start_1based"] - x["start_1based"]) >= min_sep for x in chosen):
            chosen.append(c)
        if len(chosen) >= n:
            break
    return {
        "strategy": "probe DNA untuk nicking RNase H pada circRNA",
        "temperature_c": temperature_c,
        "probe_length": probe_len,
        "n_candidates_passing_filters": len(cands),
        "probes": chosen,
        "interpretation": (
            "Lingkaran utuh + satu probe -> satu pita linear tunggal berukuran = panjang lingkaran "
            f"({nlen} nt). Concatemer atau prekursor linear menghasilkan pola pita berbeda; "
            "RNA linear tanpa situs probe tidak berubah."
        ),
    }


def design_all(asm, **kw) -> dict:
    """Jalankan seluruh desain oligo untuk satu konstruk."""
    return {
        "divergent": divergent_primers(asm, **kw.get("divergent", {})),
        "bsj": bsj_primers(asm, **kw.get("bsj", {})),
        "ivt_reverse": ivt_reverse_primer(asm.precursor, **kw.get("ivt_reverse", {})),
        "rnase_h": rnase_h_probes(asm, **kw.get("rnase_h", {})),
    }
