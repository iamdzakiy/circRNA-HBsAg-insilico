"""Pemindaian motif sekuens: situs restriksi, terminator T7, m6A RRACH, GC lokal.

Semua motif didefinisikan sebagai konstanta eksplisit agar dapat diaudit.
Sumber motif ditandai di docstring masing-masing fungsi:
  [LITERATUR] = definisi baku dari pustaka/katalog enzim
  [PERKIRAAN] = heuristik yang ditetapkan pipeline ini (ambang dapat diubah)
"""
from __future__ import annotations

import re

from .seqio import gc_windows, revcomp, to_dna

# [LITERATUR] situs pengenalan Tipe IIS yang dipakai pada perakitan Golden Gate.
RESTRICTION_SITES: dict[str, str] = {
    "BsaI": "GGTCTC",
    "BsmBI_Esp3I": "CGTCTC",
    "SapI": "GCTCTTC",
    "BbsI": "GAAGAC",
    "AarI": "CACCTGC",
}

# [LITERATUR] Terminator T7 kelas I (Tphi) — inti jepit rambut + runtun U.
# Regex di bawah mendeteksi inti terminator Tphi yang paling dikenal; motif
# longgar (jepit rambut GC diikuti >=4 T) ditangani terpisah.
T7_TERMINATOR_CORE = "CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG"
T7_TERMINATOR_SEED = "GGGGCCTCTAAACGGGTCTTGAGGGG"

# [LITERATUR] Promotor T7 kanonik (untai atas).
T7_PROMOTER = "TAATACGACTCACTATA"

# [LITERATUR] Konsensus situs metilasi m6A: RRACH (R=A/G, H=A/C/U).
RRACH = re.compile(r"(?=([AG][AG]AC[ACT]))")

# [PERKIRAAN] motif imunostimulan/destabilisasi yang lazim dihindari pada desain mRNA.
AVOID_MOTIFS: dict[str, str] = {
    "ARE_AUUUA": "ATTTA",          # elemen kaya AU, destabilisasi mRNA
    "polyU_6": "TTTTTT",           # runtun U, terminasi/pause T7
    "CpG_island_seed": "CGCGCGCG",  # runtun CpG, pemicu TLR9
    "SpliceDonor_GTAAGT": "GTAAGT",  # donor splice kriptik
    "SpliceAcceptor_AGGTRAGT": "AGGTAAGT",
    "KozakConflict_internal": "GCCACCATGG",
}


def find_all(seq: str, motif: str, *, both_strands: bool = True) -> list[dict]:
    """Semua kemunculan motif (0-based, overlap diperbolehkan)."""
    s = to_dna(seq)
    hits = []
    for strand, pat in (("+", to_dna(motif)),) + (
        (("-", revcomp(to_dna(motif))),) if both_strands and revcomp(motif) != motif else ()
    ):
        for m in re.finditer(f"(?=({re.escape(pat)}))", s):
            hits.append({"start": m.start(), "end": m.start() + len(pat), "strand": strand, "motif": motif})
    return sorted(hits, key=lambda d: d["start"])


def scan_restriction(seq: str, enzymes: list[str] | None = None) -> dict[str, list[dict]]:
    """[LITERATUR] Situs Tipe IIS internal yang akan merusak perakitan Golden Gate."""
    enzymes = enzymes or list(RESTRICTION_SITES)
    return {e: find_all(seq, RESTRICTION_SITES[e]) for e in enzymes}


def scan_t7_terminator(seq: str) -> list[dict]:
    """[LITERATUR] Terminator T7 prematur: inti Tphi atau jepit rambut GC + >=5 T."""
    s = to_dna(seq)
    hits = [dict(h, kind="Tphi_core") for h in find_all(s, T7_TERMINATOR_SEED)]
    # Pola longgar: batang kaya GC (>=6 nt) diikuti <=8 nt lengkung lalu runtun >=5 T.
    for m in re.finditer(r"([GC]{6,})([ACGT]{3,8})([GC]{6,})(T{5,})", s):
        stem1, _loop, stem2, _u = m.groups()
        if _hybridizes(stem1, stem2):
            hits.append(
                {
                    "start": m.start(),
                    "end": m.end(),
                    "strand": "+",
                    "motif": m.group(0),
                    "kind": "intrinsic_like",
                }
            )
    return sorted(hits, key=lambda d: d["start"])


def _hybridizes(stem1: str, stem2: str, min_bp: int = 5) -> bool:
    rc = revcomp(stem2)
    n = min(len(stem1), len(rc))
    match = sum(1 for a, b in zip(stem1[-n:], rc[:n], strict=False) if a == b)
    return match >= min_bp


def scan_m6a_rrach(seq: str) -> list[dict]:
    """[LITERATUR] Semua situs RRACH; A yang termetilasi ada di posisi ke-3 motif."""
    s = to_dna(seq)
    out = []
    for m in RRACH.finditer(s):
        start = m.start()
        out.append(
            {
                "start": start,
                "end": start + 5,
                "motif": m.group(1),
                "a_position": start + 2,  # 0-based posisi adenosin sasaran
            }
        )
    return out


def m6a_density(seq: str) -> float:
    """Situs RRACH per 100 nt."""
    return 100.0 * len(scan_m6a_rrach(seq)) / len(seq) if seq else 0.0


def scan_avoid_motifs(seq: str) -> dict[str, list[dict]]:
    """[PERKIRAAN] Motif yang dihindari pada ORF (lihat AVOID_MOTIFS)."""
    return {k: find_all(seq, v, both_strands=False) for k, v in AVOID_MOTIFS.items()}


def gc_window_violations(
    seq: str, window: int = 50, step: int = 10, lo: float = 0.30, hi: float = 0.70
) -> list[dict]:
    """[PERKIRAAN] Jendela dengan GC di luar rentang yang dapat disintesis."""
    return [
        {"start": i, "end": i + window, "gc": round(v, 4)}
        for i, v in gc_windows(seq, window, step)
        if v < lo or v > hi
    ]


def start_codon_context(seq: str, atg_pos: int, flank: int = 15) -> dict:
    """Konteks di sekitar start codon ORF (0-based atg_pos)."""
    s = to_dna(seq)
    return {
        "atg_position": atg_pos,
        "upstream": s[max(0, atg_pos - flank) : atg_pos],
        "codon": s[atg_pos : atg_pos + 3],
        "downstream": s[atg_pos + 3 : atg_pos + 3 + flank],
        "is_atg": s[atg_pos : atg_pos + 3] == "ATG",
    }


def upstream_atg_in_frame(seq: str, orf_start: int, *, within: int | None = None) -> list[int]:
    """ATG hulu yang sekerangka dengan ORF (risiko inisiasi prematur).

    `within` membatasi pencarian ke sejumlah nukleotida terakhir sebelum ORF.
    Ini penting untuk konstruk ber-IRES: IRES virus pikorna secara alami memuat
    banyak AUG di dalam domainnya, sehingga mencari ATG sekerangka di SELURUH
    hulu akan selalu menemukan sesuatu dan tidak informatif. Yang relevan secara
    mekanistik adalah AUG sekerangka tepat di hulu kodon awal, yang akan
    memperpanjang ujung-N produk translasi.
    """
    s = to_dna(seq)
    lo = 0 if within is None else max(0, orf_start - within)
    first = lo + ((orf_start - lo) % 3)
    return [i for i in range(first, orf_start, 3) if s[i : i + 3] == "ATG"]


def upstream_atg_out_of_frame(seq: str, orf_start: int, *, within: int = 30) -> list[int]:
    """ATG hulu TIDAK sekerangka dalam `within` nt (risiko ORF bersaing)."""
    s = to_dna(seq)
    lo = max(0, orf_start - within)
    return [
        i for i in range(lo, orf_start) if s[i : i + 3] == "ATG" and (orf_start - i) % 3 != 0
    ]


def summarize(seq: str, *, enzymes: list[str] | None = None) -> dict:
    """Ringkasan seluruh pemindaian motif untuk satu sekuens."""
    rs = scan_restriction(seq, enzymes)
    return {
        "length": len(seq),
        "restriction": {k: len(v) for k, v in rs.items()},
        "restriction_hits": rs,
        "t7_terminator": scan_t7_terminator(seq),
        "m6a_sites": len(scan_m6a_rrach(seq)),
        "m6a_density_per_100nt": round(m6a_density(seq), 3),
        "avoid_motifs": {k: len(v) for k, v in scan_avoid_motifs(seq).items()},
        "gc_violations": gc_window_violations(seq),
    }
