"""Lipatan RNA sirkular dan linear dengan ViennaRNA (lisensi: ViennaRNA License, gratis).

Fungsi di sini memakai binding Python resmi ViennaRNA (`import RNA`), yang setara
dengan `RNAfold --circ` pada antarmuka baris perintah. Tidak ada angka yang
dihardcode: semua nilai berasal dari perhitungan ViennaRNA saat dijalankan.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import RNA

from .seqio import to_rna

VIENNA_VERSION = RNA.__version__


@dataclass
class FoldResult:
    """Hasil lipatan satu sekuens pada satu suhu."""

    sequence_length: int
    temperature_c: float
    circular: bool
    mfe_structure: str
    mfe_kcal: float
    mfe_per_nt: float
    ensemble_free_energy_kcal: float
    ensemble_diversity: float
    frequency_mfe_in_ensemble: float
    pair_probabilities: dict = field(default_factory=dict, repr=False)

    def to_dict(self, *, with_matrix: bool = False) -> dict:
        d = {
            "sequence_length": self.sequence_length,
            "temperature_c": self.temperature_c,
            "circular": self.circular,
            "mfe_structure": self.mfe_structure,
            "mfe_kcal": round(self.mfe_kcal, 3),
            "mfe_per_nt": round(self.mfe_per_nt, 5),
            "ensemble_free_energy_kcal": round(self.ensemble_free_energy_kcal, 3),
            "ensemble_diversity": round(self.ensemble_diversity, 4),
            "frequency_mfe_in_ensemble": round(self.frequency_mfe_in_ensemble, 6),
            "vienna_version": VIENNA_VERSION,
        }
        if with_matrix:
            d["pair_probabilities"] = self.pair_probabilities
        return d


def _compound(seq: str, *, circular: bool, temperature_c: float) -> RNA.fold_compound:
    md = RNA.md()
    md.circ = 1 if circular else 0
    md.temperature = float(temperature_c)
    md.uniq_ML = 1
    return RNA.fold_compound(to_rna(seq), md)


def fold(seq: str, *, circular: bool = True, temperature_c: float = 37.0) -> FoldResult:
    """Lipat sekuens; `circular=True` setara `RNAfold --circ`.

    Mengembalikan MFE, energi bebas ensemble, keragaman ensemble (mean base pair
    distance), frekuensi struktur MFE dalam ensemble, dan matriks probabilitas
    pasangan basa.
    """
    rna = to_rna(seq)
    fc = _compound(rna, circular=circular, temperature_c=temperature_c)
    structure, mfe = fc.mfe()
    # Partition function harus dihitung relatif terhadap MFE agar stabil numerik.
    fc.exp_params_rescale(mfe)
    _, efe = fc.pf()
    diversity = fc.mean_bp_distance()
    freq = fc.pr_structure(structure)
    bpp = _bpp_dict(fc, len(rna))
    return FoldResult(
        sequence_length=len(rna),
        temperature_c=float(temperature_c),
        circular=circular,
        mfe_structure=structure,
        mfe_kcal=float(mfe),
        mfe_per_nt=float(mfe) / len(rna),
        ensemble_free_energy_kcal=float(efe),
        ensemble_diversity=float(diversity),
        frequency_mfe_in_ensemble=float(freq),
        pair_probabilities=bpp,
    )


def _bpp_dict(fc: RNA.fold_compound, n: int, cutoff: float = 1e-4) -> dict:
    """Probabilitas pasangan basa sebagai {(i,j) 1-based: p} di atas cutoff."""
    matrix = fc.bpp()
    out = {}
    for i in range(1, n + 1):
        row = matrix[i]
        for j in range(i + 1, n + 1):
            p = row[j]
            if p > cutoff:
                out[(i, j)] = float(p)
    return out


def unpaired_probability(result: FoldResult, positions: list[int]) -> dict[int, float]:
    """P(tak-berpasangan) tiap posisi (1-based) dari matriks probabilitas pasangan."""
    paired: dict[int, float] = {p: 0.0 for p in positions}
    wanted = set(positions)
    for (i, j), p in result.pair_probabilities.items():
        if i in wanted:
            paired[i] += p
        if j in wanted:
            paired[j] += p
    return {p: max(0.0, 1.0 - v) for p, v in paired.items()}


def junction_pairing(result: FoldResult, junction_pos: int, span: int = 10) -> dict:
    """Probabilitas pasangan basa di jendela sekitar sambungan balik (BSJ).

    `junction_pos` adalah posisi 1-based nukleotida terakhir E2 (yaitu basa di
    sisi 5' sambungan) pada lingkaran. Jendela mencakup span nt di tiap sisi.
    """
    n = result.sequence_length
    window = [((junction_pos - span - 1 + k) % n) + 1 for k in range(2 * span)]
    unpaired = unpaired_probability(result, window)
    mean_unpaired = sum(unpaired.values()) / len(unpaired)
    return {
        "junction_position": junction_pos,
        "window_span_each_side": span,
        "positions": window,
        "unpaired_probability": {k: round(v, 5) for k, v in unpaired.items()},
        "mean_unpaired_probability": round(mean_unpaired, 5),
        "mean_paired_probability": round(1.0 - mean_unpaired, 5),
        "structure_window": "".join(result.mfe_structure[p - 1] for p in window),
    }


def fold_profile(
    seq: str, temperatures: tuple[float, ...] = (37.0, 55.0), *, circular: bool = True
) -> dict[str, FoldResult]:
    """Lipat pada beberapa suhu. 37 C = fisiologis, 55 C = kondisi splicing PIE in vitro."""
    return {f"{t:g}C": fold(seq, circular=circular, temperature_c=t) for t in temperatures}


def subsequence_structure_identity(
    full_structure: str, sub_structure: str, offset: int
) -> float:
    """Fraksi posisi dengan status lipatan (berpasangan/tidak) yang sama.

    Dipakai untuk menilai apakah domain IRES mempertahankan lipatan seperti saat
    dilipat sendiri (`sub_structure`) ketika berada dalam konstruk penuh.
    """
    if not sub_structure:
        return 0.0
    same = 0
    for k, ch in enumerate(sub_structure):
        pos = offset + k
        if pos >= len(full_structure):
            break
        a = "." if full_structure[pos] == "." else "p"
        b = "." if ch == "." else "p"
        same += a == b
    return same / len(sub_structure)
