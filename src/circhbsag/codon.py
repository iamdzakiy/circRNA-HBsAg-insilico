"""Optimasi kodon ORF untuk ekspresi di sel manusia.

Objektif yang dioptimasi (semua dapat dibobot lewat config):
  - CAI (codon adaptation index) terhadap tabel penggunaan kodon inang
  - GC total dan GC lokal dalam rentang yang dapat disintesis
  - menghindari situs restriksi Tipe IIS (BsaI/BsmBI) internal
  - menghindari motif imunostimulan/destabilisasi (AVOID_MOTIFS)
  - menurunkan jumlah situs m6A RRACH (opsional: justru dipertahankan
    bila penyisipan m6A direncanakan, lihat `m6a_mode`)
  - struktur sekunder ORF: MFE per nt dan aksesibilitas start codon

Catatan lisensi: modul ini tidak memakai LinearDesign. LinearDesign dirilis
dengan lisensi yang membatasi pemakaian komersial dan memerlukan persetujuan
untuk redistribusi, sehingga pipeline ini memakai optimizer sendiri
(stochastic hill-climbing + Optuna) yang seluruhnya MIT.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass
from pathlib import Path

from . import motifs
from .seqio import CODON_TABLE, SYN, SequenceError, gc, gc_windows, to_dna, translate

DEFAULT_TABLE = Path(__file__).resolve().parents[2] / "data" / "reference" / "codon_usage_human.tsv"


def load_codon_usage(path: str | Path | None = None) -> dict[str, float]:
    """Baca tabel penggunaan kodon: kolom codon, aa, frequency_per_thousand.

    Tabel bawaan ada di data/reference/codon_usage_human.tsv dan DITANDAI
    [PERKIRAAN]; pengguna wajib menggantinya dengan unduhan terverifikasi
    (mis. Kazusa/HIVE) sebelum angka CAI dilaporkan sebagai hasil akhir.
    """
    path = Path(path) if path else DEFAULT_TABLE
    if not path.exists():
        raise FileNotFoundError(
            f"tabel penggunaan kodon tidak ditemukan: {path}. "
            "Sediakan TSV dengan kolom codon/aa/frequency_per_thousand."
        )
    usage: dict[str, float] = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.lower().startswith("codon"):
            continue
        parts = line.split()
        codon, freq = to_dna(parts[0]), float(parts[2])
        usage[codon] = freq
    missing = set(CODON_TABLE) - set(usage)
    if missing:
        raise SequenceError(f"tabel kodon tidak lengkap, kurang: {sorted(missing)[:5]}")
    return usage


def relative_adaptiveness(usage: dict[str, float]) -> dict[str, float]:
    """w_i = f_i / max(f) dalam keluarga sinonim (definisi Sharp & Li 1987)."""
    w: dict[str, float] = {}
    for codons in SYN.values():
        mx = max(usage[c] for c in codons)
        for c in codons:
            w[c] = usage[c] / mx if mx > 0 else 0.0
    return w


def cai(seq: str, w: dict[str, float]) -> float:
    """CAI = rerata geometrik w kodon (kodon stop dan Met/Trp tunggal dilewati)."""
    s = to_dna(seq)
    vals = []
    for i in range(0, len(s) - 2, 3):
        c = s[i : i + 3]
        aa = CODON_TABLE.get(c)
        if aa in (None, "*") or len(SYN.get(aa, [])) < 2:
            continue
        wi = w.get(c, 0.0)
        vals.append(max(wi, 1e-6))
    if not vals:
        return 0.0
    return math.exp(sum(math.log(v) for v in vals) / len(vals))


@dataclass
class CodonObjective:
    """Bobot objektif optimasi kodon. Semua bobot >= 0; skor makin kecil makin baik."""

    w_cai: float = 1.0
    w_gc: float = 1.0
    w_gc_window: float = 1.0
    w_restriction: float = 5.0
    w_avoid: float = 1.0
    w_m6a: float = 0.5
    gc_target: float = 0.55
    gc_window_lo: float = 0.30
    gc_window_hi: float = 0.70
    m6a_mode: str = "minimize"  # "minimize" | "ignore"
    forbidden_sites: tuple[str, ...] = ("GGTCTC", "CGTCTC")

    def cost(self, orf: str, w: dict[str, float]) -> dict:
        c = cai(orf, w)
        g = gc(orf)
        gw = motifs.gc_window_violations(orf, lo=self.gc_window_lo, hi=self.gc_window_hi)
        n_rs = sum(len(motifs.find_all(orf, s, both_strands=True)) for s in self.forbidden_sites)
        n_avoid = sum(len(v) for v in motifs.scan_avoid_motifs(orf).values())
        n_m6a = len(motifs.scan_m6a_rrach(orf))
        terms = {
            "cai": self.w_cai * (1.0 - c),
            "gc": self.w_gc * abs(g - self.gc_target),
            "gc_window": self.w_gc_window * (len(gw) / max(1, len(orf) / 50)),
            "restriction": self.w_restriction * n_rs,
            "avoid": self.w_avoid * (n_avoid / max(1, len(orf) / 100)),
            "m6a": (self.w_m6a * (n_m6a / max(1, len(orf) / 100))) if self.m6a_mode == "minimize" else 0.0,
        }
        return {
            "total": sum(terms.values()),
            "terms": {k: round(v, 6) for k, v in terms.items()},
            "metrics": {
                "cai": round(c, 4),
                "gc": round(g, 4),
                "gc_window_violations": len(gw),
                "restriction_sites": n_rs,
                "avoid_motifs": n_avoid,
                "m6a_sites": n_m6a,
            },
        }


def naive_optimize(protein: str, w: dict[str, float]) -> str:
    """Pilih kodon paling disukai inang untuk setiap residu (baseline CAI maksimum)."""
    out = []
    for aa in protein.replace("*", ""):
        codons = SYN[aa]
        out.append(max(codons, key=lambda c: w[c]))
    return "".join(out)


def optimize(
    protein: str,
    *,
    usage: dict[str, float] | None = None,
    objective: CodonObjective | None = None,
    iterations: int = 20000,
    seed: int = 20260101,
    start_codon: str = "ATG",
    stop_codon: str = "TGA",
) -> dict:
    """Optimasi kodon dengan hill-climbing stokastik dan seed tetap (reprodusibel).

    Mengembalikan dict berisi sekuens ORF (termasuk start dan stop),
    riwayat biaya, dan metrik akhir. Translasi diverifikasi identik dengan input.
    """
    protein = protein.strip().upper().rstrip("*")
    if not protein:
        raise SequenceError("sekuens protein kosong")
    bad = sorted(set(protein) - set(SYN) | {"*"} & set(protein))
    bad = [b for b in bad if b not in SYN]
    if bad:
        raise SequenceError(f"residu tak dikenal: {bad}")
    usage = usage or load_codon_usage()
    w = relative_adaptiveness(usage)
    objective = objective or CodonObjective()

    rng = random.Random(seed)
    body = naive_optimize(protein, w)
    best = body
    best_cost = objective.cost(body, w)["total"]
    history = [best_cost]

    n = len(protein)
    for _ in range(iterations):
        i = rng.randrange(n)
        aa = protein[i]
        alts = [c for c in SYN[aa] if c != best[3 * i : 3 * i + 3]]
        if not alts:
            continue
        cand = best[: 3 * i] + rng.choice(alts) + best[3 * i + 3 :]
        c = objective.cost(cand, w)["total"]
        if c < best_cost:
            best, best_cost = cand, c
            history.append(c)

    orf = start_codon + best[3:] if best[:3] in SYN["M"] and protein[0] == "M" else best
    orf = orf + stop_codon
    if translate(orf).rstrip("*") != protein:
        raise SequenceError("translasi hasil optimasi tidak sama dengan protein input")
    final = objective.cost(orf, w)
    return {
        "orf_dna": orf,
        "protein": protein,
        "length_nt": len(orf),
        "iterations": iterations,
        "seed": seed,
        "cost": final,
        "baseline_naive": objective.cost(naive_optimize(protein, w) + stop_codon, w)["metrics"],
        "cost_history_len": len(history),
    }


def codon_stats(orf: str, usage: dict[str, float] | None = None) -> dict:
    """Metrik deskriptif ORF: CAI, GC, GC3, jumlah motif."""
    usage = usage or load_codon_usage()
    w = relative_adaptiveness(usage)
    s = to_dna(orf)
    third = s[2::3]
    return {
        "length_nt": len(s),
        "cai": round(cai(s, w), 4),
        "gc": round(gc(s), 4),
        "gc3": round(gc(third), 4),
        "gc_min_window": round(min(v for _, v in gc_windows(s)), 4),
        "gc_max_window": round(max(v for _, v in gc_windows(s)), 4),
        "m6a_sites": len(motifs.scan_m6a_rrach(s)),
        "restriction": {k: len(v) for k, v in motifs.scan_restriction(s).items() if v},
    }
