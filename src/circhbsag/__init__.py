"""circhbsag — pipeline in silico untuk konstruk circRNA-HBsAg berbasis PIE.

Paket ini menyediakan Tahap 0 (keputusan GO/NO-GO sebelum masuk laboratorium)
dan Tahap 1 (desain serta optimasi in silico) untuk konstruk:

    intron 3' PIE | E2 | spacer | IRES | ORF HBsAg | spacer | E1 | intron 5' PIE

Prinsip yang dipegang seluruh paket:

  * Tidak ada sekuens yang dikarang. Semua sekuens referensi datang dari berkas
    yang disediakan pengguna; bila tidak ada, proses berhenti dengan galat.
  * Tidak ada angka hasil yang dihardcode. Setiap nilai numerik yang dilaporkan
    berasal dari perhitungan saat dijalankan (ViennaRNA, primer3, atau kode ini).
  * Setiap ambang dan bobot didokumentasikan dan dapat diubah lewat config.
  * Setiap jalankan memakai seed tetap dan mencatat versi alat.

Lisensi paket: MIT. Lisensi dependensi dicantumkan di README.md dan dapat
diperiksa saat runtime dengan `circhbsag.licenses()`.
"""
from __future__ import annotations

__version__ = "0.1.0"

from . import (  # noqa: F401
    access,
    bsj_stats,
    codon,
    construct,
    fold,
    motifs,
    offtarget,
    optimize,
    primers,
    score,
    seqio,
)


def licenses() -> dict[str, dict[str, str]]:
    """Daftar alat yang dipakai beserta lisensinya (diverifikasi saat runtime)."""
    import importlib.metadata as md

    def ver(pkg: str) -> str:
        try:
            return md.version(pkg)
        except Exception:  # pragma: no cover
            return "tidak terpasang"

    return {
        "circhbsag": {"version": __version__, "license": "MIT", "role": "pipeline ini"},
        "ViennaRNA": {
            "version": ver("ViennaRNA"),
            "license": "ViennaRNA License (gratis, non-komersial/akademik)",
            "role": "RNAfold --circ, RNAplfold, partition function",
        },
        "primer3-py": {
            "version": ver("primer3-py"),
            "license": "GPLv2",
            "role": "desain primer divergen/BSJ, Tm, hairpin",
        },
        "biopython": {
            "version": ver("biopython"),
            "license": "Biopython License (mirip MIT)",
            "role": "utilitas sekuens dan I/O",
        },
        "optuna": {
            "version": ver("optuna"),
            "license": "MIT",
            "role": "optimasi multi-objektif NSGA-II",
        },
        "snakemake": {
            "version": ver("snakemake"),
            "license": "MIT",
            "role": "orkestrasi alur kerja reprodusibel",
        },
        "numpy": {"version": ver("numpy"), "license": "BSD-3-Clause", "role": "numerik"},
        "pandas": {"version": ver("pandas"), "license": "BSD-3-Clause", "role": "tabel hasil"},
        "BLAST+": {
            "version": "lihat `blastn -version`",
            "license": "domain publik (NCBI)",
            "role": "homologi off-target (opsional)",
        },
        "STAR": {
            "version": "opsional",
            "license": "MIT",
            "role": "validasi silang BSJ (opsional)",
        },
    }
