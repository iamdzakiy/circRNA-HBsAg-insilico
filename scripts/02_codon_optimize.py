#!/usr/bin/env python3
"""Codon-optimise the antigen ORF (CAI, GC window, avoided motifs, m6A)."""
from __future__ import annotations

from _common import out_dir, parser, sidecar

from circhbsag import cli


def main() -> int:
    p = parser("Codon optimisation", orf=False)
    p.add_argument("--iterations", type=int, default=20000)
    p.add_argument("--m6a-mode", choices=["minimize", "ignore"], default="minimize")
    a = p.parse_args()
    out = out_dir(a, "tables")
    rc = cli.main(["codon", "--parts", a.parts, "--out", str(a.results), "--seed", str(a.seed),
                   "--iterations", str(a.iterations), "--m6a-mode", a.m6a_mode])
    sidecar(out / "codon_optimize", a)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
