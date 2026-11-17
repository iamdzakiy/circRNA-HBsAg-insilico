#!/usr/bin/env python3
"""List RRACH (m6A consensus) sites in the circle and the ORF."""
from __future__ import annotations

import csv

from _common import load_parts, out_dir, parser, sidecar

from circhbsag import construct, motifs


def main() -> int:
    a = parser("m6A RRACH scan").parse_args()
    parts, _ = load_parts(a)
    asm = construct.build(parts)
    t = out_dir(a, "tables") / "m6a_rrach_sites.tsv"
    with t.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["scope", "start_0based", "motif", "a_position_0based"])
        for scope, seq in (("circle", asm.circle), ("ORF", parts["ORF"].seq)):
            for s in motifs.scan_m6a_rrach(seq):
                w.writerow([scope, s["start"], s["motif"], s["a_position"]])
    print(f"circle density {motifs.m6a_density(asm.circle):.2f}/100 nt -> {t}")
    sidecar(t, a, {"label": "[LITERATUR] RRACH konsensus; situs motif bukan bukti metilasi"})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
