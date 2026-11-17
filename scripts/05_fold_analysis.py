#!/usr/bin/env python3
"""MFE, ensemble free energy and diversity of precursor and circle (ViennaRNA)."""
from __future__ import annotations

import csv

from _common import load_parts, out_dir, parser, sidecar

from circhbsag import construct, score


def _flatten(d, prefix=""):
    out = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(_flatten(v, key + "."))
        elif isinstance(v, (int, float, str, bool)):
            out[key] = v
    return out


def main() -> int:
    p = parser("Fold analysis")
    p.add_argument("--temperatures", default="37,55")
    a = p.parse_args()
    parts, _ = load_parts(a)
    asm = construct.build(parts)
    temps = tuple(float(x) for x in a.temperatures.split(","))
    m = score.measure(asm, temperatures=temps)
    t = out_dir(a, "tables") / "fold_analysis.tsv"
    flat = _flatten(m)
    with t.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        for k, v in flat.items():
            w.writerow([k, v])
    print(f"{len(flat)} metrics -> {t}")
    sidecar(t, a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
