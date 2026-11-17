#!/usr/bin/env python3
"""Epitope / HLA coverage — requires NetMHCpan/NetMHCIIpan run LOCALLY.

NetMHCpan-4.1 and NetMHCIIpan-4.0 are free for academic use but cannot be
redistributed or run from this repository. This script therefore never emits
predicted affinities on its own. Given a prediction table you produced locally
(--predictions, TSV with columns peptide, allele, rank_percent), it computes the
fraction of peptides with %Rank <= 2 per allele; without one it writes an
UNAVAILABLE report containing the exact steps.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from _common import out_dir, parser, sidecar

INSTRUCTIONS = (
    "1. Download NetMHCpan-4.1 / NetMHCIIpan-4.0 from DTU Health Tech (academic licence).\n"
    "2. Slide 9-mers (class I) / 15-mers (class II) over the HBsAg protein.\n"
    "3. Save the output as TSV: peptide, allele, rank_percent.\n"
    "4. Re-run with --predictions <file>. Allele panel must come from verified Indonesian frequency data [PERLU VERIFIKASI]."
)


def main() -> int:
    p = parser("Epitope/HLA coverage", orf=False)
    p.add_argument("--predictions", default=None)
    p.add_argument("--rank-cutoff", type=float, default=2.0)
    a = p.parse_args()
    out = out_dir(a, "tables") / "epitope_hla_coverage.json"
    if not a.predictions or not Path(a.predictions).exists():
        res = {"status": "UNAVAILABLE", "reason": "no local NetMHCpan predictions provided", "steps": INSTRUCTIONS}
    else:
        per: dict[str, list[float]] = {}
        for r in csv.DictReader(open(a.predictions), delimiter="\t"):
            per.setdefault(r["allele"], []).append(float(r["rank_percent"]))
        res = {"status": "OK", "rank_cutoff_percent": a.rank_cutoff,
               "per_allele_binder_fraction": {k: sum(x <= a.rank_cutoff for x in v) / len(v) for k, v in per.items()}}
    out.write_text(json.dumps(res, indent=2, ensure_ascii=False))
    print(res["status"], "->", out)
    sidecar(out, a)
    return 0 if res["status"] == "OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())
