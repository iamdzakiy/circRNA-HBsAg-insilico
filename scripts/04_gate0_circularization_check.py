#!/usr/bin/env python3
"""Run the nine Gate 0 tests; write reports/gate0.{json,md} and design_lock.json."""
from __future__ import annotations

from _common import parser, sidecar

from circhbsag import cli


def main() -> int:
    p = parser("Gate 0: nine tests")
    p.add_argument("--cm", default="data/reference/RF00028.cm", help="Rfam RF00028 covariance model")
    p.add_argument("--temperatures", default="37,55")
    a = p.parse_args()
    out = f"{a.results}/reports"
    argv = ["gate0", "--parts", a.parts, "--out", out, "--seed", str(a.seed), "--cm", a.cm,
            "--temperatures", a.temperatures]
    if a.orf_fasta:
        argv += ["--orf-fasta", a.orf_fasta]
    rc = cli.main(argv)
    sidecar(__import__("pathlib").Path(out) / "gate0", a)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
