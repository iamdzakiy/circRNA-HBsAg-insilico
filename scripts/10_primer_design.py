#!/usr/bin/env python3
"""Primers and probes: divergent PCR, BSJ-spanning, IVT reverse, RNase H (primer3)."""
from __future__ import annotations

import json

from _common import load_parts, out_dir, parser, sidecar

from circhbsag import construct, primers


def main() -> int:
    a = parser("Primer design").parse_args()
    parts, _ = load_parts(a)
    asm = construct.build(parts)
    res = primers.design_all(asm)
    out = out_dir(a, "tables") / "primers.json"
    out.write_text(json.dumps(res, indent=2, ensure_ascii=False, default=str))
    print(f"primer sets: {', '.join(res)} -> {out}")
    sidecar(out, a, {"label": "[HASIL IN SILICO] belum divalidasi eksperimen"})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
