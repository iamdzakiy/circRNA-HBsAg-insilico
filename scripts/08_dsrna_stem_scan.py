#!/usr/bin/env python3
"""Longest contiguous double-stranded stretch in the circle MFE (Gate 0 Uji 7)."""
from __future__ import annotations

import json

from _common import out_dir, parser, sidecar
from _gate_view import test_metrics


def main() -> int:
    a = parser("dsRNA stem scan", orf=False).parse_args()
    t = test_metrics(a.results, "uji7")
    out = out_dir(a, "tables") / "dsrna_stems.json"
    out.write_text(json.dumps(t, indent=2, ensure_ascii=False))
    print(f"Uji 7: {t['verdict']}  longest contiguous bp={t['metrics']['longest_contiguous_bp']}")
    sidecar(out, a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
