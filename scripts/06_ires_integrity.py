#!/usr/bin/env python3
"""IRES structure integrity in context vs alone (Gate 0 Uji 6)."""
from __future__ import annotations

import json

from _common import out_dir, parser, sidecar
from _gate_view import test_metrics


def main() -> int:
    a = parser("IRES integrity", orf=False).parse_args()
    t = test_metrics(a.results, "uji6")
    out = out_dir(a, "tables") / "ires_integrity.json"
    out.write_text(json.dumps(t, indent=2, ensure_ascii=False))
    print(f"Uji 6: {t['verdict']}  identity={t['metrics']['structure_identity']}")
    sidecar(out, a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
