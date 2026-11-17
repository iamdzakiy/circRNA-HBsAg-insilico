"""Small helper: read one Gate 0 test out of results/reports/gate0.json."""
from __future__ import annotations

import json
from pathlib import Path

from _common import need


def test_metrics(results: str, test_prefix: str) -> dict:
    p = need(Path(results) / "reports" / "gate0.json", "run scripts/04_gate0_circularization_check.py first")
    g = json.loads(p.read_text())["results"]["gate0"]
    t = next(x for x in g["tests"] if x["id"].startswith(test_prefix))
    return t
