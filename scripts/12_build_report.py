#!/usr/bin/env python3
"""Assemble the Gate 0 Report (Markdown) from results/ + the three mandatory tables."""
from __future__ import annotations

from pathlib import Path

from _common import ROOT, need, out_dir, parser, sidecar

from circhbsag import content


def table(head: list[str], rows: list[list[str]]) -> str:
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    out += ["| " + " | ".join(str(c).replace("\n", " ") for c in r) + " |" for r in rows]
    return "\n".join(out) + "\n"


def main() -> int:
    a = parser("Build Gate 0 report", orf=False).parse_args()
    res = Path(a.results)
    g0 = need(res / "reports" / "gate0.md", "run scripts/04_gate0_circularization_check.py first").read_text()
    it = res / "reports" / "gate0_iterations.md"
    if not it.exists():
        it = res / "demo" / "gate0_iterations.md"
    c = content.load(ROOT / "content", design_lock=res / "reports" / "design_lock.json")
    fig = [f for f in ("dG_vs_arm_length.png", "gc_profile.png") if (res / "figures" / f).exists()]
    parts = [
        "# Gate 0 Report — circRNA-HBsAg (PIE, IRES CVB3)\n",
        f"> {content.PROJECTION_LABEL}\n> {c.data_label}\n",
        g0,
        "\n## Redesign iterations (one parameter per iteration)\n",
        it.read_text() if it.exists() else "_No redesign loop run yet (python -m circhbsag.cli redesign)._\n",
        "\n## Figures\n",
        "\n".join(f"![{f}](../figures/{f})" for f in fig) or "_run scripts/11_make_figures.py_",
        "\n\n## Assumption Register\n",
        table(["ID", "Assumption", "Why", "Impact if wrong", "Verify"], content.assumptions_rows(c)),
        "\n## Open Questions\n",
        table(["ID", "Question", "Answered by", "Blocking"], content.open_questions_rows(c)),
        "\n## Software & Versions\n",
        table(["Tool", "Role", "Version", "Licence", "URL", "Status"], content.software_rows(c)),
        "\n## AI-use statement\nAI (Claude) was used to draft code, tests and text. Every number in this report was "
        "computed by the code in this repository in the run recorded in the sidecar files; literature values carry "
        "full citations in content/references.yaml. No experimental data exist yet.\n",
    ]
    out = out_dir(a, "reports") / "gate0_report.md"
    out.write_text("\n".join(parts))
    print(out)
    sidecar(out, a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
