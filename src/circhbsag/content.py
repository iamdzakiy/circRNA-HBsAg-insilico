"""Loader and validator for ``content/`` — the project's single source of truth.

Every deliverable (website, thesis document, slide script, README tables) is
generated from the YAML files in ``content/`` plus ``design_lock.json``. Nothing
downstream is allowed to hold its own copy of a number or a citation, so that a
disagreement between deliverables becomes a stale generated file — which CI can
detect — rather than an inconsistency nobody can trace.

The validator is the part that gives this teeth. It refuses content that cites a
reference key that does not exist, uses a provenance label that is not declared,
or lists a tool without a licence. A build that would produce a document with a
fabricated-looking citation fails instead.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

#: Files expected in ``content/``, mapped to the top-level key each one carries.
FILES = {
    "project.yaml": None,
    "references.yaml": "references",
    "decisions.yaml": "decisions",
    "software.yaml": "software",
    "assumptions.yaml": "assumptions",
    "open_questions.yaml": "open_questions",
}

PROJECTION_LABEL = "Proyeksi berbasis literatur — bukan data eksperimen"


class ContentError(ValueError):
    """Raised when content is internally inconsistent."""


@dataclass
class Content:
    """All project content, plus the design lock when a run has produced one."""

    root: Path
    project: dict
    references: dict
    decisions: list
    software: list
    software_not_used: list
    assumptions: list
    open_questions: list
    design_lock: dict | None = None
    raw: dict[str, Any] = None

    # -- convenience -----------------------------------------------------
    @property
    def labels(self) -> dict:
        return self.project.get("labels", {})

    @property
    def has_numbers(self) -> bool:
        return bool(self.design_lock)

    @property
    def data_label(self) -> str:
        if not self.design_lock:
            return (
                "[RENCANA] belum ada jalankan pipeline, sehingga dokumen ini belum "
                "memuat angka konstruk."
            )
        return self.design_lock.get("data_label", "")

    @property
    def is_demo(self) -> bool:
        return bool(self.design_lock and self.design_lock.get("data_status") == "DEMO-PLACEHOLDER")

    def ref(self, key: str) -> dict:
        try:
            return self.references[key]
        except KeyError as exc:  # pragma: no cover - validator catches this first
            raise ContentError(f"referensi tidak dikenal: {key}") from exc

    def citation(self, key: str) -> str:
        r = self.ref(key)
        return (
            f"{r['authors']}. {r['title']}. *{r['journal']}* {r['year']}; "
            f"{r['citation']}. https://doi.org/{r['doi']}"
        )

    def short_cite(self, key: str) -> str:
        r = self.ref(key)
        first = r["authors"].split(",")[0].split(" ")[0]
        return f"{first} dkk. {r['year']}"

    def lock(self, *path: str, default: Any = None) -> Any:
        """Read a nested value out of the design lock, or ``default``.

        Deliberately forgiving: a deliverable that wants a number which no run
        has produced yet renders a visible placeholder instead of crashing, so
        the document can be built before the pipeline has real sequences.
        """
        cur: Any = self.design_lock or {}
        for p in path:
            if not isinstance(cur, dict) or p not in cur:
                return default
            cur = cur[p]
        return cur

    def num(self, *path: str, unit: str = "", nd: int | None = None) -> str:
        """Format a design-lock number for prose, or a visible placeholder."""
        v = self.lock(*path)
        if v is None:
            return "[BELUM DIHITUNG]"
        if isinstance(v, float) and nd is not None:
            s = f"{v:.{nd}f}"
        else:
            s = f"{v}"
        return f"{s} {unit}".strip()


def load(root: str | Path = "content", *, design_lock: str | Path | None = None) -> Content:
    """Load and validate every content file. Raises :class:`ContentError`."""
    root = Path(root)
    raw: dict[str, Any] = {}
    for name in FILES:
        path = root / name
        if not path.exists():
            raise ContentError(f"berkas konten tidak ada: {path}")
        raw[name] = yaml.safe_load(path.read_text()) or {}

    lock = None
    if design_lock:
        p = Path(design_lock)
        if p.exists():
            lock = json.loads(p.read_text())

    c = Content(
        root=root,
        project=raw["project.yaml"],
        references=raw["references.yaml"]["references"],
        decisions=raw["decisions.yaml"]["decisions"],
        software=raw["software.yaml"]["software"],
        software_not_used=raw["software.yaml"].get("not_used_on_purpose", []),
        assumptions=raw["assumptions.yaml"]["assumptions"],
        open_questions=raw["open_questions.yaml"]["open_questions"],
        design_lock=lock,
        raw=raw,
    )
    errors = validate(c)
    if errors:
        raise ContentError("konten tidak konsisten:\n  - " + "\n  - ".join(errors))
    return c


def validate(c: Content) -> list[str]:
    """Return a list of content problems; empty means the content is coherent."""
    errors: list[str] = []
    labels = set(c.labels)

    if not labels:
        errors.append("project.yaml tidak mendeklarasikan daftar label provenans")

    # References must be complete enough to be checkable by a reader.
    for key, r in c.references.items():
        for field in ("authors", "title", "journal", "year", "citation", "doi"):
            if not r.get(field):
                errors.append(f"referensi {key} tidak punya {field}")
        if r.get("doi") and not str(r["doi"]).startswith("10."):
            errors.append(f"referensi {key}: DOI tidak berbentuk DOI ({r['doi']})")

    # Decisions must point at labels and references that exist.
    seen_ids = set()
    for d in c.decisions:
        if d["id"] in seen_ids:
            errors.append(f"id keputusan ganda: {d['id']}")
        seen_ids.add(d["id"])
        if d.get("label") not in labels:
            errors.append(f"keputusan {d['id']}: label tidak dikenal ({d.get('label')})")
        for key in d.get("refs", []):
            if key not in c.references:
                errors.append(f"keputusan {d['id']}: referensi tidak ada ({key})")
        if d.get("label") == "LITERATUR" and not d.get("refs"):
            errors.append(
                f"keputusan {d['id']} berlabel LITERATUR tetapi tidak menyebut referensi"
            )
        if not d.get("rationale"):
            errors.append(f"keputusan {d['id']}: tidak ada alasan")

    # Every tool must carry a licence and a URL — the point of the table.
    for s in c.software:
        for field in ("name", "role", "license", "url"):
            if not s.get(field):
                errors.append(f"software {s.get('name', '?')}: tidak punya {field}")
    for s in c.software_not_used:
        if not s.get("reason"):
            errors.append(f"software ditolak {s.get('name', '?')}: tidak ada alasan")

    for a in c.assumptions:
        for field in ("assumption", "why", "impact_if_wrong", "verify"):
            if not a.get(field):
                errors.append(f"asumsi {a.get('id', '?')}: tidak punya {field}")
    for q in c.open_questions:
        for field in ("question", "answered_by"):
            if not q.get(field):
                errors.append(f"pertanyaan {q.get('id', '?')}: tidak punya {field}")

    return errors


# --------------------------------------------------------------------------
# Shared renderings of the three mandatory tables
# --------------------------------------------------------------------------


def assumptions_rows(c: Content) -> list[list[str]]:
    return [
        [a["id"], a["assumption"], a["why"], a["impact_if_wrong"], a["verify"]]
        for a in c.assumptions
    ]


def open_questions_rows(c: Content) -> list[list[str]]:
    return [
        [q["id"], q["question"], q["answered_by"], q.get("blocking", "—")]
        for q in c.open_questions
    ]


def software_rows(c: Content) -> list[list[str]]:
    rows = []
    for s in c.software:
        rows.append(
            [
                s["name"],
                s["role"],
                str(s.get("version", "—")),
                s["license"],
                s["url"],
                "dipakai" if s.get("used") else "belum dipakai",
            ]
        )
    return rows


def detected_versions() -> dict[str, str]:
    """Versions read from the installed packages, never written by hand."""
    out: dict[str, str] = {}
    try:
        import RNA

        out["ViennaRNA"] = RNA.__version__
    except Exception:  # pragma: no cover
        pass
    for mod, label in (
        ("primer3", "primer3 / primer3-py"),
        ("Bio", "Biopython"),
        ("optuna", "Optuna"),
        ("snakemake", "Snakemake"),
        ("numpy", "NumPy / pandas / SciPy"),
        ("streamlit", "Streamlit"),
        ("docx", "python-docx"),
        ("pptx", "python-pptx"),
    ):
        try:
            m = __import__(mod)
            v = getattr(m, "__version__", None)
            if v:
                out[label] = str(v)
        except Exception:
            continue
    return out
