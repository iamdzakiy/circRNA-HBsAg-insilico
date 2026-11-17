"""Tests for the single-source-of-truth layer and its generators.

The point of these tests is that a deliverable cannot quietly disagree with the
content: the validator must reject a dangling citation, and the build must put
the design-lock numbers into the generated files rather than a copy of them.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from circhbsag import content as C

REPO = Path(__file__).resolve().parents[1]
LOCK = REPO / "results" / "demo" / "design_lock.json"


@pytest.fixture(scope="module")
def loaded():
    return C.load(REPO / "content", design_lock=LOCK if LOCK.exists() else None)


def test_content_loads_and_validates(loaded):
    assert loaded.references and loaded.decisions and loaded.software
    assert C.validate(loaded) == []


def test_every_reference_has_a_real_doi_and_pmid(loaded):
    for key, r in loaded.references.items():
        assert str(r["doi"]).startswith("10."), key
        assert str(r.get("pmid", "")).isdigit(), f"{key} tidak punya PMID"


def test_literature_decisions_must_cite_something(loaded):
    for d in loaded.decisions:
        if d["label"] == "LITERATUR":
            assert d["refs"], d["id"]
            for k in d["refs"]:
                assert k in loaded.references


def test_dangling_citation_is_rejected(tmp_path):
    src = REPO / "content"
    dst = tmp_path / "content"
    shutil.copytree(src, dst)
    data = yaml.safe_load((dst / "decisions.yaml").read_text())
    data["decisions"][0]["refs"] = ["tidak_ada_2099"]
    (dst / "decisions.yaml").write_text(yaml.safe_dump(data, allow_unicode=True))
    with pytest.raises(C.ContentError) as exc:
        C.load(dst)
    assert "tidak_ada_2099" in str(exc.value)


def test_unknown_provenance_label_is_rejected(tmp_path):
    dst = tmp_path / "content"
    shutil.copytree(REPO / "content", dst)
    data = yaml.safe_load((dst / "decisions.yaml").read_text())
    data["decisions"][0]["label"] = "KEREN"
    (dst / "decisions.yaml").write_text(yaml.safe_dump(data, allow_unicode=True))
    with pytest.raises(C.ContentError) as exc:
        C.load(dst)
    assert "label tidak dikenal" in str(exc.value)


def test_software_without_a_licence_is_rejected(tmp_path):
    dst = tmp_path / "content"
    shutil.copytree(REPO / "content", dst)
    data = yaml.safe_load((dst / "software.yaml").read_text())
    data["software"][0]["license"] = ""
    (dst / "software.yaml").write_text(yaml.safe_dump(data, allow_unicode=True))
    with pytest.raises(C.ContentError) as exc:
        C.load(dst)
    assert "license" in str(exc.value)


def test_missing_numbers_render_as_a_visible_placeholder():
    c = C.load(REPO / "content")  # no design lock
    assert c.has_numbers is False
    assert c.num("lengths_nt", "circle") == "[BELUM DIHITUNG]"
    assert "[RENCANA]" in c.data_label


def test_mandatory_tables_are_non_empty(loaded):
    assert len(C.assumptions_rows(loaded)) >= 5
    assert len(C.open_questions_rows(loaded)) >= 5
    assert len(C.software_rows(loaded)) >= 10
    for row in C.software_rows(loaded):
        assert row[3] and row[4]  # licence and URL


def test_versions_are_detected_not_written_by_hand():
    v = C.detected_versions()
    assert "ViennaRNA" in v
    assert v["ViennaRNA"][0].isdigit()


# ------------------------------------------------------------------ builders


def _build(target: str, cwd: Path = REPO) -> subprocess.CompletedProcess:
    env = {
        "PYTHONPATH": str(REPO / "src"),
        "PATH": "/usr/bin:/bin:/usr/local/bin",
        "HOME": str(cwd),
    }
    return subprocess.run(
        [sys.executable, "-m", "tools.build", target],
        cwd=REPO,
        env=env,
        capture_output=True,
        text=True,
        timeout=900,
    )


@pytest.mark.skipif(not LOCK.exists(), reason="butuh design_lock.json dari jalankan pipeline")
def test_site_pages_exist_and_quote_the_design_lock():
    lock = json.loads(LOCK.read_text())
    circle = str(lock["lengths_nt"]["circle"])
    for name in ("index.html", "gate0.html", "design.html", "methods.html", "decisions.html"):
        p = REPO / "build" / "site" / name
        assert p.exists(), name
        text = p.read_text()
        assert "DIHASILKAN OLEH tools/build.py" in text
    assert circle in (REPO / "build" / "site" / "index.html").read_text()
    assert circle in (REPO / "build" / "site" / "design.html").read_text()


@pytest.mark.skipif(not LOCK.exists(), reason="butuh design_lock.json")
def test_demo_label_reaches_every_generated_surface():
    for name in ("index.html", "gate0.html", "design.html"):
        assert "DEMO — SEKUENS PLACEHOLDER" in (REPO / "build" / "site" / name).read_text()
    assert "Proyeksi berbasis literatur" in (REPO / "build" / "deck_script.md").read_text()


@pytest.mark.skipif(not LOCK.exists(), reason="butuh design_lock.json")
def test_build_check_detects_a_hand_edited_page():
    page = REPO / "build" / "site" / "index.html"
    original = page.read_text()
    try:
        page.write_text(original.replace("</main>", "<p>disunting tangan</p></main>"))
        proc = _build("check")
        assert proc.returncode == 1
        assert "kedaluwarsa" in proc.stdout
    finally:
        page.write_text(original)
    assert _build("check").returncode == 0


@pytest.mark.skipif(not LOCK.exists(), reason="butuh design_lock.json")
def test_generated_site_has_no_external_resources():
    """A GitHub Pages site that works offline is one fewer thing that can rot."""
    for p in (REPO / "build" / "site").glob("*.html"):
        text = p.read_text()
        for bad in ("<script src", "cdn.", "unpkg", "googleapis"):
            assert bad not in text, f"{p.name} memuat sumber eksternal: {bad}"


@pytest.mark.skipif(not LOCK.exists(), reason="butuh design_lock.json")
def test_data_target_publishes_numbers_for_the_public_site():
    """docs/ is written by hand elsewhere; this target hands it JSON, not HTML."""
    d = REPO / "docs" / "data"
    lock = json.loads((d / "design_lock.json").read_text())
    assert lock["lengths_nt"]["circle"] > 0
    cj = json.loads((d / "content.json").read_text())
    assert cj["references"] and cj["decisions"] and cj["software"]
    assert cj["projection_label"] == C.PROJECTION_LABEL
    assert not list(d.glob("*.html")), "target data tidak boleh menulis HTML"
