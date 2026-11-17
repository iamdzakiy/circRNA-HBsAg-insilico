"""Tes provenans sekuens, penulisan laporan, dan CLI ujung-ke-ujung.

Tes provenans adalah pengaman utama proyek ini: ia memastikan sekuens sintetis
tidak dapat menyusup ke jalur analisis "sungguhan", dan bahwa setiap laporan
membawa label yang benar.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from circhbsag import construct, report, score
from circhbsag import parts as parts_mod
from circhbsag.seqio import SequenceError

REPO = Path(__file__).resolve().parents[1]
EXAMPLE_CFG = REPO / "config" / "parts.example.yaml"
REAL_CFG = REPO / "config" / "parts.yaml"


# ------------------------------------------------------------- provenance


def test_example_config_loads_and_is_flagged_synthetic():
    parts, meta = parts_mod.load_parts(EXAMPLE_CFG)
    assert meta["provenance_policy"] == "allow_synthetic"
    assert meta["synthetic_parts"], "konfigurasi contoh harus menandai bagian sintetis"
    assert meta["results_may_be_reported_as_scientific"] is False
    for name in ("intron_3half", "E2", "IRES", "E1", "intron_5half"):
        assert name in parts and parts[name].seq


def test_real_config_is_strict_and_refuses_missing_references():
    """config/parts.yaml sengaja menunjuk berkas yang belum ada.

    Pipeline harus BERHENTI dengan galat yang menyebut berkas mana yang kurang,
    bukan melanjutkan dengan sekuens pengganti.
    """
    with pytest.raises(SequenceError) as exc:
        parts_mod.load_parts(REAL_CFG)
    msg = str(exc.value)
    assert "data/reference" in msg or "tidak ada" in msg


def test_strict_policy_rejects_synthetic_parts(tmp_path):
    cfg = tmp_path / "p.yaml"
    cfg.write_text(
        "provenance_policy: strict\n"
        "parts:\n"
        "  E1:\n"
        "    inline: ACGTACGTACGT\n"
        "    synthetic: true\n"
    )
    with pytest.raises(parts_mod.ProvenanceError) as exc:
        parts_mod.load_parts(cfg)
    assert "synthetic" in str(exc.value)


def test_part_provenance_records_hash_and_length():
    _parts, meta = parts_mod.load_parts(EXAMPLE_CFG)
    p = meta["parts"]["intron_3half"]
    assert p["length"] == 180
    assert len(p["sha256"]) == 64
    assert p["provenance_status"] == "SINTETIS-UJI"


def test_inline_parts_have_provenance_status():
    _parts, meta = parts_mod.load_parts(EXAMPLE_CFG)
    assert meta["parts"]["T7_promoter"]["provenance_status"] == "DISEDIAKAN-PENGGUNA"


def test_load_protein_from_example():
    protein, pmeta = parts_mod.load_protein(EXAMPLE_CFG)
    assert protein.startswith("M")
    assert pmeta["length_aa"] == 226  # panjang disamakan dengan SHBs
    assert pmeta["synthetic"] is True


def test_missing_part_file_names_the_file(tmp_path):
    cfg = tmp_path / "p.yaml"
    cfg.write_text("provenance_policy: strict\nparts:\n  E1:\n    file: nope/missing.fa\n")
    with pytest.raises(SequenceError) as exc:
        parts_mod.load_parts(cfg)
    assert "missing.fa" in str(exc.value)


# ------------------------------------------------------------- report


@pytest.fixture(scope="module")
def built_report():
    parts, meta = parts_mod.load_parts(EXAMPLE_CFG)
    asm = construct.build(parts)
    m = score.measure(asm, temperatures=(37.0,))
    dec = score.score(m)
    payload = {"measurements": m, "decision": dec.to_dict()}
    return report.envelope(payload, provenance=meta, seed=123)


def test_envelope_metadata(built_report):
    r = built_report
    assert r["seed"] == 123
    assert r["circhbsag_version"]
    assert "ViennaRNA" in r["tool_licenses"]
    assert r["tool_licenses"]["circhbsag"]["license"] == "MIT"


def test_synthetic_input_is_labelled_as_software_test(built_report):
    assert "[DEMO — SEKUENS PLACEHOLDER]" in built_report["data_label"]
    assert "tidak boleh dikutip" in built_report["data_label"]


def test_user_supplied_input_is_labelled_in_silico_result():
    r = report.envelope({"x": 1}, provenance={"synthetic_parts": []}, seed=1)
    assert "[HASIL IN SILICO]" in r["data_label"]


def test_markdown_contains_decision_and_limitations(built_report):
    md = report.markdown(built_report)
    dec = built_report["results"]["decision"]
    assert dec["verdict"] in md
    assert "Gerbang keras" in md
    assert "Kriteria berbobot" in md
    assert "Catatan keterbatasan" in md
    assert "pseudoknot" in md
    # Setiap kriteria harus muncul bersama penalarannya.
    for c in dec["criteria"]:
        assert c["name"] in md
        assert c["rationale"][:40] in md


def test_html_is_self_contained(built_report):
    h = report.html(built_report)
    assert h.startswith("<!doctype html>")
    assert "prefers-color-scheme" in h
    assert "<script" not in h.lower()        # tidak butuh JS
    assert "http://" not in h and "https://" not in h  # tidak ada sumber eksternal
    assert built_report["results"]["decision"]["verdict"] in h


def test_write_all_creates_three_files(built_report, tmp_path):
    paths = report.write_all(built_report, tmp_path, stem="uji")
    assert set(paths) == {"json", "markdown", "html"}
    for p in paths.values():
        assert p.exists() and p.stat().st_size > 0
    back = json.loads(paths["json"].read_text())
    assert back["results"]["decision"]["verdict"] == built_report["results"]["decision"]["verdict"]


# ------------------------------------------------------------- CLI


def _run(args: list[str], cwd: Path) -> subprocess.CompletedProcess:
    env = {"PYTHONPATH": str(REPO / "src"), "PATH": "/usr/bin:/bin:/usr/local/bin",
           "HOME": str(cwd)}
    return subprocess.run(
        [sys.executable, "-m", "circhbsag.cli", *args],
        cwd=REPO, env=env, capture_output=True, text=True, timeout=1800,
    )


def test_cli_licenses_is_valid_json():
    proc = _run(["licenses"], REPO)
    assert proc.returncode == 0, proc.stderr
    data = json.loads(proc.stdout)
    assert data["circhbsag"]["license"] == "MIT"
    assert "license" in data["ViennaRNA"]


def test_cli_tahap0_writes_reports(tmp_path):
    out = tmp_path / "out"
    proc = _run(["tahap0", "--parts", str(EXAMPLE_CFG), "--out", str(out),
                 "--temperatures", "37", "--sensitivity", "0"], tmp_path)
    assert proc.returncode in (0, 2), proc.stderr   # 2 = NO-GO, tetap menulis laporan
    for name in ("tahap0.json", "tahap0.md", "tahap0.html"):
        assert (out / name).exists()
    data = json.loads((out / "tahap0.json").read_text())
    assert data["results"]["decision"]["verdict"] in {"GO", "GO-BERSYARAT", "NO-GO"}
    assert "[DEMO — SEKUENS PLACEHOLDER]" in data["data_label"]


def test_cli_codon_writes_fasta(tmp_path):
    out = tmp_path / "out"
    proc = _run(["codon", "--parts", str(EXAMPLE_CFG), "--out", str(out),
                 "--iterations", "300"], tmp_path)
    assert proc.returncode == 0, proc.stderr
    fa = (out / "orf_optimized.fa").read_text()
    assert fa.startswith(">")
    seq = "".join(fa.splitlines()[1:])
    assert seq.startswith("ATG") and len(seq) % 3 == 0


def test_cli_rejects_strict_config_with_missing_files(tmp_path):
    proc = _run(["tahap0", "--parts", str(REAL_CFG), "--out", str(tmp_path / "o")], tmp_path)
    assert proc.returncode != 0
    assert "reference" in (proc.stderr + proc.stdout) or "tidak ada" in (proc.stderr + proc.stdout)


def test_cli_codon_table_from_fasta(tmp_path):
    cds = tmp_path / "cds.fa"
    cds.write_text(">a\nATGGCCGGCTAA\n>b\nATGGCAGGATAA\n")
    out = tmp_path / "table.tsv"
    proc = _run(["codon-table", "--fasta", str(cds), "--out", str(out)], tmp_path)
    assert proc.returncode == 0, proc.stderr
    text = out.read_text()
    assert "[HASIL IN SILICO]" in text
    rows = [ln for ln in text.splitlines() if not ln.startswith("#") and "\t" in ln]
    assert len(rows) == 65  # header + 64 kodon
