"""Smoke tests for scripts/: honest failure modes matter more than happy paths."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def run(script: str, *args: str, cwd: Path = REPO):
    return subprocess.run([sys.executable, str(REPO / "scripts" / script), *args],
                          capture_output=True, text=True, cwd=cwd, timeout=300)


def test_every_script_has_help():
    for p in sorted((REPO / "scripts").glob("[0-9][0-9]_*.py")):
        r = run(p.name, "--help")
        assert r.returncode == 0, p.name
    assert len(list((REPO / "scripts").glob("[0-9][0-9]_*.py"))) == 13


def test_epitope_script_refuses_to_invent_predictions(tmp_path):
    r = run("09_epitope_hla_coverage.py", "--results", str(tmp_path))
    assert r.returncode == 1
    out = json.loads((tmp_path / "tables" / "epitope_hla_coverage.json").read_text())
    assert out["status"] == "UNAVAILABLE"
    assert "per_allele_binder_fraction" not in out
    assert (tmp_path / "tables" / "epitope_hla_coverage.json.params.json").exists()


def test_epitope_script_computes_from_supplied_table(tmp_path):
    tsv = tmp_path / "pred.tsv"
    tsv.write_text("peptide\tallele\trank_percent\nAAA\tHLA-A*02:01\t0.5\nCCC\tHLA-A*02:01\t5.0\n")
    r = run("09_epitope_hla_coverage.py", "--results", str(tmp_path), "--predictions", str(tsv))
    assert r.returncode == 0
    out = json.loads((tmp_path / "tables" / "epitope_hla_coverage.json").read_text())
    assert out["per_allele_binder_fraction"]["HLA-A*02:01"] == 0.5


def test_fetch_without_accessions_marks_rows_not_substitutes(tmp_path):
    acc = tmp_path / "acc.tsv"
    acc.write_text("part\taccession\tdb\nX\t\tnuccore\n")
    r = run("01_fetch_sequences.py", "--accessions", str(acc), "--out", str(tmp_path / "o"))
    assert r.returncode == 1
    prov = (tmp_path / "o" / "provenance.tsv").read_text()
    assert "[PERLU VERIFIKASI]" in prov and "NO-ACCESSION" in prov
    assert not list((tmp_path / "o").glob("X.fa"))


def test_downstream_script_explains_missing_gate0(tmp_path):
    r = run("06_ires_integrity.py", "--results", str(tmp_path))
    assert r.returncode != 0
    assert "04_gate0_circularization_check.py" in (r.stderr + r.stdout)


def test_workbook_export_flattens_existing_results_and_skips_missing(tmp_path):
    import shutil

    import pytest

    src = REPO / "results" / "demo"
    if not (src / "gate0.json").exists():
        pytest.skip("butuh results/demo")
    for f in ("gate0.json", "gate0_iterations.json", "optimization.json", "codon.json", "primers.json"):
        if (src / f).exists():
            shutil.copy(src / f, tmp_path / f)
    r = run("13_export_workbook.py", "--results", str(tmp_path))
    assert r.returncode == 0, r.stderr
    import pandas as pd

    sheets = pd.read_excel(tmp_path / "tables" / "circhbsag_hasil_lengkap.xlsx", sheet_name=None)
    assert {"Ringkasan", "Gate0_uji", "Varian_Optuna", "Primer_pasangan"} <= set(sheets)
    assert len(sheets["Gate0_uji"]) == 9
    # missing sources are reported, never invented
    assert "bsj_simulation.json" in " ".join(sheets["Ringkasan"]["Nilai"].astype(str))
    assert (tmp_path / "tables" / "csv" / "Varian_Optuna.csv").exists()
