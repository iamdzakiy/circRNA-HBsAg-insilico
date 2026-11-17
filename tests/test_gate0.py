"""Tests for the nine Gate 0 tests.

These tests check the two properties that make the report trustworthy: a test
whose tool is missing must report that it is missing rather than return a number,
and ``UNAVAILABLE`` must never be counted as a pass.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from circhbsag import construct, gate0, score
from circhbsag import parts as parts_mod

REPO = Path(__file__).resolve().parents[1]
EXAMPLE_CFG = REPO / "config" / "parts.example.yaml"


@pytest.fixture(scope="module")
def gate():
    parts, meta = parts_mod.load_parts(EXAMPLE_CFG)
    asm = construct.build(parts)
    m = score.measure(asm, temperatures=(37.0,))
    dec = score.score(m)
    g = gate0.run(asm, m, arm_lengths=(20, 30, 40))
    return g, asm, m, dec, meta


def test_nine_tests_are_present_and_ordered(gate):
    g, *_ = gate
    ids = [t["id"] for t in g["tests"]]
    assert len(ids) == 9
    assert ids[0].startswith("uji1") and ids[-1].startswith("uji9")
    assert len(set(ids)) == 9


def test_every_test_has_a_threshold_and_its_justification(gate):
    g, *_ = gate
    for t in g["tests"]:
        assert t["threshold"], t["id"]
        assert t["threshold_label"] in {"LITERATUR", "PERKIRAAN"}, t["id"]
        assert len(t["threshold_rationale"]) > 60, t["id"]
        # A literature threshold must name the paper it came from.
        if t["threshold_label"] == "LITERATUR":
            assert t["refs"], t["id"]


def test_verdicts_are_from_the_closed_set(gate):
    g, *_ = gate
    allowed = {gate0.PASS, gate0.REVISE, gate0.FAIL, gate0.UNAVAILABLE}
    assert {t["verdict"] for t in g["tests"]} <= allowed


def test_missing_infernal_is_reported_not_guessed(gate):
    """Uji 2 must say the tool is missing and hand back the command to run."""
    g, *_ = gate
    t = next(x for x in g["tests"] if x["id"] == "uji2_ribozyme_integrity")
    assert t["verdict"] == gate0.UNAVAILABLE
    assert t["metrics"]["status"] == gate0.UNAVAILABLE
    assert "cmsearch" in t["metrics"]["command"]
    assert "RF00028" in t["metrics"]["command"]
    assert "bit_score" not in str(t["metrics"]), "tidak boleh ada angka tebakan"


def test_unavailable_is_not_counted_as_pass(gate):
    g, *_ = gate
    assert g["counts"][gate0.UNAVAILABLE] >= 1
    assert g["decision"] != "PASS"
    assert g["decision"] == "BELUM-LENGKAP"


def test_cmsearch_without_model_never_invents_hits(tmp_path):
    out = gate0.cmsearch("ACGT" * 50, str(tmp_path / "nope.cm"))
    assert out["status"] == gate0.UNAVAILABLE
    assert out["covariance_model_found"] is False
    assert "hits" not in out


def test_arm_duplex_is_more_negative_for_longer_arms():
    """A longer perfect complement must not bind less tightly."""
    arm5 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    arm3 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    scan = gate0.arm_length_scan(arm5, arm3, lengths=(20, 30, 40), temperature_c=55.0)
    dg = [p["duplex_kcal"] for p in scan if p["available"]]
    assert len(dg) == 3
    assert dg[0] >= dg[1] >= dg[2]


def test_cofold_separates_binding_from_monomer_structure():
    a = "GGGGGGGGGGCCCCCCCCCC"
    out = gate0.cofold_dg(a, a, temperature_c=37.0)
    # The hairpin each monomer forms on its own must be subtracted out.
    assert out["binding_kcal"] == pytest.approx(
        out["dimer_kcal"] - out["monomer_a_kcal"] - out["monomer_b_kcal"], abs=1e-6
    )
    assert out["monomer_a_kcal"] <= 0.0


def test_polyu_runs_finds_only_runs_at_or_above_the_limit():
    assert polyu_lengths("ATTTTTA", 5) == [5]
    assert polyu_lengths("ATTTTA", 5) == []
    assert polyu_lengths("TTTTTTT", 5) == [7]


def polyu_lengths(seq: str, n: int) -> list[int]:
    return [r["length"] for r in gate0.polyu_runs(seq, min_run=n)]


def test_cryptic_splice_scan_requires_a_plausible_intron_length():
    # A GT and an AG two nucleotides apart is not a candidate intron.
    assert gate0.cryptic_splice_sites("GTAG", min_intron=70) == []
    seq = "GT" + "A" * 100 + "AG"
    hits = gate0.cryptic_splice_sites(seq, min_intron=70, max_intron=2000)
    assert hits and hits[0]["donor_pos"] == 1


def test_design_lock_carries_the_numbers_other_deliverables_quote(gate):
    g, asm, m, dec, meta = gate
    lock = gate0.design_lock(g, m, dec.to_dict(), provenance=meta)
    assert lock["data_status"] == "DEMO-PLACEHOLDER"
    assert "tidak boleh dikutip" in lock["data_label"]
    assert lock["projection_label"] == "Proyeksi berbasis literatur — bukan data eksperimen"
    assert lock["lengths_nt"]["precursor"] == len(asm.precursor)
    assert lock["lengths_nt"]["circle"] == len(asm.circle)
    assert lock["gate0"]["decision"] == g["decision"]
    assert lock["composite_verdict"] == dec.verdict
    assert lock["ribozyme_check"] == gate0.UNAVAILABLE


def test_markdown_report_states_the_decision_and_every_threshold(gate):
    g, *_ = gate
    md = gate0.markdown(g)
    assert g["decision"] in md
    for t in g["tests"]:
        assert t["name"] in md
        assert t["threshold"][:30] in md
    assert "BELUM BISA DIJALANKAN" in md
