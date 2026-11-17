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


#: A path that cannot exist, so the "tool or model missing" branch is exercised
#: deterministically whether or not Infernal happens to be installed here.
NO_CM = "/nonexistent/RF00028.cm"


@pytest.fixture(scope="module")
def gate():
    """Gate 0 run with no covariance model available."""
    parts, meta = parts_mod.load_parts(EXAMPLE_CFG)
    asm = construct.build(parts)
    m = score.measure(asm, temperatures=(37.0,))
    dec = score.score(m)
    g = gate0.run(asm, m, cm_path=NO_CM, arm_lengths=(20, 30, 40))
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


def test_missing_model_is_reported_not_guessed(gate):
    """Uji 2 must say what is missing and hand back the command to run."""
    g, *_ = gate
    t = next(x for x in g["tests"] if x["id"] == "uji2_ribozyme_integrity")
    assert t["verdict"] == gate0.UNAVAILABLE
    assert t["metrics"]["status"] == gate0.UNAVAILABLE
    assert "cmsearch" in t["metrics"]["command"]
    assert "RF00028" in t["metrics"]["command"]
    assert "bit_score" not in str(t["metrics"]), "tidak boleh ada angka tebakan"


def test_unavailable_is_not_counted_as_pass(gate):
    """A test that could not run must never be rounded up to a pass."""
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


# ------------------------------------------------- Uji 2 with Infernal present

CM = REPO / "data" / "reference" / "RF00028.cm"
POSITIVE = REPO / "tests" / "data" / "RF00028_positive_control.fa"
HAVE_INFERNAL = __import__("shutil").which("cmsearch") is not None


@pytest.mark.skipif(not (HAVE_INFERNAL and CM.exists()), reason="butuh Infernal dan RF00028.cm")
def test_uji2_detects_a_real_group_I_intron():
    """The positive control: a gate that always fails would be worthless.

    This is the test that gives Uji 2 its meaning. Without it, a FAIL on the
    synthetic construct could just as easily mean the covariance model is
    broken, the sequence was passed in the wrong alphabet, or the threshold is
    unreachable. With it, a FAIL means what it claims to mean.
    """
    seq = "".join(
        line.strip() for line in POSITIVE.read_text().splitlines() if not line.startswith(">")
    )
    out = gate0.cmsearch(seq, str(CM), e_value=10)
    assert out["status"] == "OK"
    assert out["n_hits"] >= 1, "model gagal mengenali intron grup I sungguhan"
    best = out["hits"][0]
    assert best["e_value"] < 1e-5
    assert best["bit_score"] > 50


@pytest.mark.skipif(not (HAVE_INFERNAL and CM.exists()), reason="butuh Infernal dan RF00028.cm")
def test_uji2_rejects_a_sequence_that_only_bears_the_name():
    """The synthetic 'intron' in the demo data must not be recognised."""
    parts, _ = parts_mod.load_parts(EXAMPLE_CFG)
    fake_intron = parts["intron_5half"].seq + parts["intron_3half"].seq
    out = gate0.cmsearch(fake_intron, str(CM), e_value=10)
    assert out["status"] == "OK"
    assert out["n_hits"] == 0


@pytest.mark.skipif(not (HAVE_INFERNAL and CM.exists()), reason="butuh Infernal dan RF00028.cm")
def test_uji2_fail_is_not_patchable():
    """A FAIL on Uji 2 must route to the backup circularisation path."""
    parts, meta = parts_mod.load_parts(EXAMPLE_CFG)
    asm = construct.build(parts)
    m = score.measure(asm, temperatures=(37.0,))
    g = gate0.run(asm, m, cm_path=str(CM), arm_lengths=(20, 40))
    t = next(x for x in g["tests"] if x["id"] == "uji2_ribozyme_integrity")
    assert t["verdict"] == gate0.FAIL
    assert g["decision"] == "FAIL-PINDAH-JALUR"
    assert "cadangan" in g["statement"]
