"""Tes desain oligo dan filter statistik BSJ."""
from __future__ import annotations

import pytest

from circhbsag import bsj_stats, construct, primers
from circhbsag.seqio import revcomp

from .test_construct_score import mk_parts


@pytest.fixture(scope="module")
def asm():
    return construct.build(mk_parts())


# ------------------------------------------------------------- primers


def test_divergent_primers_are_found(asm):
    res = primers.divergent_primers(asm, n=3)
    assert res["pairs"], res["explain"]
    for p in res["pairs"]:
        assert p["left_seq"] in asm.circle or p["left_seq"] in asm.circle * 2
        assert 55.0 <= p["left_tm"] <= 65.0
        assert 55.0 <= p["right_tm"] <= 65.0
        assert 1 <= p["left_start_in_circle"] <= len(asm.circle)
        # left/right dari primer3 berorientasi konvergen (produk pendek juga ada pada linear);
        # pasangan divergen sejati = komplemen terbalik keduanya.
        assert p["convergent_product_nt"] == p["product_size"]
        assert p["divergent_fwd_5to3"] == primers._rc(p["right_seq"])
        assert p["divergent_rev_5to3"] == primers._rc(p["left_seq"])
        assert "konvergen" in p["expected_product_on_linear"]


def test_divergent_primers_sit_opposite_the_junction(asm):
    res = primers.divergent_primers(asm, n=3)
    n = len(asm.circle)
    anti = res["antipodal_position"]
    # Titik antipodal harus kira-kira setengah lingkaran dari BSJ.
    d = min(abs(anti - asm.bsj_position), n - abs(anti - asm.bsj_position))
    assert abs(d - n // 2) <= 1


def test_bsj_primers_span_junction(asm):
    res = primers.bsj_primers(asm, n=3, overlap=8)
    jsp = res["junction_spanning_primers"]
    assert jsp
    for p in jsp:
        assert p["nt_on_E1_side"] == p["nt_on_E2_side"]
        assert len(p["seq"]) == p["length"]
        # Separuh kiri berasal dari ujung E1, separuh kanan dari awal E2.
        half = p["length"] // 2
        assert asm.parts["E1"].seq.endswith(p["seq"][:half])
        assert asm.parts["E2"].seq.startswith(p["seq"][half:])


def test_bsj_primer_absent_from_precursor(asm):
    """Primer yang menumpang sambungan TIDAK boleh ada pada prekursor linear."""
    res = primers.bsj_primers(asm, n=1, overlap=10)
    for p in res["junction_spanning_primers"]:
        assert p["seq"] not in asm.precursor
        assert revcomp(p["seq"]) not in asm.precursor


def test_ivt_reverse_primer_is_revcomp_of_template_end(asm):
    res = primers.ivt_reverse_primer(asm.precursor, length=20, two_o_methyl_nt=2)
    assert res["primer_reverse_5to3"] == revcomp(asm.precursor)[:20]
    assert res["with_modifications"].startswith("m")
    assert res["with_modifications"].count("m") == 2
    assert res["length"] == 20


def test_rnase_h_probes_are_revcomp_of_target(asm):
    res = primers.rnase_h_probes(asm, n=3, probe_len=20)
    assert res["probes"]
    for p in res["probes"]:
        assert p["probe_dna_5to3"] == revcomp(p["rna_target"])
        assert asm.circle[p["start_1based"] - 1 : p["end_1based"]] == p["rna_target"]
        assert 55.0 <= p["tm"] <= 70.0
        assert 0.40 <= p["gc"] <= 0.60
        assert p["distance_to_bsj_nt"] >= 40
        assert p["expected_linear_product_nt"] == len(asm.circle)


def test_rnase_h_probes_are_separated(asm):
    res = primers.rnase_h_probes(asm, n=3, probe_len=20, min_sep=60)
    starts = sorted(p["start_1based"] for p in res["probes"])
    for a, b in zip(starts, starts[1:], strict=False):
        assert b - a >= 60


def test_design_all_keys(asm):
    res = primers.design_all(asm)
    assert set(res) == {"divergent", "bsj", "ivt_reverse", "rnase_h"}


# ------------------------------------------------------------- BSJ stats


def test_simulate_reads_counts_and_origins():
    circle = "ACGT" * 100
    linear = "TTTT" + circle + "AAAA"
    reads = bsj_stats.simulate_reads(circle, linear, n_circular=50, n_linear=50,
                                     read_len=100, error_rate=0.0, seed=1)
    assert len(reads) == 100
    assert sum(r.origin == "circular" for r in reads) == 50
    assert all(r.n_errors == 0 for r in reads)


def test_simulation_is_deterministic():
    circle, linear = "ACGTTGCA" * 50, "ACGTTGCA" * 60
    a = bsj_stats.simulate_reads(circle, linear, n_circular=20, n_linear=20, seed=9)
    b = bsj_stats.simulate_reads(circle, linear, n_circular=20, n_linear=20, seed=9)
    assert [r.seq for r in a] == [r.seq for r in b]


def test_error_rate_introduces_mismatches():
    circle, linear = "ACGTTGCA" * 50, "ACGTTGCA" * 60
    reads = bsj_stats.simulate_reads(circle, linear, n_circular=200, n_linear=0,
                                     read_len=150, error_rate=0.05, seed=3)
    assert sum(r.n_errors for r in reads) > 0


def test_read_posterior_increases_with_anchor_overlap():
    base = {
        "read": "r", "overlap_left_nt": 20, "overlap_right_nt": 20,
        "anchor_left_mismatch": 0, "anchor_right_mismatch": 0,
        "maps_to_linear_contiguously": False,
    }
    short = bsj_stats.read_posterior(base)
    longer = bsj_stats.read_posterior({**base, "overlap_left_nt": 40, "overlap_right_nt": 40})
    assert longer["log10_bayes_factor"] > short["log10_bayes_factor"]


def test_read_posterior_decreases_with_mismatches():
    base = {
        "read": "r", "overlap_left_nt": 30, "overlap_right_nt": 30,
        "anchor_left_mismatch": 0, "anchor_right_mismatch": 0,
        "maps_to_linear_contiguously": False,
    }
    clean = bsj_stats.read_posterior(base)
    noisy = bsj_stats.read_posterior({**base, "anchor_left_mismatch": 3})
    assert noisy["log10_bayes_factor"] < clean["log10_bayes_factor"]


def test_read_posterior_penalises_linear_explanation():
    """Bacaan yang dapat dijelaskan referensi linear harus kalah telak.

    Dibandingkan lewat log10 Bayes factor, bukan posterior: dengan jangkar
    30+30 nt kedua posterior sudah tersaturasi pada 1,0 sehingga tidak dapat
    membedakan keduanya.
    """
    base = {
        "read": "r", "overlap_left_nt": 30, "overlap_right_nt": 30,
        "anchor_left_mismatch": 0, "anchor_right_mismatch": 0,
        "maps_to_linear_contiguously": False,
    }
    uniq = bsj_stats.read_posterior(base)
    lin = bsj_stats.read_posterior({**base, "maps_to_linear_contiguously": True})
    assert lin["log10_bayes_factor"] < uniq["log10_bayes_factor"] - 5


def test_linear_explained_reads_are_excluded_by_the_caller():
    """Uji keunikan: kandidat yang cocok ke referensi linear tidak dihitung."""
    unit = "ACGTTGCAAC"
    circle = unit * 20                      # sangat berulang
    linear = "TTTT" + circle * 2 + "AAAA"   # sambungan palsu ada di templat linear
    reads = bsj_stats.simulate_reads(circle, linear, n_circular=60, n_linear=60,
                                     read_len=60, error_rate=0.0, seed=4)
    calls = bsj_stats.call_bsj(reads, circle, linear, anchor=20, max_mismatch=0)
    assert calls["n_excluded_linear_explained"] > 0
    assert all(not c["maps_to_linear_contiguously"] for c in calls["candidates"])


def test_posterior_respects_prior():
    base = {
        "read": "r", "overlap_left_nt": 25, "overlap_right_nt": 25,
        "anchor_left_mismatch": 0, "anchor_right_mismatch": 0,
        "maps_to_linear_contiguously": False,
    }
    low = bsj_stats.read_posterior(base, prior_circular=1e-6)
    high = bsj_stats.read_posterior(base, prior_circular=0.5)
    assert high["posterior_circular"] >= low["posterior_circular"]


def test_end_to_end_simulation_detects_circle(asm):
    res = bsj_stats.run_simulation(asm.circle, asm.precursor, n_circular=300,
                                   n_linear=300, read_len=150, anchor=20, seed=11)
    sup = res["support"]
    assert sup["n_candidates"] > 0
    conf = sup["confusion_on_simulated_origin"]
    # Bacaan dari prekursor linear tidak memuat sambungan balik, sehingga
    # tidak boleh lolos sebagai pendukung BSJ.
    assert conf["false_positive_from_linear"] == 0
    assert conf["true_positive_circular"] > 0
    assert sup["joint_posterior_circular"] > 0.5


def test_linear_only_reads_give_no_support(asm):
    res = bsj_stats.run_simulation(asm.circle, asm.precursor, n_circular=0,
                                   n_linear=400, read_len=150, anchor=20, seed=5)
    assert res["support"]["n_passing"] == 0
    assert res["support"]["joint_posterior_circular"] == 0.0


def test_external_tools_reported_not_assumed():
    av = bsj_stats.external_available()
    assert set(av) == set(bsj_stats.EXTERNAL_TOOLS)
    # Nilainya boleh None (tidak terpasang) — yang penting dilaporkan apa adanya.
    assert all(v is None or isinstance(v, str) for v in av.values())


def test_star_command_contains_chimeric_options():
    cmd = bsj_stats.write_star_cmd("reads.fq", "idx", "out/")
    assert "--chimSegmentMin" in cmd and "STAR" in cmd
