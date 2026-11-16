"""Tes pemindaian motif dan optimasi kodon."""
from __future__ import annotations

import pytest

from circhbsag import codon, motifs, seqio

# ---------------------------------------------------------------- motifs


def test_find_all_both_strands():
    seq = "AAAGGTCTCAAA"
    hits = motifs.find_all(seq, "GGTCTC")
    assert [h["start"] for h in hits if h["strand"] == "+"] == [3]
    # Komplemen balik GGTCTC = GAGACC, tidak ada di sekuens ini.
    assert not [h for h in hits if h["strand"] == "-"]
    rc_seq = "AAAGAGACCAAA"
    assert any(h["strand"] == "-" for h in motifs.find_all(rc_seq, "GGTCTC"))


def test_scan_restriction_detects_bsai_and_bsmbi():
    seq = "TTT" + motifs.RESTRICTION_SITES["BsaI"] + "TTT" + motifs.RESTRICTION_SITES["BsmBI_Esp3I"]
    rs = motifs.scan_restriction(seq)
    assert len(rs["BsaI"]) >= 1
    assert len(rs["BsmBI_Esp3I"]) >= 1
    clean = "ACAACAACAAACAACAAACAACAA"
    assert all(not v for v in motifs.scan_restriction(clean).values())


def test_t7_terminator_core_is_detected():
    seq = "AAAA" + motifs.T7_TERMINATOR_CORE + "AAAA"
    hits = motifs.scan_t7_terminator(seq)
    assert hits
    assert any(h["kind"] in ("Tphi_core", "intrinsic_like") for h in hits)


def test_t7_terminator_absent_in_clean_sequence():
    assert motifs.scan_t7_terminator("ACAACAACAAACAACAAACAACAAACAACAA") == []


def test_rrach_sites_and_position():
    # GGACT adalah RRACH; adenosin sasaran ada di indeks ke-2 motif.
    seq = "TTGGACTTT"
    sites = motifs.scan_m6a_rrach(seq)
    assert len(sites) == 1
    s = sites[0]
    assert s["motif"] == "GGACT"
    assert seq[s["a_position"]] == "A"


def test_rrach_rejects_non_consensus():
    assert motifs.scan_m6a_rrach("TTTTCACTTTT") == []


def test_m6a_density_scales_with_length():
    one = "TTGGACTTT"
    assert motifs.m6a_density(one) == pytest.approx(100.0 / len(one))


def test_gc_window_violations():
    assert motifs.gc_window_violations("GC" * 60, window=50, step=25) != []
    assert motifs.gc_window_violations("ACAACAAACAACAAACAACAAACAACAAACAACAAACAACAAACAACAAACAACAA" * 2,
                                       window=50, step=25) != []
    balanced = ("ACGT" * 30)
    assert motifs.gc_window_violations(balanced, window=50, step=25) == []


def test_upstream_atg_in_frame_window():
    # ATG sekerangka 30 nt di hulu ORF.
    seq = "ATG" + "CCC" * 9 + "ATGGCCTGA"
    orf_start = len(seq) - 9
    assert motifs.upstream_atg_in_frame(seq, orf_start, within=60) == [0]
    assert motifs.upstream_atg_in_frame(seq, orf_start, within=10) == []


def test_upstream_atg_out_of_frame():
    seq = "CC" + "ATG" + "CCCCC" + "ATGGCCTGA"
    orf_start = len(seq) - 9
    hits = motifs.upstream_atg_out_of_frame(seq, orf_start, within=30)
    assert hits == [2]


def test_summarize_keys():
    s = "ACAACAAACAACAAACAACAAACAACAAACAACAA"
    out = motifs.summarize(s)
    for key in ("length", "restriction", "t7_terminator", "m6a_sites",
                "m6a_density_per_100nt", "avoid_motifs", "gc_violations"):
        assert key in out


# ---------------------------------------------------------------- codon


@pytest.fixture(scope="module")
def usage():
    return codon.load_codon_usage()


def test_codon_table_complete(usage):
    assert len(usage) == 64
    assert all(c in usage for c in seqio.CODON_TABLE)


def test_relative_adaptiveness_bounds(usage):
    w = codon.relative_adaptiveness(usage)
    assert all(0.0 <= v <= 1.0 for v in w.values())
    # Tiap keluarga sinonim harus punya tepat satu kodon dengan w = 1.
    for codons in seqio.SYN.values():
        assert max(w[c] for c in codons) == pytest.approx(1.0)


def test_cai_of_optimal_sequence_is_one(usage):
    w = codon.relative_adaptiveness(usage)
    protein = "MAKELVGWY"
    best = codon.naive_optimize(protein, w)
    assert codon.cai(best, w) == pytest.approx(1.0)


def test_cai_of_worst_sequence_is_low(usage):
    w = codon.relative_adaptiveness(usage)
    protein = "LLLRRRSSS"
    worst = "".join(min(seqio.SYN[aa], key=lambda c: w[c]) for aa in protein)
    assert codon.cai(worst, w) < 0.5


def test_optimize_preserves_protein():
    protein = "MAKELVGWYCDEFHIKLMNPQRSTV"
    res = codon.optimize(protein, iterations=800, seed=1)
    assert seqio.translate(res["orf_dna"]).rstrip("*") == protein
    assert res["orf_dna"].startswith("ATG")
    assert res["orf_dna"].endswith("TGA")
    assert res["length_nt"] == 3 * (len(protein) + 1)


def test_optimize_is_deterministic_for_fixed_seed():
    protein = "MAKELVGWYCDEFHIKLMNPQ"
    a = codon.optimize(protein, iterations=500, seed=42)["orf_dna"]
    b = codon.optimize(protein, iterations=500, seed=42)["orf_dna"]
    assert a == b


def test_optimize_removes_restriction_sites():
    # Protein yang kodon paling disukainya mudah membentuk GGTCTC/CGTCTC
    protein = "MGLSAGLSAGLSAGLSAGLSA"
    res = codon.optimize(protein, iterations=4000, seed=7)
    assert res["cost"]["metrics"]["restriction_sites"] == 0


def test_optimize_reduces_m6a_versus_naive():
    protein = "MKKNKNKDDEEGGAAKKNKNKDDEEGGAAKKNKNK"
    res = codon.optimize(protein, iterations=6000, seed=3)
    assert res["cost"]["metrics"]["m6a_sites"] <= res["baseline_naive"]["m6a_sites"]


def test_optimize_rejects_bad_protein():
    with pytest.raises(seqio.SequenceError):
        codon.optimize("MAKEXLV", iterations=10)
    with pytest.raises(seqio.SequenceError):
        codon.optimize("", iterations=10)


def test_codon_stats_keys():
    res = codon.optimize("MAKELVGWY", iterations=300, seed=5)
    st = codon.codon_stats(res["orf_dna"])
    for key in ("length_nt", "cai", "gc", "gc3", "m6a_sites"):
        assert key in st
    assert 0.0 <= st["cai"] <= 1.0
    assert 0.0 <= st["gc"] <= 1.0
