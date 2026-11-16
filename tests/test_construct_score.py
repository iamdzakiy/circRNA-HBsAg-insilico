"""Tes perakitan konstruk, skor komposit, dan keputusan GO/NO-GO."""
from __future__ import annotations

import pytest

from circhbsag import construct, score
from circhbsag.seqio import SequenceError

P = construct.Part


def mk_parts(**over) -> dict:
    """Bagian minimal yang sah untuk perakitan (sekuens uji, bukan biologis)."""
    base = {
        "T7_promoter": P("T7_promoter", "TAATACGACTCACTATA"),
        "gg_leader": P("gg_leader", "GGG"),
        "HA5": P("HA5", "ACAACAAACAACAAACAACAAACAACAAAC"),
        "intron_3half": P("intron_3half", "CCGGCAACTAACAATCAACTACAATCAACCGG" * 3,
                          features={"IGS": (90, 96)}),
        "E2": P("E2", "ACAACTACAATCAACTACAATCAACTACAATC"),
        "spacer5": P("spacer5", "ACTAACTACAATCAACTACAATC"),
        "IRES": P("IRES", "GGGCCGCAACTAAACCAATCAACTACAATCAACCGGCC" * 6),
        "ORF": P("ORF", "ATG" + "GCCACC" * 20 + "TGA"),
        "spacer3": P("spacer3", "AAACAAAGAAACAAAGAAACAAAG"),
        "E1": P("E1", "CAACTACAATCAACTACAATCAACTACAATCA"),
        "intron_5half": P("intron_5half", "GGTGGT" + "CAACTACAATCAACTACAATCAACCGGCC" * 3,
                          features={"IGS": (0, 6)}),
        "HA3": P("HA3", "GTTTGTTGTTTGTTGTTTGTTGTTTGTTGT"),
    }
    base.update(over)
    return base


# ------------------------------------------------------------- construct


def test_build_orders_parts_and_tracks_coords():
    asm = construct.build(mk_parts())
    # Lingkaran = E2 + spacer5 + IRES + ORF + spacer3 + E1
    expected = sum(
        len(asm.parts[k].seq) for k in ("E2", "spacer5", "IRES", "ORF", "spacer3", "E1")
    )
    assert len(asm.circle) == expected
    assert asm.circle.startswith(asm.parts["E2"].seq)
    assert asm.circle.endswith(asm.parts["E1"].seq)
    assert asm.bsj_position == len(asm.circle)
    # Prekursor memuat lingkaran dalam urutan tersisip di antara intron.
    assert asm.parts["intron_3half"].seq in asm.precursor
    assert asm.precursor.index(asm.parts["intron_3half"].seq) < asm.precursor.index(
        asm.parts["E2"].seq
    )


def test_coordinates_match_subsequences():
    asm = construct.build(mk_parts())
    for name in ("E1", "E2", "IRES", "ORF"):
        assert asm.subseq("circle", name) == asm.parts[name].seq
        assert asm.subseq("precursor", name) == asm.parts[name].seq


def test_feature_coordinates_are_offset():
    asm = construct.build(mk_parts())
    a, b = asm.precursor_coords["intron_5half:IGS"]
    assert asm.precursor[a - 1 : b] == "GGTGGT"


def test_build_requires_core_parts():
    parts = mk_parts()
    del parts["IRES"]
    with pytest.raises(SequenceError) as exc:
        construct.build(parts)
    assert "IRES" in str(exc.value)


def test_build_non_strict_allows_missing():
    parts = mk_parts()
    del parts["IRES"]
    asm = construct.build(parts, strict=False)
    assert "IRES" not in asm.circle_coords


def test_permuted_template_checks_t7_start():
    asm = construct.build(mk_parts())
    tpl = construct.permuted_dna_template(asm)
    assert tpl["t7_promoter_found"] is True
    assert tpl["first_three_nt"] == "GGG"
    assert tpl["notes"] == []
    # Tanpa gg_leader, pipeline harus memperingatkan.
    parts = mk_parts(gg_leader=P("gg_leader", ""))
    tpl2 = construct.permuted_dna_template(construct.build(parts))
    assert tpl2["first_three_nt"] != "GGG"
    assert tpl2["notes"]


def test_bsj_amplicon_centres_the_junction():
    asm = construct.build(mk_parts())
    amp = construct.bsj_amplicon(asm, flank=60)
    assert len(amp["amplicon"]) == 120
    assert amp["junction_offset_in_amplicon"] == 60
    # Sisi kiri sambungan berasal dari ujung E1, sisi kanan dari awal E2.
    assert amp["left_of_junction"].endswith(asm.parts["E1"].seq[-10:])
    assert amp["right_of_junction"].startswith(asm.parts["E2"].seq[:10])


def test_circle_rotations_preserve_length_and_content():
    asm = construct.build(mk_parts())
    rots = construct.circle_rotations(asm, step=97)
    assert all(len(r) == len(asm.circle) for r in rots)
    assert all(sorted(r) == sorted(asm.circle) for r in rots)


def test_guide_sequence_pairs_reports_igs():
    g = construct.guide_sequence_pairs(construct.build(mk_parts()))
    assert g["P1"]["igs_seq"] == "GGTGGT"
    assert "P10" in g


# ------------------------------------------------------------- score


@pytest.fixture(scope="module")
def measured():
    asm = construct.build(mk_parts())
    return asm, score.measure(asm, temperatures=(37.0, 55.0))


def test_measure_structure(measured):
    _asm, m = measured
    for key in ("assembly", "folds", "junction", "accessibility", "ires_integrity",
                "dsrna", "start_codon", "motifs", "orf_checks", "gc", "homopolymer"):
        assert key in m
    assert set(m["folds"]) == {"circle_37C", "circle_55C", "precursor_37C", "precursor_55C"}
    assert m["folds"]["circle_37C"]["circular"] is True
    assert m["folds"]["precursor_37C"]["circular"] is False


def test_measure_accessibility_includes_splice_sites(measured):
    _asm, m = measured
    regions = m["accessibility"]["37C"]
    assert "splice5_site" in regions and "splice3_site" in regions
    assert "intron_5half_IGS" in regions and "intron_3half_IGS" in regions


def test_ires_integrity_between_zero_and_one(measured):
    _asm, m = measured
    for v in m["ires_integrity"].values():
        assert 0.0 <= v["identity"] <= 1.0


def test_score_returns_normalised_subscores(measured):
    _asm, m = measured
    dec = score.score(m)
    assert dec.verdict in {"GO", "GO-BERSYARAT", "NO-GO"}
    assert 0.0 <= dec.composite <= 1.0
    assert len(dec.criteria) == len(score.DEFAULT_WEIGHTS)
    for c in dec.criteria:
        assert 0.0 <= c.subscore <= 1.0
        assert c.rationale  # setiap kriteria wajib menjelaskan dirinya
    assert set(dec.hard_gates) == set(score.HARD_GATES)


def test_composite_is_weighted_mean(measured):
    _asm, m = measured
    dec = score.score(m)
    total_w = sum(c.weight for c in dec.criteria)
    manual = sum(c.weight * c.subscore for c in dec.criteria) / total_w
    assert dec.composite == pytest.approx(manual)


def test_weights_sum_to_one():
    assert sum(score.DEFAULT_WEIGHTS.values()) == pytest.approx(1.0)


def test_hard_gate_failure_forces_no_go(measured):
    """Situs BsaI internal harus memaksa NO-GO walau kriteria lain bagus."""
    parts = mk_parts(spacer5=P("spacer5", "ACTAAGGTCTCACTACAATCAACTACAATC"))
    asm = construct.build(parts)
    m = score.measure(asm, temperatures=(37.0,))
    dec = score.score(m)
    assert "no_internal_typeIIS" in [k for k, v in dec.hard_gates.items() if not v["passed"]]
    assert dec.verdict == "NO-GO"


def test_premature_t7_terminator_forces_no_go():
    from circhbsag.motifs import T7_TERMINATOR_CORE

    parts = mk_parts(spacer3=P("spacer3", T7_TERMINATOR_CORE))
    asm = construct.build(parts)
    dec = score.score(score.measure(asm, temperatures=(37.0,)))
    assert dec.verdict == "NO-GO"
    assert not dec.hard_gates["no_premature_T7_terminator"]["passed"]


def test_oversized_circle_forces_no_go(measured):
    _asm, m = measured
    dec = score.score(m, thresholds={"max_circle_nt": 10})
    assert not dec.hard_gates["circle_length_within_limit"]["passed"]
    assert dec.verdict == "NO-GO"


def test_in_frame_upstream_atg_near_start_forces_no_go():
    # ATG sekerangka 9 nt di hulu ORF, pada spacer5 (dalam jendela 60 nt).
    parts = mk_parts(spacer5=P("spacer5", "ACTAACTACAATCAACTACAATC"))
    asm = construct.build(parts)
    # Suntikkan ATG tepat di hulu kodon awal, di ujung 3' wilayah IRES.
    ires = asm.parts["IRES"].seq
    ires = ires[:-9] + "ATG" + ires[-6:]
    parts2 = mk_parts(IRES=P("IRES", ires))
    asm2 = construct.build(parts2)
    start2 = asm2.circle_coords["ORF"][0] - 1
    if (start2 - (start2 - 9)) % 3 == 0:  # pastikan sekerangka
        dec = score.score(score.measure(asm2, temperatures=(37.0,)))
        gate = dec.hard_gates["orf_intact_and_unique_start"]
        assert gate["passed"] is False
        assert dec.verdict == "NO-GO"


def test_thresholds_shift_verdict(measured):
    _asm, m = measured
    strict = score.score(m, thresholds={"go": 1.01, "conditional": 1.0})
    assert strict.verdict == "NO-GO"
    lenient = score.score(m, thresholds={"go": 0.0, "conditional": 0.0})
    if not lenient.to_dict()["hard_gates_failed"]:
        assert lenient.verdict == "GO"


def test_sensitivity_reports_distribution(measured):
    _asm, m = measured
    s = score.sensitivity(m, n=25, seed=1)
    assert s["n"] == 25
    assert s["composite_min"] <= s["composite_median"] <= s["composite_max"]
    assert sum(s["verdict_counts"].values()) == 25


def test_ramp_helpers():
    assert score._ramp(0.0, 0.0, 1.0) == 0.0
    assert score._ramp(1.0, 0.0, 1.0) == 1.0
    assert score._ramp(0.5, 0.0, 1.0) == pytest.approx(0.5)
    assert score._ramp_down(0.0, 0.0, 1.0) == 1.0
    assert score._band(0.5, 0.4, 0.6, 0.1) == 1.0
    assert score._band(0.3, 0.4, 0.6, 0.1) == pytest.approx(0.0)
    assert score._band(0.35, 0.4, 0.6, 0.1) == pytest.approx(0.5)
