"""Tes lipatan RNA dan aksesibilitas.

Tes di sini memeriksa SIFAT yang harus benar (invarian), bukan angka hasil
tertentu, agar tetap sahih bila versi ViennaRNA berubah. Satu-satunya angka
yang dikunci adalah yang dihitung ulang di dalam tes itu sendiri.
"""
from __future__ import annotations

import pytest

from circhbsag import access, fold

# Jepit rambut sederhana dengan batang 10 bp dan lengkung 6 nt.
STEM = "GGGGCCCCGG"
HAIRPIN = STEM + "AAAAAA" + "CCGGGGCCCC"
RANDOMISH = "ACUAACUACAAUCAACUACAAUCAACUAACAAUCAACUACAACUAACUACAAUCAACUACAAU"


def test_fold_returns_consistent_shapes():
    r = fold.fold(HAIRPIN, circular=False)
    assert r.sequence_length == len(HAIRPIN)
    assert len(r.mfe_structure) == len(HAIRPIN)
    assert set(r.mfe_structure) <= set("().")
    assert r.mfe_kcal < 0
    assert r.mfe_per_nt == pytest.approx(r.mfe_kcal / len(HAIRPIN))
    # Energi bebas ensemble tidak boleh lebih tinggi daripada MFE.
    assert r.ensemble_free_energy_kcal <= r.mfe_kcal + 1e-6
    assert 0.0 <= r.frequency_mfe_in_ensemble <= 1.0


def test_circular_mode_changes_the_model():
    """Mode sirkular harus benar-benar mengubah model energi, bukan hanya label.

    Catatan: MFE sirkular TIDAK selalu lebih rendah daripada MFE linear. Untuk
    sekuens tanpa struktur, ViennaRNA membebankan energi lengkung eksterior pada
    lingkaran sehingga MFE sirkular dapat bernilai positif, sementara MFE linear
    mendekati nol. Karena itu tes ini memeriksa bahwa hasil kedua mode berbeda
    dan keduanya sah, bukan mengurutkan keduanya.
    """
    lin = fold.fold(RANDOMISH, circular=False)
    circ = fold.fold(RANDOMISH, circular=True)
    assert lin.circular is False and circ.circular is True
    assert (circ.mfe_kcal, circ.mfe_structure) != (lin.mfe_kcal, lin.mfe_structure)
    for r in (lin, circ):
        assert len(r.mfe_structure) == len(RANDOMISH)
        assert r.ensemble_free_energy_kcal <= r.mfe_kcal + 1e-6


def test_pair_probabilities_span_the_whole_molecule():
    """Pasangan jarak-jauh terdeteksi dan terindeks dengan benar.

    Sekuens disusun agar 12 nt pertama komplementer dengan 12 nt terakhir,
    dipisahkan wilayah tak-berstruktur. Tes ini memeriksa penanganan indeks
    pasangan jarak-jauh (i kecil, j besar) pada mode sirkular; ia TIDAK
    mengklaim pasangan tersebut khas lingkaran, karena lipatan linear pun dapat
    memasangkan kedua ujung ini.
    """
    left = "GGGCGCGGCCGC"
    middle = "ACAACAAACAACAAACAACAAACAACAAACAACAA"
    right = "GCGGCCGCGCCC"  # komplemen balik dari `left`
    seq = left + middle + right
    circ = fold.fold(seq, circular=True)
    n = len(seq)
    # Pasangan yang menghubungkan basa dekat ujung 3' dengan basa dekat ujung 5'
    across = [
        (i, j) for (i, j), p in circ.pair_probabilities.items()
        if p > 0.5 and i <= 12 and j >= n - 12
    ]
    assert across, f"tidak ada pasangan melintasi batas pada struktur {circ.mfe_structure}"


def test_higher_temperature_destabilizes():
    cold = fold.fold(HAIRPIN, circular=False, temperature_c=37.0)
    hot = fold.fold(HAIRPIN, circular=False, temperature_c=75.0)
    assert hot.mfe_kcal > cold.mfe_kcal


def test_pair_probabilities_are_valid():
    r = fold.fold(HAIRPIN, circular=False)
    for (i, j), p in r.pair_probabilities.items():
        assert 1 <= i < j <= r.sequence_length
        assert 0.0 < p <= 1.0 + 1e-9
    # Probabilitas total berpasangan per posisi tidak boleh > 1.
    unpaired = fold.unpaired_probability(r, list(range(1, r.sequence_length + 1)))
    assert all(-1e-6 <= v <= 1.0 + 1e-6 for v in unpaired.values())


def test_junction_pairing_window_wraps_circularly():
    r = fold.fold(RANDOMISH, circular=True)
    j = fold.junction_pairing(r, junction_pos=2, span=5)
    assert len(j["positions"]) == 10
    # Jendela harus melingkar melewati batas sekuens.
    assert max(j["positions"]) > r.sequence_length - 5
    assert 0.0 <= j["mean_unpaired_probability"] <= 1.0
    assert j["mean_unpaired_probability"] + j["mean_paired_probability"] == pytest.approx(1.0, abs=1e-4)


def test_fold_profile_two_temperatures():
    prof = fold.fold_profile(RANDOMISH, (37.0, 55.0))
    assert set(prof) == {"37C", "55C"}
    assert prof["55C"].mfe_kcal >= prof["37C"].mfe_kcal


def test_subsequence_identity_perfect_and_zero():
    s = "((((....))))"
    assert fold.subsequence_structure_identity(s, s, 0) == 1.0
    assert fold.subsequence_structure_identity("............", "((((((((((((", 0) == 0.0


def test_helices_from_structure():
    hel = access.helices_from_structure("((((....))))")
    assert hel[0]["length"] == 4
    assert hel[0]["i5"] == 1 and hel[0]["j3"] == 12
    # Bulge memotong heliks menjadi dua.
    hel2 = access.helices_from_structure("(((.(((....))))))")
    assert max(h["length"] for h in hel2) == 3


def test_longest_dsrna_and_violations():
    long_stem = "(" * 35 + "...." + ")" * 35
    d = access.longest_dsrna(long_stem)
    assert d["longest_bp"] == 35
    assert len(access.dsrna_violations(long_stem, 30)) == 1
    assert access.dsrna_violations(long_stem, 40) == []
    assert access.longest_dsrna("....")["longest_bp"] == 0


def test_plfold_accessibility_ranges():
    seq = RANDOMISH * 3
    acc = access.accessibility(seq, {"awal": (5, 25), "akhir": (len(seq) - 25, len(seq) - 5)},
                               ulength=8, window=60)
    for v in acc.values():
        assert v["n"] > 0
        assert 0.0 <= v["min"] <= v["mean"] <= v["max"] <= 1.0


def test_plfold_hairpin_less_accessible_than_unstructured():
    structured = (STEM + "AAAAAA" + "CCGGGGCCCC") * 3
    loose = "ACAACAACAAACAACAAACAACAAACAACAAACAACAAACAACAAACAACAAACAACAAACAACAA"
    a_struct = access.accessibility(structured, {"mid": (20, 40)}, ulength=8, window=60)["mid"]["mean"]
    a_loose = access.accessibility(loose, {"mid": (20, 40)}, ulength=8, window=60)["mid"]["mean"]
    assert a_loose > a_struct
