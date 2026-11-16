"""Tes utilitas sekuens."""
from __future__ import annotations

import pytest

from circhbsag import seqio


def test_clean_dna_and_rna():
    assert seqio.clean("acg u\n12") == "ACGT"
    assert seqio.clean("ACGT", alphabet="rna") == "ACGU"


def test_clean_rejects_invalid():
    with pytest.raises(seqio.SequenceError):
        seqio.clean("ACGNT")
    with pytest.raises(seqio.SequenceError):
        seqio.clean("")


def test_revcomp_roundtrip():
    s = "ACGTTGCA"
    assert seqio.revcomp(seqio.revcomp(s)) == s
    assert seqio.revcomp("ACGT") == "ACGT"


def test_gc_and_windows():
    assert seqio.gc("GCGC") == 1.0
    assert seqio.gc("ATAT") == 0.0
    w = seqio.gc_windows("GC" * 50, window=20, step=20)
    assert all(v == 1.0 for _, v in w)


def test_translate_and_errors():
    assert seqio.translate("ATGGCCTGA") == "MA*"
    with pytest.raises(seqio.SequenceError):
        seqio.translate("ATGGC")


def test_rotate_is_circular():
    s = "ABCDE"
    assert seqio.rotate(s, 2) == "CDEAB"
    assert seqio.rotate(s, len(s)) == s
    assert seqio.rotate(s, len(s) + 1) == seqio.rotate(s, 1)


def test_longest_homopolymer():
    assert seqio.longest_homopolymer("AAACGGGGT") == ("G", 4)


def test_find_orfs_finds_longest_first():
    orf = "ATG" + "GCC" * 40 + "TGA"
    seq = "TTTT" + orf + "TTTT"
    hits = seqio.find_orfs(seq, min_aa=10)
    assert hits
    assert hits[0]["aa"] == 41
    assert seq[hits[0]["start"] : hits[0]["start"] + 3] == "ATG"


def test_fasta_roundtrip(tmp_path):
    recs = [seqio.Record("a", "ACGT", "desc satu"), seqio.Record("b", "TTTTGGGG")]
    p = tmp_path / "x.fa"
    seqio.write_fasta(recs, p)
    back = seqio.read_fasta(p)
    assert [r.seq_id for r in back] == ["a", "b"]
    assert [r.seq for r in back] == ["ACGT", "TTTTGGGG"]
    assert back[0].description == "desc satu"
    assert back[0].sha256 == recs[0].sha256
