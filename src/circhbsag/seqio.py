"""Pembacaan, validasi, dan utilitas dasar sekuens.

Semua fungsi di modul ini deterministik dan tidak mengambil data dari internet.
Sekuens referensi WAJIB disediakan pengguna (lihat data/reference/README.md);
paket ini tidak pernah mengarang sekuens.
"""
from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path

DNA = set("ACGT")
RNA = set("ACGU")

AA1 = (
    "ACDEFGHIKLMNPQRSTVWY"
)

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

SYN = {}
for _c, _a in CODON_TABLE.items():
    SYN.setdefault(_a, []).append(_c)


class SequenceError(ValueError):
    """Sekuens tidak valid untuk analisis."""


@dataclass
class Record:
    """Satu rekaman FASTA."""

    seq_id: str
    seq: str
    description: str = ""
    meta: dict = field(default_factory=dict)

    def __len__(self) -> int:  # pragma: no cover - trivial
        return len(self.seq)

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.seq.encode()).hexdigest()


def clean(seq: str, *, alphabet: str = "dna") -> str:
    """Normalkan sekuens: huruf besar, buang spasi/nomor, U<->T sesuai alfabet."""
    s = re.sub(r"[\s\d]", "", seq).upper()
    if alphabet == "dna":
        s = s.replace("U", "T")
        allowed = DNA
    elif alphabet == "rna":
        s = s.replace("T", "U")
        allowed = RNA
    else:  # pragma: no cover
        raise ValueError(f"alfabet tidak dikenal: {alphabet}")
    bad = sorted(set(s) - allowed)
    if bad:
        raise SequenceError(
            f"karakter non-{alphabet.upper()} ditemukan: {bad[:10]}"
        )
    if not s:
        raise SequenceError("sekuens kosong")
    return s


def read_fasta(path: str | Path, *, alphabet: str = "dna") -> list[Record]:
    """Baca FASTA (multi-rekaman). Tidak toleran terhadap karakter ambigu."""
    path = Path(path)
    text = path.read_text()
    records: list[Record] = []
    header, chunks = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append(_mk(header, chunks, alphabet))
            header, chunks = line[1:].strip(), []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        records.append(_mk(header, chunks, alphabet))
    if not records:
        raise SequenceError(f"tidak ada rekaman FASTA di {path}")
    return records


def _mk(header: str, chunks: list[str], alphabet: str) -> Record:
    parts = header.split(None, 1)
    seq_id = parts[0] if parts else "unnamed"
    desc = parts[1] if len(parts) > 1 else ""
    return Record(seq_id=seq_id, seq=clean("".join(chunks), alphabet=alphabet), description=desc)


def write_fasta(records: list[Record], path: str | Path, width: int = 60) -> None:
    lines = []
    for r in records:
        lines.append(f">{r.seq_id} {r.description}".rstrip())
        lines += [r.seq[i : i + width] for i in range(0, len(r.seq), width)]
    Path(path).write_text("\n".join(lines) + "\n")


def to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def to_dna(seq: str) -> str:
    return seq.upper().replace("U", "T")


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTUacgtu", "TGCAATGCAA")
    return seq.translate(table)[::-1]


def gc(seq: str) -> float:
    """Fraksi GC (0-1)."""
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s) if s else 0.0


def gc_windows(seq: str, window: int = 50, step: int = 10) -> list[tuple[int, float]]:
    """GC lokal per jendela; mengembalikan (posisi_awal_0based, fraksi_GC)."""
    out = []
    if len(seq) < window:
        return [(0, gc(seq))]
    for i in range(0, len(seq) - window + 1, step):
        out.append((i, gc(seq[i : i + window])))
    return out


def translate(seq: str, *, stop_at_stop: bool = True) -> str:
    """Translasi kerangka 1. Panjang bukan kelipatan 3 memicu SequenceError."""
    s = to_dna(clean(seq))
    if len(s) % 3:
        raise SequenceError(f"panjang {len(s)} bukan kelipatan 3")
    aa = []
    for i in range(0, len(s), 3):
        res = CODON_TABLE[s[i : i + 3]]
        if res == "*" and stop_at_stop:
            aa.append("*")
            break
        aa.append(res)
    return "".join(aa)


def rotate(seq: str, k: int) -> str:
    """Rotasi sirkular sekuens sebesar k nt (untuk templat terpermutasi)."""
    if not seq:
        return seq
    k %= len(seq)
    return seq[k:] + seq[:k]


def longest_homopolymer(seq: str) -> tuple[str, int]:
    best_base, best = "", 0
    for m in re.finditer(r"(.)\1*", seq):
        if len(m.group(0)) > best:
            best, best_base = len(m.group(0)), m.group(1)
    return best_base, best


def find_orfs(seq: str, min_aa: int = 30) -> list[dict]:
    """Cari ORF ATG..stop di tiga kerangka untai maju (koordinat 0-based)."""
    s = to_dna(seq)
    orfs = []
    for frame in range(3):
        i = frame
        while i < len(s) - 2:
            if s[i : i + 3] == "ATG":
                j = i
                while j < len(s) - 2:
                    if CODON_TABLE.get(s[j : j + 3]) == "*":
                        if (j - i) // 3 >= min_aa:
                            orfs.append({"start": i, "end": j + 3, "frame": frame, "aa": (j - i) // 3})
                        break
                    j += 3
                i = j + 3 if j > i else i + 3
            else:
                i += 3
    return sorted(orfs, key=lambda d: -d["aa"])
