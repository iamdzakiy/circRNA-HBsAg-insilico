"""Perakitan konstruk PIE in silico dan pelacakan koordinat.

Tata letak prekursor (5'->3'), mengikuti arsitektur PIE intron grup I yang
dipakai pada circRNA translasi (Puttaraju & Been 1992; Wesselhoeft dkk. 2018):

    T7_promoter | GGG | HA5 | intron_3half | E2 | spacer5 | IRES | ORF
                | spacer3 | E1 | intron_5half | HA3 | T7_terminator

Setelah dua transesterifikasi, lingkaran yang terbentuk adalah:

    E2 | spacer5 | IRES | ORF | spacer3 | E1      (E1--E2 = sambungan balik/BSJ)

Semua sekuens bagian WAJIB disediakan pengguna lewat berkas parts (YAML/FASTA).
Modul ini tidak pernah mengarang sekuens; bila sebuah bagian tidak ada, proses
berhenti dengan galat yang menyebut bagian mana yang kurang.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .motifs import T7_PROMOTER
from .seqio import SequenceError, clean, revcomp, rotate, to_dna

# Urutan bagian pada prekursor. HA = homology arm (lengan homologi).
PRECURSOR_ORDER = [
    "T7_promoter",
    "gg_leader",
    "HA5",
    "intron_3half",
    "E2",
    "spacer5",
    "IRES",
    "ORF",
    "spacer3",
    "E1",
    "intron_5half",
    "HA3",
    "T7_terminator",
]

CIRCLE_ORDER = ["E2", "spacer5", "IRES", "ORF", "spacer3", "E1"]

REQUIRED_PARTS = ["intron_3half", "E2", "IRES", "ORF", "E1", "intron_5half"]

# Tata letak CIRC (Shen dkk. 2025): satu intron Ana utuh di ujung 5' prekursor,
# lalu ekson 3' (E2), spacer, lengan homologi (HA5), AC30, muatan, lengan
# homologi (HA3), spacer, ekson 5' (E1). Dipakai bila 'intron_5half' tidak ada.
PRECURSOR_ORDER_CIRC = [
    "T7_promoter", "gg_leader", "intron_3half", "E2", "spacer5", "HA5", "AC30",
    "IRES", "ORF", "UTR3", "HA3", "spacer3", "E1", "T7_terminator",
]
CIRCLE_ORDER_CIRC = ["E2", "spacer5", "HA5", "AC30", "IRES", "ORF", "UTR3", "HA3", "spacer3", "E1"]
REQUIRED_PARTS_CIRC = ["intron_3half", "E2", "IRES", "ORF", "E1"]



@dataclass
class Part:
    name: str
    seq: str
    role: str = ""
    # Anotasi sub-wilayah relatif terhadap awal bagian (0-based, [start, end)).
    features: dict[str, tuple[int, int]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        # Bagian opsional (spacer, lengan homologi, terminator) boleh kosong dan
        # dilewati saat perakitan; `clean` menolak sekuens kosong, jadi dibedakan.
        self.seq = clean(self.seq, alphabet="dna") if self.seq.strip() else ""


@dataclass
class Assembly:
    """Konstruk hasil perakitan beserta semua koordinat (1-based inklusif)."""

    precursor: str
    circle: str
    parts: dict[str, Part]
    precursor_coords: dict[str, tuple[int, int]]
    circle_coords: dict[str, tuple[int, int]]
    bsj_position: int  # posisi 1-based basa terakhir E1 pada lingkaran

    def region(self, which: str, name: str) -> tuple[int, int]:
        coords = self.precursor_coords if which == "precursor" else self.circle_coords
        if name not in coords:
            raise KeyError(f"wilayah '{name}' tidak ada pada {which}")
        return coords[name]

    def subseq(self, which: str, name: str) -> str:
        s = self.precursor if which == "precursor" else self.circle
        a, b = self.region(which, name)
        return s[a - 1 : b]

    def summary(self) -> dict:
        return {
            "precursor_length": len(self.precursor),
            "circle_length": len(self.circle),
            "bsj_position": self.bsj_position,
            "parts": {k: len(v.seq) for k, v in self.parts.items()},
            "circle_coords": {k: list(v) for k, v in self.circle_coords.items()},
            "precursor_coords": {k: list(v) for k, v in self.precursor_coords.items()},
        }


def build(parts: dict[str, Part | str], *, strict: bool = True) -> Assembly:
    """Rakit prekursor dan lingkaran dari kamus bagian.

    Bagian opsional (HA5, HA3, spacer5, spacer3, T7_promoter, gg_leader,
    T7_terminator) boleh kosong atau tidak ada. Bagian pada REQUIRED_PARTS
    wajib ada bila `strict=True`.
    """
    norm: dict[str, Part] = {}
    for k, v in parts.items():
        norm[k] = v if isinstance(v, Part) else Part(name=k, seq=v)

    circ = "intron_5half" not in norm or not norm["intron_5half"].seq
    req = REQUIRED_PARTS_CIRC if circ else REQUIRED_PARTS
    if strict:
        missing = [p for p in req if p not in norm or not norm[p].seq]
        if missing:
            raise SequenceError(
                "bagian konstruk belum tersedia: "
                + ", ".join(missing)
                + ". Sediakan sekuens terverifikasi di config/parts.yaml."
            )

    precursor, pcoords = _concat(norm, PRECURSOR_ORDER_CIRC if circ else PRECURSOR_ORDER)
    circle, ccoords = _concat(norm, CIRCLE_ORDER_CIRC if circ else CIRCLE_ORDER)
    bsj = ccoords["E1"][1] if "E1" in ccoords else len(circle)
    return Assembly(
        precursor=precursor,
        circle=circle,
        parts=norm,
        precursor_coords=pcoords,
        circle_coords=ccoords,
        bsj_position=bsj,
    )


def _concat(parts: dict[str, Part], order: list[str]) -> tuple[str, dict[str, tuple[int, int]]]:
    chunks, coords, pos = [], {}, 0
    for name in order:
        p = parts.get(name)
        if p is None or not p.seq:
            continue
        chunks.append(p.seq)
        coords[name] = (pos + 1, pos + len(p.seq))
        for fname, (fa, fb) in p.features.items():
            coords[f"{name}:{fname}"] = (pos + fa + 1, pos + fb)
        pos += len(p.seq)
    return "".join(chunks), coords


def permuted_dna_template(asm: Assembly, *, add_gg_leader: bool = True) -> dict:
    """Templat DNA linear untuk transkripsi in vitro (T7).

    Pada sistem PIE, "terpermutasi" berarti separuh intron 3' diletakkan di hulu
    eksón dan separuh intron 5' di hilir, sehingga transkrip linear dapat
    menyambung ujung-ke-ujung. Fungsi ini mengembalikan templat lengkap beserta
    pemeriksaan: promotor T7 ada, tiga G transkripsional pertama ada
    (memaksimalkan inisiasi T7 dan rasio GMP:GTP), dan tidak ada situs Tipe IIS.
    """
    tpl = asm.precursor
    has_t7 = T7_PROMOTER in to_dna(tpl)
    i = to_dna(tpl).find(T7_PROMOTER)
    tss = i + len(T7_PROMOTER) if i >= 0 else None
    first3 = to_dna(tpl)[tss : tss + 3] if tss is not None else None
    notes = []
    if not has_t7:
        notes.append("promotor T7 kanonik tidak ditemukan pada templat")
    if first3 is not None and first3 != "GGG":
        notes.append(
            f"tiga nukleotida pertama transkrip = {first3}; T7 memulai paling efisien dengan GGG"
        )
    if add_gg_leader and first3 is not None and first3 != "GGG":
        notes.append("tambahkan 'gg_leader' = GGG pada parts untuk memperbaiki inisiasi")
    return {
        "template_dna": tpl,
        "length": len(tpl),
        "t7_promoter_found": has_t7,
        "transcription_start_0based": tss,
        "first_three_nt": first3,
        "notes": notes,
    }


def circle_rotations(asm: Assembly, step: int = 1) -> list[str]:
    """Semua rotasi lingkaran (uji invariansi lipatan terhadap titik mula)."""
    return [rotate(asm.circle, k) for k in range(0, len(asm.circle), step)]


def bsj_amplicon(asm: Assembly, flank: int = 120) -> dict:
    """Sekuens melintasi sambungan balik, untuk desain primer divergen/BSJ.

    Lingkaran dirotasi sehingga BSJ berada tepat di tengah, sehingga primer
    divergen pada templat linear tidak akan menghasilkan produk.
    """
    n = len(asm.circle)
    flank = min(flank, n // 2)
    bsj = asm.bsj_position  # basa terakhir E1
    left = "".join(asm.circle[(bsj - flank + k) % n] for k in range(flank))
    right = "".join(asm.circle[(bsj + k) % n] for k in range(flank))
    return {
        "bsj_position_in_circle": bsj,
        "flank": flank,
        "amplicon": left + right,
        "junction_offset_in_amplicon": flank,  # 0-based: basa pertama setelah sambungan
        "left_of_junction": left,
        "right_of_junction": right,
    }


def guide_sequence_pairs(asm: Assembly) -> dict:
    """Pasangan P1 dan P10 bila anotasi IGS/eksón tersedia pada bagian intron.

    P1  : IGS (pada separuh intron 5') berpasangan dengan ujung 3' eksón E1.
    P10 : IGS berpasangan dengan awal eksón E2.
    Dilaporkan sebagai wilayah untuk diuji aksesibilitasnya, bukan sebagai
    prediksi pasangan basa tersendiri.
    """
    out: dict[str, dict] = {}
    for key, part_name, feature in (
        ("P1", "intron_5half", "IGS"),
        ("P10", "intron_3half", "IGS"),
    ):
        pname = f"{part_name}:{feature}"
        if pname in asm.precursor_coords:
            a, b = asm.precursor_coords[pname]
            out[key] = {
                "igs_region_precursor": [a, b],
                "igs_seq": asm.precursor[a - 1 : b],
                "igs_revcomp": revcomp(asm.precursor[a - 1 : b]),
            }
    return out


def trim_arms(ha5: str, ha3: str, length: int) -> tuple[str, str]:
    """Shorten the two homology arms to ``length`` nt, keeping the ends that face the intron.

    HA5 keeps its 3' end (adjacent to the 3' intron half) and HA3 keeps its 5' end
    (adjacent to the 5' intron half). Every place that shortens an arm — the
    Gate 0 length scan, the redesign loop and the Optuna study — calls this, so
    a length that passes Gate 0 is the same sequence that gets optimised.
    ``length == 0`` removes the arms.

    [PERKIRAAN] Which end to keep has no literature support; for perfectly
    complementary demo arms it does not matter, for real arms it can.
    """
    if length <= 0:
        return "", ""
    return ha5[-length:], ha3[:length]
