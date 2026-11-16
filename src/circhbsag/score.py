"""Skor kelayakan komposit Tahap 0 dan keputusan GO / NO-GO.

Struktur keputusan dua lapis:

  1. GERBANG KERAS (hard gates). Pelanggaran satu gerbang keras memaksa NO-GO
     tanpa memandang skor, karena kegagalannya tidak bisa dikompensasi kriteria
     lain (mis. situs BsaI internal membuat perakitan Golden Gate gagal).

  2. SKOR KOMPOSIT berbobot dari tujuh kriteria kontinu, masing-masing
     dinormalkan ke [0, 1] dengan fungsi yang didokumentasikan di bawah.

Semua bobot dan ambang ada di satu tempat (DEFAULT_WEIGHTS, DEFAULT_THRESHOLDS)
dan dapat ditimpa lewat config/config.yaml, sehingga sensitivitas keputusan
terhadap pilihan bobot dapat diuji (lihat tests/test_score.py dan
`circhbsag.score.sensitivity`).

Nilai yang dipakai pada normalisasi ditandai:
  [LITERATUR] bila berasal dari pustaka
  [PERKIRAAN] bila ditetapkan pipeline ini sebagai ambang kerja yang dapat diubah
"""
from __future__ import annotations

from dataclasses import dataclass, field

from . import access, fold, motifs
from .seqio import gc, longest_homopolymer

# --------------------------------------------------------------------------
# Bobot dan ambang
# --------------------------------------------------------------------------

DEFAULT_WEIGHTS: dict[str, float] = {
    # Aksesibilitas situs sambungan + P1/P10: penentu langsung sirkularisasi.
    "junction_accessibility": 0.25,
    # Integritas struktur IRES di dalam konstruk penuh.
    "ires_integrity": 0.20,
    # Beban dsRNA beruntun (imunogenisitas intrinsik, aktivasi PKR/RIG-I).
    "dsrna_burden": 0.15,
    # Aksesibilitas start codon ORF.
    "start_codon_accessibility": 0.15,
    # Kestabilan dan kekokohan lipatan lingkaran (MFE/nt, keragaman ensemble).
    "circle_stability": 0.10,
    # Kebersihan sekuens (situs enzim, GC jendela, homopolimer, terminator).
    "sequence_cleanliness": 0.10,
    # Beban motif m6A/imunostimulan.
    "motif_burden": 0.05,
}

DEFAULT_THRESHOLDS: dict[str, float] = {
    "go": 0.70,          # [PERKIRAAN] skor >= 0.70 dan tanpa gerbang keras gagal -> GO
    "conditional": 0.55,  # [PERKIRAAN] 0.55 <= skor < 0.70 -> GO BERSYARAT (optimasi dulu)
    # Parameter normalisasi
    "junction_unpaired_target": 0.50,   # [PERKIRAAN] P(tak-berpasangan) ideal di jendela BSJ
    "ires_identity_floor": 0.70,        # [PERKIRAAN] identitas status lipatan IRES minimum
    "ires_identity_target": 0.90,       # [PERKIRAAN] dianggap utuh penuh
    "dsrna_max_bp": 30,                 # [LITERATUR] dsRNA beruntun >~30 bp memicu PKR/RIG-I
    "dsrna_hard_bp": 45,                # [PERKIRAAN] gerbang keras
    "start_unpaired_target": 0.35,      # [PERKIRAAN] aksesibilitas start codon memadai
    "mfe_per_nt_lo": -0.50,             # [PERKIRAAN] lingkaran terlalu longgar bila > -0.25
    "mfe_per_nt_hi": -0.25,
    "ensemble_diversity_per_nt_max": 0.12,  # [PERKIRAAN] ensemble terlalu beragam = lipatan tak kokoh
    "m6a_per_100nt_max": 2.0,           # [PERKIRAAN]
    "gc_lo": 0.35,
    "gc_hi": 0.70,
    "homopolymer_max": 6,               # [PERKIRAAN] sintesis gen umumnya membatasi <= 6-7
    "max_circle_nt": 3000,              # [PERKIRAAN] batas praktis dari draf TA (<= 3,0 kb)
    # Jendela hulu kodon awal yang diperiksa untuk AUG sekerangka. IRES
    # pikornavirus memuat banyak AUG di dalam domainnya, sehingga pemeriksaan
    # dibatasi ke wilayah yang mekanistik relevan tepat di hulu AUG.
    "upstream_atg_window_nt": 60,       # [PERKIRAAN]
}

HARD_GATES = (
    "no_internal_typeIIS",
    "no_premature_T7_terminator",
    "orf_intact_and_unique_start",
    "dsrna_below_hard_limit",
    "circle_length_within_limit",
)


def _ramp(x: float, lo: float, hi: float) -> float:
    """Normalisasi linear naik: <=lo -> 0, >=hi -> 1."""
    if hi == lo:
        return 1.0 if x >= hi else 0.0
    return max(0.0, min(1.0, (x - lo) / (hi - lo)))


def _ramp_down(x: float, lo: float, hi: float) -> float:
    """Normalisasi linear turun: <=lo -> 1, >=hi -> 0."""
    return 1.0 - _ramp(x, lo, hi)


def _band(x: float, lo: float, hi: float, slack: float) -> float:
    """1 di dalam [lo, hi], turun linear sampai 0 pada jarak `slack` di luar band."""
    if lo <= x <= hi:
        return 1.0
    d = (lo - x) if x < lo else (x - hi)
    return max(0.0, 1.0 - d / slack) if slack > 0 else 0.0


@dataclass
class Criterion:
    name: str
    weight: float
    raw: dict
    subscore: float
    rationale: str

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "weight": self.weight,
            "subscore": round(self.subscore, 4),
            "weighted": round(self.weight * self.subscore, 4),
            "raw": self.raw,
            "rationale": self.rationale,
        }


@dataclass
class Decision:
    verdict: str
    composite: float
    criteria: list[Criterion]
    hard_gates: dict[str, dict]
    weights: dict[str, float] = field(default_factory=lambda: dict(DEFAULT_WEIGHTS))
    thresholds: dict[str, float] = field(default_factory=lambda: dict(DEFAULT_THRESHOLDS))

    def to_dict(self) -> dict:
        return {
            "verdict": self.verdict,
            "composite_score": round(self.composite, 4),
            "thresholds": {"go": self.thresholds["go"], "conditional": self.thresholds["conditional"]},
            "hard_gates": self.hard_gates,
            "hard_gates_failed": [k for k, v in self.hard_gates.items() if not v["passed"]],
            "criteria": [c.to_dict() for c in self.criteria],
            "weights": self.weights,
            "all_thresholds": self.thresholds,
        }


# --------------------------------------------------------------------------
# Pengukuran Tahap 0
# --------------------------------------------------------------------------


def measure(
    asm,
    *,
    temperatures: tuple[float, ...] = (37.0, 55.0),
    thresholds: dict | None = None,
    plfold_window: int = 80,
    plfold_ulength: int = 10,
) -> dict:
    """Jalankan seluruh pengukuran in silico Tahap 0 pada satu Assembly.

    Mengembalikan dict pengukuran mentah (tanpa keputusan) agar dapat
    disimpan, diperiksa ulang, dan diberi skor dengan bobot berbeda.
    """
    th = {**DEFAULT_THRESHOLDS, **(thresholds or {})}
    circle = asm.circle
    precursor = asm.precursor

    folds = {
        f"circle_{t:g}C": fold.fold(circle, circular=True, temperature_c=t) for t in temperatures
    }
    folds.update(
        {f"precursor_{t:g}C": fold.fold(precursor, circular=False, temperature_c=t) for t in temperatures}
    )

    # --- aksesibilitas sambungan balik pada lingkaran
    junction = {
        k: fold.junction_pairing(v, asm.bsj_position, span=10)
        for k, v in folds.items()
        if k.startswith("circle_")
    }

    # --- aksesibilitas situs splice dan P1/P10 pada prekursor (RNAplfold)
    regions: dict[str, tuple[int, int]] = {}
    for name in ("E1", "E2", "intron_5half", "intron_3half", "IRES", "ORF", "HA5", "HA3"):
        if name in asm.precursor_coords:
            regions[name] = asm.precursor_coords[name]
    # Situs splice 5' = batas E1/intron_5half; situs splice 3' = batas intron_3half/E2.
    if "E1" in asm.precursor_coords:
        e1a, e1b = asm.precursor_coords["E1"]
        regions["splice5_site"] = (max(1, e1b - 9), min(len(precursor), e1b + 10))
    if "E2" in asm.precursor_coords:
        e2a, _ = asm.precursor_coords["E2"]
        regions["splice3_site"] = (max(1, e2a - 10), min(len(precursor), e2a + 9))
    for key in ("intron_5half:IGS", "intron_3half:IGS"):
        if key in asm.precursor_coords:
            regions[key.replace(":", "_")] = asm.precursor_coords[key]

    acc = {
        f"{t:g}C": access.accessibility(
            precursor,
            regions,
            ulength=plfold_ulength,
            window=plfold_window,
            temperature_c=t,
        )
        for t in temperatures
    }

    # --- integritas IRES: bandingkan status lipatan IRES sendiri vs dalam konstruk
    ires_integrity = {}
    if "IRES" in asm.circle_coords:
        ia, ib = asm.circle_coords["IRES"]
        ires_seq = circle[ia - 1 : ib]
        for t in temperatures:
            alone = fold.fold(ires_seq, circular=False, temperature_c=t)
            in_circle = folds[f"circle_{t:g}C"]
            ident = fold.subsequence_structure_identity(
                in_circle.mfe_structure, alone.mfe_structure, ia - 1
            )
            ires_integrity[f"{t:g}C"] = {
                "identity": round(ident, 4),
                "ires_alone_mfe_kcal": round(alone.mfe_kcal, 3),
                "ires_alone_mfe_per_nt": round(alone.mfe_per_nt, 5),
                "ires_length": len(ires_seq),
            }

    # --- beban dsRNA pada lingkaran dan prekursor
    dsrna = {
        k: {
            **access.longest_dsrna(v.mfe_structure),
            "violations": access.dsrna_violations(v.mfe_structure, int(th["dsrna_max_bp"])),
        }
        for k, v in folds.items()
    }

    # --- aksesibilitas start codon ORF (pada lingkaran)
    start = {}
    if "ORF" in asm.circle_coords:
        oa, _ = asm.circle_coords["ORF"]
        win = (max(1, oa - 6), min(len(circle), oa + 14))
        for t in temperatures:
            a = access.accessibility(
                circle, {"start_codon": win}, ulength=plfold_ulength, window=plfold_window, temperature_c=t
            )
            jp = fold.unpaired_probability(folds[f"circle_{t:g}C"], list(range(win[0], win[1] + 1)))
            start[f"{t:g}C"] = {
                "window": list(win),
                "plfold": a["start_codon"],
                "mean_unpaired_global_fold": round(sum(jp.values()) / len(jp), 5),
                "structure": folds[f"circle_{t:g}C"].mfe_structure[win[0] - 1 : win[1]],
            }

    # --- pemindaian motif pada prekursor (relevan untuk IVT/kloning) dan lingkaran
    motif_scan = {
        "precursor": motifs.summarize(precursor),
        "circle": motifs.summarize(circle),
    }
    orf_seq = asm.subseq("circle", "ORF") if "ORF" in asm.circle_coords else ""
    orf_checks = {}
    if orf_seq:
        orfs = __import__("circhbsag.seqio", fromlist=["find_orfs"]).find_orfs(circle * 2, min_aa=50)
        orf_checks = {
            "orf_length_nt": len(orf_seq),
            "orf_starts_with_atg": orf_seq[:3] == "ATG",
            "orf_in_frame_length_ok": len(orf_seq) % 3 == 0,
            "upstream_in_frame_atg_near": motifs.upstream_atg_in_frame(
                circle, asm.circle_coords["ORF"][0] - 1, within=int(th["upstream_atg_window_nt"])
            ),
            "upstream_out_of_frame_atg_near": motifs.upstream_atg_out_of_frame(
                circle, asm.circle_coords["ORF"][0] - 1, within=30
            ),
            "upstream_in_frame_atg_all": motifs.upstream_atg_in_frame(
                circle, asm.circle_coords["ORF"][0] - 1
            ),
            "longest_orf_on_doubled_circle_aa": orfs[0]["aa"] if orfs else 0,
        }

    return {
        "assembly": asm.summary(),
        "folds": {k: v.to_dict() for k, v in folds.items()},
        "junction": junction,
        "accessibility": acc,
        "ires_integrity": ires_integrity,
        "dsrna": dsrna,
        "start_codon": start,
        "motifs": motif_scan,
        "orf_checks": orf_checks,
        "gc": {"circle": round(gc(circle), 4), "precursor": round(gc(precursor), 4)},
        "homopolymer": {
            "circle": longest_homopolymer(circle),
            "precursor": longest_homopolymer(precursor),
        },
        "thresholds_used": th,
        "temperatures_c": list(temperatures),
    }


# --------------------------------------------------------------------------
# Penilaian
# --------------------------------------------------------------------------


def score(m: dict, *, weights: dict | None = None, thresholds: dict | None = None) -> Decision:
    """Hitung gerbang keras, subskor tiap kriteria, dan keputusan GO/NO-GO."""
    w = {**DEFAULT_WEIGHTS, **(weights or {})}
    th = {**DEFAULT_THRESHOLDS, **(thresholds or {}), **m.get("thresholds_used", {})}
    if thresholds:
        th.update(thresholds)

    temps = m["temperatures_c"]
    worst = lambda vals: min(vals) if vals else 0.0  # noqa: E731

    # ---------------- gerbang keras ----------------
    pre = m["motifs"]["precursor"]
    n_typeIIS = pre["restriction"].get("BsaI", 0) + pre["restriction"].get("BsmBI_Esp3I", 0)
    gates: dict[str, dict] = {}
    gates["no_internal_typeIIS"] = {
        "passed": n_typeIIS == 0,
        "value": n_typeIIS,
        "limit": 0,
        "why": "situs BsaI/BsmBI internal akan dipotong saat perakitan Golden Gate",
    }
    n_term = len(pre["t7_terminator"])
    gates["no_premature_T7_terminator"] = {
        "passed": n_term == 0,
        "value": n_term,
        "limit": 0,
        "why": "terminator T7 prematur memotong transkrip prekursor sebelum intron 5'",
    }
    oc = m.get("orf_checks", {})
    orf_ok = bool(oc) and oc.get("orf_starts_with_atg") and oc.get("orf_in_frame_length_ok")
    orf_ok = bool(orf_ok) and not oc.get("upstream_in_frame_atg_near")
    gates["orf_intact_and_unique_start"] = {
        "passed": bool(orf_ok),
        "value": {
            "starts_with_atg": oc.get("orf_starts_with_atg"),
            "length_mod3_zero": oc.get("orf_in_frame_length_ok"),
            "upstream_in_frame_atg_near": oc.get("upstream_in_frame_atg_near"),
            "window_nt": th["upstream_atg_window_nt"],
            "upstream_out_of_frame_atg_near_count": len(
                oc.get("upstream_out_of_frame_atg_near") or []
            ),
            "upstream_in_frame_atg_anywhere_count": len(
                oc.get("upstream_in_frame_atg_all") or []
            ),
        },
        "limit": (
            f"ATG di awal, panjang kelipatan 3, tanpa ATG hulu sekerangka dalam "
            f"{th['upstream_atg_window_nt']:.0f} nt terakhir"
        ),
        "why": (
            "ATG sekerangka tepat di hulu kodon awal memperpanjang ujung-N produk. "
            "AUG di dalam domain IRES dihitung terpisah sebagai informasi, bukan "
            "pelanggaran, karena IRES pikornavirus memang memuat banyak AUG."
        ),
    }
    ds_max = max(m["dsrna"][k]["longest_bp"] for k in m["dsrna"])
    gates["dsrna_below_hard_limit"] = {
        "passed": ds_max < th["dsrna_hard_bp"],
        "value": ds_max,
        "limit": th["dsrna_hard_bp"],
        "why": "dsRNA sangat panjang memicu PKR/RIG-I dan menekan translasi",
    }
    clen = m["assembly"]["circle_length"]
    gates["circle_length_within_limit"] = {
        "passed": clen <= th["max_circle_nt"],
        "value": clen,
        "limit": th["max_circle_nt"],
        "why": "batas praktis sintesis gen dan efisiensi sirkularisasi",
    }

    # ---------------- kriteria berbobot ----------------
    crit: list[Criterion] = []

    # 1. aksesibilitas sambungan: gabungan BSJ pada lingkaran + situs splice/P1/P10 pada prekursor
    bsj_vals = [m["junction"][f"circle_{t:g}C"]["mean_unpaired_probability"] for t in temps]
    splice_vals = []
    for t in temps:
        a = m["accessibility"][f"{t:g}C"]
        for key in ("splice5_site", "splice3_site", "intron_5half_IGS", "intron_3half_IGS"):
            if key in a and a[key]["mean"] is not None:
                splice_vals.append(a[key]["mean"])
    bsj_sub = _ramp(worst(bsj_vals), 0.05, th["junction_unpaired_target"])
    splice_sub = _ramp(worst(splice_vals), 0.01, 0.30) if splice_vals else 0.5
    sub = 0.6 * bsj_sub + 0.4 * splice_sub
    crit.append(
        Criterion(
            "junction_accessibility",
            w["junction_accessibility"],
            {
                "bsj_mean_unpaired_per_temp": dict(zip([f"{t:g}C" for t in temps], bsj_vals, strict=False)),
                "bsj_worst": round(worst(bsj_vals), 5),
                "splice_igs_worst_accessibility": round(worst(splice_vals), 6) if splice_vals else None,
                "splice_regions_measured": len(splice_vals),
            },
            sub,
            "60% dari P(tak-berpasangan) jendela BSJ pada lingkaran (suhu terburuk), "
            "40% dari aksesibilitas RNAplfold situs splice 5'/3' dan IGS (P1/P10) pada prekursor. "
            "Bila anotasi IGS tidak tersedia, komponen kedua diberi nilai netral 0,5.",
        )
    )

    # 2. integritas IRES
    if m["ires_integrity"]:
        idents = [v["identity"] for v in m["ires_integrity"].values()]
        sub = _ramp(worst(idents), th["ires_identity_floor"], th["ires_identity_target"])
        raw = {"identity_per_temp": {k: v["identity"] for k, v in m["ires_integrity"].items()},
               "worst_identity": round(worst(idents), 4)}
        rat = (
            "Fraksi posisi IRES yang mempertahankan status berpasangan/tak-berpasangan "
            f"sama seperti saat IRES dilipat sendiri. 0 pada <= {th['ires_identity_floor']}, "
            f"1 pada >= {th['ires_identity_target']}. Proksi struktural, bukan uji fungsi translasi."
        )
    else:
        sub, raw, rat = 0.0, {"error": "wilayah IRES tidak dianotasi"}, "IRES tidak ada pada konstruk"
    crit.append(Criterion("ires_integrity", w["ires_integrity"], raw, sub, rat))

    # 3. beban dsRNA
    sub = _ramp_down(ds_max, th["dsrna_max_bp"] * 0.6, th["dsrna_hard_bp"])
    crit.append(
        Criterion(
            "dsrna_burden",
            w["dsrna_burden"],
            {
                "longest_bp_any_fold": ds_max,
                "per_fold": {k: v["longest_bp"] for k, v in m["dsrna"].items()},
                "n_violations_over_threshold": {k: len(v["violations"]) for k, v in m["dsrna"].items()},
            },
            sub,
            f"1 bila heliks beruntun terpanjang <= {th['dsrna_max_bp'] * 0.6:.0f} bp, "
            f"0 pada >= {th['dsrna_hard_bp']:.0f} bp (gerbang keras).",
        )
    )

    # 4. aksesibilitas start codon
    if m["start_codon"]:
        svals = [v["plfold"]["mean"] for v in m["start_codon"].values() if v["plfold"]["mean"] is not None]
        gvals = [v["mean_unpaired_global_fold"] for v in m["start_codon"].values()]
        sub = 0.5 * _ramp(worst(svals), 0.01, th["start_unpaired_target"]) + 0.5 * _ramp(
            worst(gvals), 0.05, 0.60
        )
        raw = {
            "plfold_mean_per_temp": {k: v["plfold"]["mean"] for k, v in m["start_codon"].items()},
            "global_unpaired_per_temp": {k: v["mean_unpaired_global_fold"] for k, v in m["start_codon"].items()},
        }
    else:
        sub, raw = 0.0, {"error": "ORF tidak dianotasi"}
    crit.append(
        Criterion(
            "start_codon_accessibility",
            w["start_codon_accessibility"],
            raw,
            sub,
            "Rerata dari aksesibilitas lokal RNAplfold dan P(tak-berpasangan) lipatan global "
            "pada jendela -6..+14 relatif AUG; keduanya diambil pada suhu terburuk.",
        )
    )

    # 5. kestabilan lingkaran
    mfe_pnt = [m["folds"][f"circle_{t:g}C"]["mfe_per_nt"] for t in temps]
    div_pnt = [
        m["folds"][f"circle_{t:g}C"]["ensemble_diversity"] / m["assembly"]["circle_length"]
        for t in temps
    ]
    band_sub = min(_band(v, th["mfe_per_nt_lo"], th["mfe_per_nt_hi"], 0.15) for v in mfe_pnt)
    div_sub = min(_ramp_down(v, 0.0, th["ensemble_diversity_per_nt_max"]) for v in div_pnt)
    sub = 0.5 * band_sub + 0.5 * div_sub
    crit.append(
        Criterion(
            "circle_stability",
            w["circle_stability"],
            {
                "mfe_per_nt_per_temp": dict(zip([f"{t:g}C" for t in temps], [round(v, 5) for v in mfe_pnt], strict=False)),
                "ensemble_diversity_per_nt": dict(
                    zip([f"{t:g}C" for t in temps], [round(v, 5) for v in div_pnt], strict=False)
                ),
            },
            sub,
            f"Separuh dari MFE/nt berada dalam band [{th['mfe_per_nt_lo']}, {th['mfe_per_nt_hi']}] "
            "(terlalu longgar = lingkaran tidak kompak; terlalu rapat = risiko dsRNA dan "
            "situs tertutup), separuh dari keragaman ensemble per nt (lipatan kokoh).",
        )
    )

    # 6. kebersihan sekuens
    gcv = len(pre["gc_violations"])
    hp = m["homopolymer"]["precursor"][1]
    gcc = m["gc"]["circle"]
    sub = (
        0.4 * (1.0 if n_typeIIS == 0 else 0.0)
        + 0.2 * _ramp_down(gcv, 0, 10)
        + 0.2 * _ramp_down(hp, th["homopolymer_max"], th["homopolymer_max"] + 4)
        + 0.2 * _band(gcc, th["gc_lo"], th["gc_hi"], 0.10)
    )
    crit.append(
        Criterion(
            "sequence_cleanliness",
            w["sequence_cleanliness"],
            {
                "typeIIS_sites": n_typeIIS,
                "gc_window_violations": gcv,
                "longest_homopolymer": hp,
                "circle_gc": gcc,
                "t7_terminator_hits": n_term,
            },
            sub,
            "40% tanpa situs Tipe IIS, 20% jendela GC 30-70%, 20% homopolimer pendek, "
            "20% GC lingkaran dalam band yang dapat disintesis.",
        )
    )

    # 7. beban motif m6A/imunostimulan
    m6a = m["motifs"]["circle"]["m6a_density_per_100nt"]
    n_avoid = sum(m["motifs"]["circle"]["avoid_motifs"].values())
    sub = 0.6 * _ramp_down(m6a, th["m6a_per_100nt_max"] * 0.5, th["m6a_per_100nt_max"] * 2) + 0.4 * _ramp_down(
        n_avoid / max(1.0, m["assembly"]["circle_length"] / 1000), 0, 10
    )
    crit.append(
        Criterion(
            "motif_burden",
            w["motif_burden"],
            {
                "m6a_per_100nt": m6a,
                "m6a_sites": m["motifs"]["circle"]["m6a_sites"],
                "avoid_motif_hits": m["motifs"]["circle"]["avoid_motifs"],
            },
            sub,
            "60% kerapatan RRACH (m6A) dan 40% motif yang dihindari per kb. "
            "Catatan: bila penyisipan m6A 1% direncanakan sebagai fitur, kriteria ini "
            "harus dibobot ulang, bukan diminimalkan.",
        )
    )

    total_w = sum(c.weight for c in crit)
    composite = sum(c.weight * c.subscore for c in crit) / total_w if total_w else 0.0

    failed = [k for k, v in gates.items() if not v["passed"]]
    if failed:
        verdict = "NO-GO"
    elif composite >= th["go"]:
        verdict = "GO"
    elif composite >= th["conditional"]:
        verdict = "GO-BERSYARAT"
    else:
        verdict = "NO-GO"

    return Decision(verdict=verdict, composite=composite, criteria=crit, hard_gates=gates,
                    weights=w, thresholds=th)


def sensitivity(m: dict, *, n: int = 200, seed: int = 7, jitter: float = 0.5) -> dict:
    """Uji sensitivitas keputusan terhadap bobot (Dirichlet jitter di sekitar bobot default).

    Melaporkan distribusi skor komposit dan fraksi keputusan GO, sehingga
    kesimpulan GO/NO-GO tidak bergantung pada satu set bobot saja.
    """
    import random as _r

    rng = _r.Random(seed)
    base = DEFAULT_WEIGHTS
    scores, verdicts = [], {}
    for _ in range(n):
        wts = {}
        for k, v in base.items():
            lo, hi = v * (1 - jitter), v * (1 + jitter)
            wts[k] = rng.uniform(lo, hi)
        d = score(m, weights=wts)
        scores.append(d.composite)
        verdicts[d.verdict] = verdicts.get(d.verdict, 0) + 1
    scores.sort()
    return {
        "n": n,
        "jitter_fraction": jitter,
        "composite_min": round(scores[0], 4),
        "composite_p05": round(scores[int(0.05 * n)], 4),
        "composite_median": round(scores[n // 2], 4),
        "composite_p95": round(scores[min(n - 1, int(0.95 * n))], 4),
        "composite_max": round(scores[-1], 4),
        "verdict_counts": verdicts,
    }
