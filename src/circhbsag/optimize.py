"""Optimasi multi-objektif desain konstruk dengan Optuna (lisensi: MIT).

Variabel keputusan
------------------
  spacer5_len, spacer3_len   : panjang spacer antara eksón dan muatan (nt)
  spacer_gc_target           : target GC spacer yang dibangkitkan
  spacer_motif               : pilihan pustaka spacer (unstructured / polyA-rich / poly-purine)
  ha_len                     : panjang lengan homologi yang dipakai (dipangkas dari ujung)
  codon_seed                 : seed optimizer kodon (menghasilkan ORF sinonim berbeda)
  codon_w_gc, codon_w_m6a    : bobot objektif optimasi kodon

Objektif (semua diubah menjadi "maximize")
------------------------------------------
  f1 = P(tak-berpasangan) jendela BSJ pada lingkaran, suhu terburuk
  f2 = aksesibilitas situs splice/IGS pada prekursor, suhu terburuk
  f3 = identitas struktur IRES (utuh), suhu terburuk
  f4 = -panjang dsRNA beruntun terpanjang
  f5 = aksesibilitas start codon
  f6 = -|GC_lingkaran - target|

Dipakai NSGA-II (Optuna `NSGAIISampler`) sehingga hasilnya adalah front Pareto,
bukan satu titik; pemilihan akhir memakai skor komposit Tahap 0 (score.score).
Seluruh proses memakai seed tetap agar reprodusibel.
"""
from __future__ import annotations

import random
from dataclasses import dataclass

from . import codon, construct, score
from .seqio import gc as gc_frac

SPACER_LIBRARIES = {
    # Pustaka spacer. Pilihan jenisnya mengikuti [LITERATUR]:
    #
    #   Wesselhoeft dkk. (2018, Nat Commun 9:2629) menunjukkan bahwa spacer yang
    #   diletakkan antara situs splice PIE 3' dan IRES menaikkan efisiensi
    #   splicing dari 46% menjadi 87%, sedangkan spacer yang dirancang MERUSAK
    #   struktur kedua separuh intron menghilangkan splicing sama sekali.
    #   Prinsipnya: spacer harus membiarkan ribozim intron dan IRES melipat
    #   secara independen. Mereka juga menemukan bahwa spacer polyA maupun
    #   polyAC tak-berstruktur sama-sama MENAIKKAN ekspresi protein, diduga
    #   karena memberi pemisahan tak-berstruktur yang lebih lega di sekitar IRES.
    #
    # Karena itu ketiga pustaka di bawah adalah varian "tak-berstruktur" dengan
    # komposisi berbeda, bukan spacer berstruktur. Peringkat akhir tetap
    # ditentukan pengukuran aksesibilitas, bukan oleh pilihan pustaka.
    "unstructured": "ACTAACTACAATCAACTACAATCAACTAACAATCAACTACA",
    # polyA: dilaporkan menaikkan ekspresi, kemungkinan lewat pengikatan PABP.
    "polyA": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    # polyAC: kontrol spacer tak-berstruktur pada studi yang sama.
    "polyAC": "ACACACACACACACACACACACACACACACACACACACACA",
}

# Panjang lengan homologi yang diuji. 9 nt ("weak") dan 19 nt ("strong") adalah
# nilai yang diuji Wesselhoeft dkk. (2018): penambahan lengan homologi menaikkan
# efisiensi splicing dari 0% menjadi 16% (9 nt) dan 48% (19 nt). Nilai 0, 30,
# dan 40 disertakan sebagai kontrol dan ekstrapolasi. [LITERATUR] + [PERKIRAAN]
HOMOLOGY_ARM_LENGTHS = (0, 9, 19, 30, 40)


def make_spacer(length: int, *, library: str = "unstructured", gc_target: float = 0.35,
                seed: int = 0) -> str:
    """Bangkitkan spacer dengan panjang dan GC yang diminta dari pustaka terpilih.

    Spacer dibangun dengan mengulang motif pustaka lalu menukar basa secara
    deterministik (seed tetap) sampai GC mendekati target, sambil menolak
    homopolimer >= 4 dan runtun GG >= 3 (menekan struktur dan pause T7).
    """
    base = SPACER_LIBRARIES[library]
    s = (base * (length // len(base) + 1))[:length]
    if library == "polyA":
        # polyA adalah spacer yang didefinisikan komposisinya; GC tidak diubah.
        return s
    rng = random.Random(seed)
    chars = list(s)
    for _ in range(600):
        cur = gc_frac("".join(chars))
        if abs(cur - gc_target) <= 0.02:
            break
        i = rng.randrange(len(chars))
        chars[i] = rng.choice("GC" if cur < gc_target else "AT")
        cand = "".join(chars)
        if _bad_spacer(cand):
            chars = list(s)  # mulai ulang dari motif dasar
    out = "".join(chars)
    return out if not _bad_spacer(out) else s


def _bad_spacer(s: str) -> bool:
    import re

    # Homopolimer >= 4 dihindari untuk sintesis dan untuk menekan pause T7,
    # KECUALI pada pustaka polyA yang memang homopolimer secara desain.
    if set(s) == {"A"}:
        return False
    if re.search(r"(.)\1{3,}", s):
        return True
    if "GGG" in s:
        return True
    return False


@dataclass
class SearchSpace:
    spacer5_len: tuple[int, int] = (10, 60)
    spacer3_len: tuple[int, int] = (10, 60)
    ha_lengths: tuple[int, ...] = HOMOLOGY_ARM_LENGTHS
    gc_target: tuple[float, float] = (0.25, 0.50)
    codon_seeds: tuple[int, ...] = (1, 2, 3, 4, 5, 6, 7, 8)
    libraries: tuple[str, ...] = tuple(SPACER_LIBRARIES)


def build_variant(
    base_parts: dict,
    *,
    protein: str,
    spacer5_len: int,
    spacer3_len: int,
    spacer_library: str,
    spacer_gc: float,
    ha_len: int,
    codon_seed: int,
    codon_w_gc: float = 1.0,
    codon_w_m6a: float = 0.5,
    codon_iterations: int = 3000,
    usage: dict | None = None,
) -> construct.Assembly:
    """Bangun satu varian konstruk dari parameter desain."""
    parts = dict(base_parts)
    parts["spacer5"] = make_spacer(spacer5_len, library=spacer_library, gc_target=spacer_gc, seed=codon_seed)
    parts["spacer3"] = make_spacer(spacer3_len, library=spacer_library, gc_target=spacer_gc, seed=codon_seed + 1)
    if "HA5" in parts and "HA3" in parts:
        seq5 = parts["HA5"].seq if hasattr(parts["HA5"], "seq") else parts["HA5"]
        seq3 = parts["HA3"].seq if hasattr(parts["HA3"], "seq") else parts["HA3"]
        parts["HA5"], parts["HA3"] = construct.trim_arms(seq5, seq3, ha_len)
    obj = codon.CodonObjective(w_gc=codon_w_gc, w_m6a=codon_w_m6a)
    opt = codon.optimize(
        protein, usage=usage, objective=obj, iterations=codon_iterations, seed=codon_seed
    )
    parts["ORF"] = opt["orf_dna"]
    return construct.build(parts)


def objectives(asm: construct.Assembly, *, temperatures=(37.0, 55.0), gc_target: float = 0.50) -> tuple:
    """Hitung keenam objektif untuk satu varian (semua diarahkan ke 'maximize')."""
    m = score.measure(asm, temperatures=temperatures)
    temps = [f"{t:g}C" for t in temperatures]
    f1 = min(m["junction"][f"circle_{t}"]["mean_unpaired_probability"] for t in temps)
    splice = []
    for t in temps:
        a = m["accessibility"][t]
        for key in ("splice5_site", "splice3_site", "intron_5half_IGS", "intron_3half_IGS"):
            if key in a and a[key]["mean"] is not None:
                splice.append(a[key]["mean"])
    f2 = min(splice) if splice else 0.0
    f3 = min((v["identity"] for v in m["ires_integrity"].values()), default=0.0)
    f4 = -max(v["longest_bp"] for v in m["dsrna"].values())
    f5 = min(
        (v["plfold"]["mean"] for v in m["start_codon"].values() if v["plfold"]["mean"] is not None),
        default=0.0,
    )
    f6 = -abs(m["gc"]["circle"] - gc_target)
    return (f1, f2, f3, float(f4), f5, f6), m


def run_study(
    base_parts: dict,
    protein: str,
    *,
    n_trials: int = 40,
    seed: int = 20260101,
    space: SearchSpace | None = None,
    temperatures=(37.0, 55.0),
    gc_target: float = 0.50,
    usage: dict | None = None,
    codon_iterations: int = 2000,
    storage: str | None = None,
    study_name: str = "circhbsag_pie",
):
    """Jalankan studi Optuna NSGA-II. Mengembalikan (study, records).

    `records` menyimpan metrik lengkap tiap trial sehingga front Pareto dapat
    diperingkat ulang dengan skor komposit Tahap 0 tanpa menjalankan ulang.
    """
    import optuna

    optuna.logging.set_verbosity(optuna.logging.WARNING)
    space = space or SearchSpace()
    records: list[dict] = []

    def objective(trial: optuna.Trial):
        p = dict(
            spacer5_len=trial.suggest_int("spacer5_len", *space.spacer5_len, step=5),
            spacer3_len=trial.suggest_int("spacer3_len", *space.spacer3_len, step=5),
            spacer_library=trial.suggest_categorical("spacer_library", list(space.libraries)),
            spacer_gc=trial.suggest_float("spacer_gc", *space.gc_target, step=0.05),
            ha_len=trial.suggest_categorical("ha_len", list(space.ha_lengths)),
            codon_seed=trial.suggest_categorical("codon_seed", list(space.codon_seeds)),
            codon_w_gc=trial.suggest_float("codon_w_gc", 0.5, 2.0, step=0.5),
            codon_w_m6a=trial.suggest_float("codon_w_m6a", 0.0, 1.5, step=0.5),
        )
        asm = build_variant(
            base_parts, protein=protein, codon_iterations=codon_iterations, usage=usage, **p
        )
        vals, m = objectives(asm, temperatures=temperatures, gc_target=gc_target)
        d = score.score(m)
        records.append(
            {
                "trial": trial.number,
                "params": p,
                "objectives": dict(
                    zip(
                        ["bsj_unpaired", "splice_access", "ires_identity", "neg_dsrna_bp",
                         "start_access", "neg_gc_dev"],
                        [round(v, 6) for v in vals],
                        strict=False,
                    )
                ),
                "composite_score": round(d.composite, 4),
                "verdict": d.verdict,
                "hard_gates_failed": [k for k, v in d.hard_gates.items() if not v["passed"]],
                "circle_length": len(asm.circle),
                "orf_dna": asm.parts["ORF"].seq,
                "spacer5": asm.parts["spacer5"].seq if "spacer5" in asm.parts else "",
                "spacer3": asm.parts["spacer3"].seq if "spacer3" in asm.parts else "",
            }
        )
        return vals

    study = optuna.create_study(
        directions=["maximize"] * 6,
        sampler=optuna.samplers.NSGAIISampler(seed=seed),
        storage=storage,
        study_name=study_name,
        load_if_exists=bool(storage),
    )
    study.optimize(objective, n_trials=n_trials, gc_after_trial=True)
    return study, records


def rank_records(records: list[dict]) -> list[dict]:
    """Peringkat varian: gerbang keras lulus dulu, lalu skor komposit menurun."""
    return sorted(
        records,
        key=lambda r: (len(r["hard_gates_failed"]) > 0, -r["composite_score"]),
    )
