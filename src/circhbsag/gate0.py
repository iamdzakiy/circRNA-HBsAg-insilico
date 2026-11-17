"""Gate 0: the nine-test circularisation feasibility check.

This module implements the nine tests required before any wet-lab protocol is
written. Each test returns measured numbers, a threshold, the justification for
that threshold, and a verdict in ``{PASS, REVISE, FAIL, UNAVAILABLE}``.

Two rules are deliberate:

* A threshold that comes from the literature is tagged ``LITERATUR`` and carries
  a reference key from ``content/references.yaml``. A threshold that is our own
  reasoning is tagged ``PERKIRAAN`` and says so in the report. No threshold is
  presented as established when it is not.
* A test whose tool is missing returns ``UNAVAILABLE`` together with the exact
  command that would produce the number. It never returns a guessed value, and
  ``UNAVAILABLE`` is not counted as ``PASS`` by the decision rule.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from typing import Any

import RNA

from . import fold
from .seqio import clean, gc, gc_windows, longest_homopolymer, translate

PASS, REVISE, FAIL, UNAVAILABLE = "PASS", "REVISE", "FAIL", "UNAVAILABLE"

#: Tests whose failure forbids patching the construct: the spec requires moving
#: to the backup circularisation route instead (BLOK 3, decision rule).
NO_PATCH_ON_FAIL = ("uji2_ribozyme_integrity", "uji6_ires_integrity")


@dataclass
class Test:
    """One Gate 0 test."""

    id: str
    name: str
    question: str
    verdict: str
    metrics: dict[str, Any]
    threshold: str
    threshold_label: str
    threshold_rationale: str
    refs: list[str] = field(default_factory=list)
    note: str = ""

    def to_dict(self) -> dict:
        return {
            "id": self.id,
            "name": self.name,
            "question": self.question,
            "verdict": self.verdict,
            "metrics": self.metrics,
            "threshold": self.threshold,
            "threshold_label": self.threshold_label,
            "threshold_rationale": self.threshold_rationale,
            "refs": self.refs,
            "note": self.note,
        }


# --------------------------------------------------------------------------
# Test 3 helpers: homology-arm duplex
# --------------------------------------------------------------------------


def duplex_dg(a: str, b: str, *, temperature_c: float = 37.0) -> float:
    """Free energy (kcal/mol) of the duplex between two sequences.

    Uses ``RNA.duplexfold``, the Python binding of RNAduplex: it considers only
    intermolecular pairs, which is what a homology arm pair is supposed to form.
    """
    prev = RNA.cvar.temperature
    try:
        RNA.cvar.temperature = float(temperature_c)
        d = RNA.duplexfold(_rna(a), _rna(b))
        return float(d.energy)
    finally:
        RNA.cvar.temperature = prev


def cofold_dg(a: str, b: str, *, temperature_c: float = 37.0) -> dict:
    """Co-folding energies for two sequences (RNAcofold).

    Returns the energy of the dimer and of each monomer, so that the binding
    energy can be separated from the energy each arm would have on its own.
    A duplex that looks strong only because one arm is already structured is
    not a duplex that will bring the splice sites together.
    """
    prev = RNA.cvar.temperature
    try:
        RNA.cvar.temperature = float(temperature_c)
        ab = RNA.cofold(f"{_rna(a)}&{_rna(b)}")[1]
        aa = RNA.fold(_rna(a))[1]
        bb = RNA.fold(_rna(b))[1]
        return {
            "dimer_kcal": round(float(ab), 3),
            "monomer_a_kcal": round(float(aa), 3),
            "monomer_b_kcal": round(float(bb), 3),
            "binding_kcal": round(float(ab) - float(aa) - float(bb), 3),
        }
    finally:
        RNA.cvar.temperature = prev


def arm_length_scan(
    arm5: str,
    arm3: str,
    *,
    lengths: tuple[int, ...] = (20, 25, 30, 35, 40),
    temperature_c: float = 55.0,
) -> list[dict]:
    """ΔG of the arm duplex as a function of arm length.

    Both arms are truncated from the end that faces the construct's interior,
    because that is what shortening a homology arm does in the template design:
    the 5' arm keeps its 3' end and the 3' arm keeps its 5' end.
    """
    out = []
    for n in lengths:
        if len(arm5) < n or len(arm3) < n:
            out.append({"length_nt": n, "available": False})
            continue
        a, b = arm5[-n:], arm3[:n]
        out.append(
            {
                "length_nt": n,
                "available": True,
                "duplex_kcal": round(duplex_dg(a, b, temperature_c=temperature_c), 3),
                "duplex_kcal_per_nt": round(duplex_dg(a, b, temperature_c=temperature_c) / n, 4),
                **cofold_dg(a, b, temperature_c=temperature_c),
            }
        )
    return out


# --------------------------------------------------------------------------
# Test 2 helper: Infernal / Rfam RF00028
# --------------------------------------------------------------------------


def cmsearch(sequence: str, cm_path: str | None = None, *, e_value: float = 1e-3) -> dict:
    """Run Infernal ``cmsearch`` against the group I intron covariance model.

    Returns ``{"status": "UNAVAILABLE", ...}`` with the command to run when
    either ``cmsearch`` or the covariance model is missing. The caller must not
    substitute a number for this: a missing ribozyme check is the single most
    diagnostic gap in a PIE design, and reporting it as absent is the point.
    """
    cm = cm_path or os.environ.get("RFAM_RF00028_CM", "data/reference/RF00028.cm")
    have_tool = shutil.which("cmsearch") is not None
    have_cm = os.path.exists(cm)
    command = (
        "# once, on macOS:\n"
        "mamba install -c bioconda infernal\n"
        "curl -sS https://rfam.org/family/RF00028/cm -o data/reference/RF00028.cm\n"
        "cmpress data/reference/RF00028.cm\n"
        "# then:\n"
        f"cmsearch --cut_ga --notrunc --tblout results/gate0/rf00028.tbl {cm} precursor.fa"
    )
    if not (have_tool and have_cm):
        return {
            "status": UNAVAILABLE,
            "cmsearch_found": have_tool,
            "covariance_model_found": have_cm,
            "covariance_model_path": cm,
            "command": command,
            "why": (
                "Infernal belum terpasang"
                if not have_tool
                else f"model kovarians tidak ada di {cm}"
            ),
        }
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "precursor.fa")
        tbl = os.path.join(td, "hits.tbl")
        with open(fa, "w") as fh:
            fh.write(">precursor\n" + clean(sequence) + "\n")
        proc = subprocess.run(
            ["cmsearch", "--notrunc", "-E", str(e_value), "--tblout", tbl, cm, fa],
            capture_output=True,
            text=True,
            timeout=900,
        )
        if proc.returncode != 0:
            return {"status": UNAVAILABLE, "command": command, "why": proc.stderr[-800:]}
        hits = []
        for line in open(tbl):
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 18:
                continue
            hits.append(
                {
                    "model": f[2],
                    "seq_from": int(f[7]),
                    "seq_to": int(f[8]),
                    "strand": f[9],
                    "bit_score": float(f[14]),
                    "e_value": float(f[15]),
                }
            )
    hits.sort(key=lambda h: h["e_value"])
    return {"status": "OK", "n_hits": len(hits), "hits": hits[:10], "command": command}


# --------------------------------------------------------------------------
# Test 9 helpers: sequence hygiene
# --------------------------------------------------------------------------

#: Canonical splice-site dinucleotides. A cryptic GT...AG pair inside the insert
#: is only worth reporting when both halves are present in the same order.
def cryptic_splice_sites(seq: str, *, min_intron: int = 70, max_intron: int = 2000) -> list[dict]:
    """Find canonical GT...AG pairs that could act as a cryptic intron.

    This is a motif-level screen, not a splice predictor: it reports candidate
    pairs so they can be inspected, and it deliberately does not claim that any
    of them would actually be spliced. Only pairs separated by a plausible
    intron length are reported, since an isolated GT is meaningless.
    """
    s = clean(seq)
    donors = [i for i in range(len(s) - 1) if s[i : i + 2] == "GT"]
    acceptors = {i for i in range(len(s) - 1) if s[i : i + 2] == "AG"}
    out = []
    for d in donors:
        for a in range(d + min_intron, min(len(s) - 1, d + max_intron)):
            if a in acceptors:
                out.append({"donor_pos": d + 1, "acceptor_pos": a + 1, "span_nt": a - d + 2})
                break
    return out


def polyu_runs(seq: str, *, min_run: int = 5) -> list[dict]:
    """Runs of >= min_run consecutive T (U in the RNA) — T7 pausing/termination risk."""
    s = clean(seq)
    out, i = [], 0
    while i < len(s):
        if s[i] == "T":
            j = i
            while j < len(s) and s[j] == "T":
                j += 1
            if j - i >= min_run:
                out.append({"start": i + 1, "end": j, "length": j - i})
            i = j
        else:
            i += 1
    return out


def _rna(s: str) -> str:
    return clean(s).replace("T", "U")


# --------------------------------------------------------------------------
# The nine tests
# --------------------------------------------------------------------------


def run(
    asm,
    measurements: dict,
    *,
    cm_path: str | None = None,
    arm_lengths: tuple[int, ...] = (20, 25, 30, 35, 40),
    splice_temperature_c: float = 55.0,
) -> dict:
    """Run the nine Gate 0 tests and return the full report payload.

    ``measurements`` is the output of :func:`circhbsag.score.measure`, so the
    folding work is not repeated: Gate 0 re-reads the same numbers the composite
    score was built from. That is intentional — the two views of the construct
    must never be able to disagree.
    """
    circle, precursor = asm.circle, asm.precursor
    temps = measurements.get("temperatures_c", [37.0])
    t_main = 37.0 if 37.0 in temps else float(temps[0])
    tests: list[Test] = []

    # ---------------------------------------------------------------- Test 1
    n_circle, n_prec = len(circle), len(precursor)
    size_verdict = PASS if 300 <= n_circle <= 5000 else REVISE
    tests.append(
        Test(
            id="uji1_size",
            name="Kelayakan ukuran",
            question="Apakah panjang prekursor dan lingkaran berada pada rentang yang terbukti bisa disirkularisasi?",
            verdict=size_verdict,
            metrics={
                "precursor_nt": n_prec,
                "circle_nt": n_circle,
                "bsj_position": asm.bsj_position,
                "circle_fraction_of_precursor": round(n_circle / n_prec, 4),
            },
            threshold="lingkaran 300–5000 nt",
            threshold_label="LITERATUR",
            threshold_rationale=(
                "Wesselhoeft et al. (2018) menunjukkan intron PIE yang direkayasa "
                "mensirkularisasi RNA hingga 5 kb secara efisien in vitro, sehingga 5 kb "
                "dipakai sebagai batas atas yang terbukti. Batas bawah 300 nt adalah "
                "[PERKIRAAN]: lingkaran yang jauh lebih pendek dari itu tidak cukup "
                "memuat IRES bersama ORF."
            ),
            refs=["wesselhoeft2018"],
        )
    )

    # ---------------------------------------------------------------- Test 2
    cm = cmsearch(precursor, cm_path)
    if cm["status"] == UNAVAILABLE:
        v2 = UNAVAILABLE
    else:
        best = cm["hits"][0] if cm.get("hits") else None
        v2 = PASS if best and best["e_value"] <= 1e-5 else FAIL
    tests.append(
        Test(
            id="uji2_ribozyme_integrity",
            name="Integritas ribozim intron grup I",
            question="Apakah permutasi intron masih dikenali sebagai ribozim grup I yang utuh?",
            verdict=v2,
            metrics=cm,
            threshold="E-value cmsearch <= 1e-5 terhadap Rfam RF00028, dengan P1–P10 lengkap",
            threshold_label="PERKIRAAN",
            threshold_rationale=(
                "Ambang E-value ini ketat tetapi tidak distandarkan oleh Rfam untuk kasus "
                "intron ter-permutasi; nilainya [PERLU VERIFIKASI] terhadap panduan "
                "Infernal/Rfam. Yang tidak perlu diverifikasi adalah arahnya: intron yang "
                "tidak dikenali modelnya tidak boleh dilanjutkan."
            ),
            note=(
                "Sesi cloud ini tidak bisa menjangkau Rfam (EBI ditolak proxy) dan Infernal "
                "belum terpasang, jadi uji ini WAJIB dijalankan di macOS sebelum Gate 0 "
                "dinyatakan lulus. Perintahnya ada di metrics.command."
            ),
        )
    )

    # ---------------------------------------------------------------- Test 3
    arm5 = asm.subseq("precursor", "HA5") if "HA5" in asm.precursor_coords else ""
    arm3 = asm.subseq("precursor", "HA3") if "HA3" in asm.precursor_coords else ""
    if arm5 and arm3:
        scan = arm_length_scan(
            arm5, arm3, lengths=arm_lengths, temperature_c=splice_temperature_c
        )
        full = cofold_dg(arm5, arm3, temperature_c=splice_temperature_c)
        dg = full["binding_kcal"]
        if -25.0 <= dg <= -8.0:
            v3 = PASS
        elif dg > -8.0:
            v3 = REVISE
        else:
            v3 = REVISE
        metrics3 = {
            "arm5_nt": len(arm5),
            "arm3_nt": len(arm3),
            "temperature_c": splice_temperature_c,
            "full_arms": full,
            "duplex_kcal_full": round(duplex_dg(arm5, arm3, temperature_c=splice_temperature_c), 3),
            "length_scan": scan,
        }
    else:
        v3, metrics3 = UNAVAILABLE, {"why": "konstruk ini tidak memuat HA5/HA3"}
    tests.append(
        Test(
            id="uji3_homology_arm_duplex",
            name="Pasangan lengan homologi",
            question="Apakah kedua lengan homologi berpasangan cukup kuat pada suhu splicing, tanpa menjebak RNA?",
            verdict=v3,
            metrics=metrics3,
            threshold=f"ΔG pengikatan −8 sampai −25 kcal/mol pada {splice_temperature_c:g} °C",
            threshold_label="PERKIRAAN",
            threshold_rationale=(
                "Rentang ini penalaran sendiri, bukan nilai dari literatur: di atas "
                "−8 kcal/mol pasangan tidak bertahan pada suhu reaksi, dan jauh di bawah "
                "−25 kcal/mol pasangan berisiko mengunci konformasi. Yang berdasar "
                "literatur adalah bahwa panjang lengan menentukan hasil: efisiensi "
                "sirkularisasi 0% tanpa lengan, 16% dengan 9 nt, 48% dengan 19 nt "
                "(Wesselhoeft et al. 2018). Nilai ambang ΔG-nya [PERLU VERIFIKASI]."
            ),
            refs=["wesselhoeft2018"],
        )
    )

    # ---------------------------------------------------------------- Test 4
    pf = measurements["folds"].get(f"precursor_{t_main:g}C", {})
    arm_bpp = {}
    if arm5 and arm3:
        f_prec = fold.fold(precursor, circular=False, temperature_c=t_main)
        for label, name in (("HA5", "HA5"), ("HA3", "HA3"), ("intron_5half_IGS", "intron_5half:IGS")):
            if name in asm.precursor_coords:
                a, b = asm.precursor_coords[name]
                up = fold.unpaired_probability(f_prec, list(range(a, b + 1)))
                arm_bpp[label] = {
                    "mean_unpaired": round(sum(up.values()) / len(up), 4),
                    "mean_paired": round(1 - sum(up.values()) / len(up), 4),
                }
    paired5 = arm_bpp.get("HA5", {}).get("mean_paired", 0.0)
    paired3 = arm_bpp.get("HA3", {}).get("mean_paired", 0.0)
    v4 = PASS if min(paired5, paired3) >= 0.30 else REVISE if arm_bpp else UNAVAILABLE
    tests.append(
        Test(
            id="uji4_precursor_fold",
            name="Lipatan prekursor linier",
            question="Apakah lipatan global prekursor mendekatkan kedua ujung, atau justru mengubur lengan homologi?",
            verdict=v4,
            metrics={
                "mfe_kcal": pf.get("mfe_kcal"),
                "mfe_per_nt": pf.get("mfe_per_nt"),
                "ensemble_free_energy_kcal": pf.get("efe_kcal"),
                "ensemble_diversity": pf.get("ensemble_diversity"),
                "mfe_structure_frequency": pf.get("mfe_frequency"),
                "temperature_c": t_main,
                "region_pairing": arm_bpp,
            },
            threshold="rata-rata probabilitas berpasangan di HA5 dan HA3 >= 0.30",
            threshold_label="PERKIRAAN",
            threshold_rationale=(
                "Lengan yang hampir tidak pernah berpasangan tidak bisa mendekatkan situs "
                "splice. Angka 0,30 adalah pilihan sendiri sebagai batas minimal yang "
                "masih bermakna, dan [PERLU VERIFIKASI] terhadap data sirkularisasi."
            ),
        )
    )

    # ---------------------------------------------------------------- Test 5
    cf = measurements["folds"].get(f"circle_{t_main:g}C", {})
    lin_per_nt = pf.get("mfe_per_nt")
    circ_per_nt = cf.get("mfe_per_nt")
    if lin_per_nt is not None and circ_per_nt is not None and lin_per_nt != 0:
        strain = (circ_per_nt - lin_per_nt) / abs(lin_per_nt)
        v5 = PASS if abs(strain) <= 0.25 else REVISE
    else:
        strain, v5 = None, UNAVAILABLE
    tests.append(
        Test(
            id="uji5_circle_fold",
            name="Lipatan produk sirkular",
            question="Apakah penutupan lingkaran menimbulkan tegangan struktural yang tidak masuk akal?",
            verdict=v5,
            metrics={
                "circle_mfe_kcal": cf.get("mfe_kcal"),
                "circle_mfe_per_nt": circ_per_nt,
                "precursor_mfe_per_nt": lin_per_nt,
                "relative_strain": None if strain is None else round(strain, 4),
                "circle_ensemble_diversity": cf.get("ensemble_diversity"),
                "temperature_c": t_main,
            },
            threshold="|MFE/nt lingkaran − MFE/nt prekursor| / |MFE/nt prekursor| <= 0.25",
            threshold_label="PERKIRAAN",
            threshold_rationale=(
                "Perbandingan per-nukleotida dipakai karena kedua molekul berbeda panjang. "
                "Catatan penting: model sirkular ViennaRNA membebankan energi loop eksterior "
                "pada lingkaran tak berstruktur, sehingga MFE sirkular dapat lebih TINGGI "
                "daripada linear tanpa berarti konstruknya buruk. Ambang 0,25 adalah "
                "[PERKIRAAN] untuk menandai perbedaan yang perlu dilihat manusia."
            ),
        )
    )

    # ---------------------------------------------------------------- Test 6
    ii = measurements.get("ires_integrity", {}).get(f"{t_main:g}C", {})
    ident = ii.get("identity")
    if ident is None:
        v6 = UNAVAILABLE
    elif ident >= 0.80:
        v6 = PASS
    elif ident >= 0.60:
        v6 = REVISE
    else:
        v6 = FAIL
    tests.append(
        Test(
            id="uji6_ires_integrity",
            name="Integritas struktur IRES dalam konteks sirkular",
            question="Apakah domain IRES tetap terlipat seperti saat terisolasi, setelah berada di dalam lingkaran?",
            verdict=v6,
            metrics={
                "structure_identity": ident,
                "ires_length_nt": ii.get("ires_length"),
                "ires_alone_mfe_kcal": ii.get("ires_alone_mfe_kcal"),
                "all_temperatures": measurements.get("ires_integrity", {}),
                "temperature_c": t_main,
            },
            threshold="fraksi posisi struktur IRES yang terpelihara >= 0.80 (REVISE 0.60–0.80, FAIL < 0.60)",
            threshold_label="PERKIRAAN",
            threshold_rationale=(
                "Ambang 0,80 adalah pilihan sendiri dan [PERLU VERIFIKASI]. Dasar "
                "kepentingannya dari literatur: efikasi IRES selain CVB3 berubah drastis "
                "ketika ORF di hilirnya diganti (Wesselhoeft et al. 2018), yang berarti "
                "konteks struktural memang menentukan. Perbandingan dilakukan pada struktur "
                "MFE, bukan matriks probabilitas penuh — keterbatasan yang dicatat."
            ),
            refs=["wesselhoeft2018"],
        )
    )

    # ---------------------------------------------------------------- Test 7
    ds = measurements.get("dsrna", {}).get(f"circle_{t_main:g}C", {})
    longest = ds.get("longest_bp", ds.get("length_bp"))
    if longest is None:
        v7 = UNAVAILABLE
    elif longest <= 30:
        v7 = PASS
    elif longest <= 45:
        v7 = REVISE
    else:
        v7 = FAIL
    tests.append(
        Test(
            id="uji7_dsrna_risk",
            name="Pemindaian risiko dsRNA",
            question="Apakah ada heliks berpasangan panjang yang bisa dikenali sensor imun bawaan?",
            verdict=v7,
            metrics={
                "longest_contiguous_bp": longest,
                "violations": ds.get("violations", []),
                "detail": {k: v for k, v in ds.items() if k != "violations"},
                "temperature_c": t_main,
            },
            threshold="heliks menerus terpanjang <= 30 bp (REVISE 31–45 bp, FAIL > 45 bp)",
            threshold_label="PERKIRAAN",
            threshold_rationale=(
                "Angka ~30 bp dipakai luas sebagai panjang duplex minimal yang mengaktivasi "
                "PKR, tetapi saya belum memverifikasi sumber primernya di sesi ini, jadi "
                "ambang ini [PERLU VERIFIKASI] — cari: 'minimal dsRNA length PKR activation'. "
                "Yang terverifikasi adalah relevansi kemurnian: sedikit RNA linear sudah "
                "memicu respons imun seluler yang kuat (Wesselhoeft et al. 2019)."
            ),
            refs=["wesselhoeft2019"],
        )
    )

    # ---------------------------------------------------------------- Test 8
    oc = measurements.get("orf_checks", {})
    orf_seq = asm.subseq("circle", "ORF") if "ORF" in asm.circle_coords else ""
    protein = translate(orf_seq) if orf_seq else ""
    internal_stop = "*" in protein[:-1] if protein else None
    ends_with_stop = protein.endswith("*") if protein else None
    ok8 = bool(
        orf_seq
        and oc.get("orf_starts_with_atg")
        and oc.get("orf_in_frame_length_ok")
        and ends_with_stop
        and not internal_stop
    )
    tests.append(
        Test(
            id="uji8_orf_integrity",
            name="Integritas ORF melewati sambungan balik",
            question="Apakah ORF utuh, satu kerangka, berhenti sebelum sambungan, dan tanpa stop prematur?",
            verdict=PASS if ok8 else FAIL,
            metrics={
                "orf_length_nt": len(orf_seq),
                "protein_length_aa": max(0, len(protein) - 1),
                "starts_with_atg": oc.get("orf_starts_with_atg"),
                "length_multiple_of_three": oc.get("orf_in_frame_length_ok"),
                "ends_with_stop_codon": ends_with_stop,
                "has_internal_stop": internal_stop,
                "upstream_in_frame_atg_near_start": oc.get("upstream_in_frame_atg_near"),
                "upstream_out_of_frame_atg_near_start": oc.get("upstream_out_of_frame_atg_near"),
                "longest_orf_on_doubled_circle_aa": oc.get("longest_orf_on_doubled_circle_aa"),
                "design": "ORF tunggal yang berakhir sebelum sambungan balik (bukan rolling/infinite ORF)",
            },
            threshold="semua pemeriksaan boolean benar; tanpa stop internal",
            threshold_label="PERKIRAAN",
            threshold_rationale=(
                "Ini bukan ambang numerik melainkan syarat biner yang diverifikasi dengan "
                "kode, bukan dengan membaca sekuens. Desainnya sengaja ORF tunggal: "
                "produk proteinnya satu spesies dengan panjang tertentu, sedangkan rolling "
                "ORF menghasilkan konkatamer yang tidak diinginkan untuk antigen vaksin."
            ),
        )
    )

    # ---------------------------------------------------------------- Test 9
    insert = circle
    typeIIS = measurements.get("motifs", {}).get("precursor", {}).get("restriction", {})
    t7term = measurements.get("motifs", {}).get("precursor", {}).get("t7_terminator", [])
    # gc_windows yields (start, gc) tuples; keep the ones outside the band.
    windows = gc_windows(insert, window=50)
    bad_windows = [
        {"start": int(start), "gc": round(float(value), 4)}
        for start, value in windows
        if value < 0.35 or value > 0.70
    ]
    gc_all = gc(insert)
    pu = polyu_runs(insert, min_run=5)
    cryptic = cryptic_splice_sites(insert)
    hp = longest_homopolymer(insert)
    hygiene_fail = bool(
        any(typeIIS.values() if isinstance(typeIIS, dict) else []) or t7term
    )
    hygiene_revise = bool(
        bad_windows or pu or not (0.45 <= gc_all <= 0.60) or (hp or {}).get("length", 0) >= 7
    )
    v9 = FAIL if hygiene_fail else (REVISE if hygiene_revise else PASS)
    tests.append(
        Test(
            id="uji9_sequence_hygiene",
            name="Higienitas sekuens",
            question="Apakah sekuens bebas situs kloning internal, motif terminasi T7, poly-U, dan GC ekstrem?",
            verdict=v9,
            metrics={
                "type_IIS_sites": typeIIS,
                "t7_terminator_motifs": t7term,
                "gc_global": round(gc_all, 4),
                "gc_window_50nt_out_of_range": bad_windows[:20],
                "n_gc_windows_out_of_range": len(bad_windows),
                "polyu_runs_ge5": pu,
                "longest_homopolymer": hp,
                "cryptic_gt_ag_pairs": cryptic[:20],
                "n_cryptic_gt_ag_pairs": len(cryptic),
            },
            threshold=(
                "FAIL bila ada situs Tipe IIS internal atau motif terminator T7; "
                "REVISE bila GC global di luar 45–60%, ada jendela 50 nt di luar 35–70%, "
                "poly-U >= 5, atau homopolimer >= 7"
            ),
            threshold_label="PERKIRAAN",
            threshold_rationale=(
                "Rentang GC 45–60% dan larangan poly-U >= 5 serta poly-A >= 7 diambil dari "
                "spesifikasi proyek ini sendiri (BLOK 4.2), bukan dari satu sumber primer, "
                "jadi [PERLU VERIFIKASI]. Pasangan GT–AG dilaporkan sebagai kandidat untuk "
                "diperiksa manusia, BUKAN sebagai prediksi splicing: tanpa prediktor splice "
                "site sungguhan, keberadaan motif saja tidak cukup untuk menyimpulkan apa pun."
            ),
        )
    )

    return _assemble(tests, asm, temps)


def _assemble(tests: list[Test], asm, temps) -> dict:
    counts = {v: sum(1 for t in tests if t.verdict == v) for v in (PASS, REVISE, FAIL, UNAVAILABLE)}
    failed_no_patch = [t.id for t in tests if t.verdict == FAIL and t.id in NO_PATCH_ON_FAIL]
    if failed_no_patch:
        decision = "FAIL-PINDAH-JALUR"
        statement = (
            "Gate 0 TIDAK lulus pada uji yang tidak boleh ditambal "
            f"({', '.join(failed_no_patch)}); desain harus berpindah ke jalur sirkularisasi "
            "cadangan dan Gate 0 dijalankan ulang dari awal untuk arsitektur itu."
        )
    elif counts[FAIL]:
        decision = "FAIL"
        statement = "Gate 0 TIDAK lulus; konstruk harus dirancang ulang sebelum protokol basah ditulis."
    elif counts[UNAVAILABLE]:
        decision = "BELUM-LENGKAP"
        statement = (
            f"Gate 0 belum bisa dinyatakan lulus: {counts[UNAVAILABLE]} uji belum bisa dijalankan "
            "di lingkungan ini dan harus dijalankan di macOS terlebih dahulu."
        )
    elif counts[REVISE]:
        decision = "REVISE"
        statement = (
            f"Gate 0 memerlukan siklus perbaikan desain: {counts[REVISE]} uji berstatus REVISE. "
            "Ubah satu parameter per iterasi dan jalankan ulang uji yang terdampak."
        )
    else:
        decision = "PASS"
        statement = "Gate 0 lulus pada kesembilan uji; konstruk boleh dibawa ke tahap berikutnya."
    return {
        "gate": "Gate 0",
        "decision": decision,
        "statement": statement,
        "counts": counts,
        "tests": [t.to_dict() for t in tests],
        "construct": {
            "precursor_nt": len(asm.precursor),
            "circle_nt": len(asm.circle),
            "bsj_position": asm.bsj_position,
            "precursor_coords": {k: list(v) for k, v in asm.precursor_coords.items()},
            "circle_coords": {k: list(v) for k, v in asm.circle_coords.items()},
        },
        "temperatures_c": list(temps),
    }


def design_lock(gate: dict, measurements: dict, decision_payload: dict, *, provenance: dict) -> dict:
    """The one set of numbers every other deliverable must quote.

    Written to ``design_lock.json`` so the thesis document, the slide deck and
    the website read the same values instead of copying them by hand. Any
    disagreement between deliverables then becomes a stale file, which CI can
    detect, rather than a number nobody can trace.
    """
    synthetic = bool(provenance.get("synthetic_parts"))
    t = gate["temperatures_c"][0]
    by_id = {x["id"]: x for x in gate["tests"]}
    arm = by_id.get("uji3_homology_arm_duplex", {}).get("metrics", {})
    return {
        "schema": "circhbsag/design-lock/1",
        "data_status": "DEMO-PLACEHOLDER" if synthetic else "SEKUENS-PENGGUNA",
        "data_label": (
            "[DEMO — SEKUENS PLACEHOLDER] angka di bawah mengilustrasikan alur kerja, "
            "bukan sifat konstruk sebenarnya, dan tidak boleh dikutip sebagai hasil"
            if synthetic
            else "[HASIL IN SILICO] dihitung dari sekuens yang disediakan pengguna"
        ),
        "projection_label": "Proyeksi berbasis literatur — bukan data eksperimen",
        "ires": {"name": "CVB3", "length_nt": by_id.get("uji6_ires_integrity", {}).get("metrics", {}).get("ires_length_nt")},
        "circularisation_method": "PIE (permuted intron–exon), intron grup I Anabaena pre-tRNA-Leu",
        "lengths_nt": {
            "precursor": gate["construct"]["precursor_nt"],
            "circle": gate["construct"]["circle_nt"],
            "bsj_position": gate["construct"]["bsj_position"],
            "orf": by_id.get("uji8_orf_integrity", {}).get("metrics", {}).get("orf_length_nt"),
            "protein_aa": by_id.get("uji8_orf_integrity", {}).get("metrics", {}).get("protein_length_aa"),
        },
        "homology_arms": {
            "arm5_nt": arm.get("arm5_nt"),
            "arm3_nt": arm.get("arm3_nt"),
            "binding_dg_kcal": (arm.get("full_arms") or {}).get("binding_kcal"),
            "duplex_dg_kcal": arm.get("duplex_kcal_full"),
            "temperature_c": arm.get("temperature_c"),
            "length_scan": arm.get("length_scan"),
        },
        "folding": {
            "temperature_c": t,
            "circle_mfe_kcal": measurements["folds"].get(f"circle_{t:g}C", {}).get("mfe_kcal"),
            "circle_mfe_per_nt": measurements["folds"].get(f"circle_{t:g}C", {}).get("mfe_per_nt"),
            "circle_ensemble_diversity": measurements["folds"].get(f"circle_{t:g}C", {}).get("ensemble_diversity"),
            "precursor_mfe_kcal": measurements["folds"].get(f"precursor_{t:g}C", {}).get("mfe_kcal"),
            "precursor_ensemble_diversity": measurements["folds"].get(f"precursor_{t:g}C", {}).get("ensemble_diversity"),
        },
        "gate0": {
            "decision": gate["decision"],
            "statement": gate["statement"],
            "counts": gate["counts"],
            "verdicts": {x["id"]: x["verdict"] for x in gate["tests"]},
        },
        "composite_score": decision_payload.get("composite_score"),
        "composite_verdict": decision_payload.get("verdict"),
        "sequence_provenance": {
            "synthetic_parts": provenance.get("synthetic_parts", []),
            "accessions": {
                k: v.get("accession")
                for k, v in (provenance.get("parts") or {}).items()
                if isinstance(v, dict)
            },
        },
        "ribozyme_check": by_id.get("uji2_ribozyme_integrity", {}).get("metrics", {}).get("status"),
    }


def markdown(gate: dict) -> str:
    """Render the Gate 0 report as Markdown."""
    sym = {PASS: "PASS", REVISE: "REVISE", FAIL: "FAIL", UNAVAILABLE: "BELUM BISA DIJALANKAN"}
    out = [
        "# Laporan Gate 0 — kelayakan sirkularisasi",
        "",
        f"**Keputusan: {gate['decision']}** — {gate['statement']}",
        "",
        f"Ringkasan: {gate['counts'][PASS]} PASS, {gate['counts'][REVISE]} REVISE, "
        f"{gate['counts'][FAIL]} FAIL, {gate['counts'][UNAVAILABLE]} belum bisa dijalankan.",
        "",
        f"Prekursor {gate['construct']['precursor_nt']} nt; lingkaran "
        f"{gate['construct']['circle_nt']} nt; sambungan balik pada posisi "
        f"{gate['construct']['bsj_position']}.",
        "",
        "| # | Uji | Verdict | Ambang | Dasar ambang |",
        "|---|---|---|---|---|",
    ]
    for i, t in enumerate(gate["tests"], 1):
        out.append(
            f"| {i} | {t['name']} | **{sym[t['verdict']]}** | {t['threshold']} | {t['threshold_label']} |"
        )
    out += ["", "## Rincian tiap uji", ""]
    for i, t in enumerate(gate["tests"], 1):
        out += [
            f"### Uji {i} — {t['name']}",
            "",
            f"*Pertanyaan:* {t['question']}",
            "",
            f"*Verdict:* **{sym[t['verdict']]}**",
            "",
            f"*Ambang* ({t['threshold_label']}): {t['threshold']}",
            "",
            f"*Dasar ambang:* {t['threshold_rationale']}",
            "",
        ]
        if t["note"]:
            out += [f"*Catatan:* {t['note']}", ""]
        out += [
            "```json",
            json.dumps(t["metrics"], indent=2, ensure_ascii=False)[:4000],
            "```",
            "",
        ]
    return "\n".join(out)
