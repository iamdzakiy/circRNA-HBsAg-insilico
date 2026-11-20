"""Rubrik skor v2 (revisi setelah temuan bahwa kontrol negatif v1 skor tertinggi).

Perubahan terhadap v1 (semuanya dapat ditelusuri ke data atau alasan mekanistik,
bukan ke penurunan ambang; ambang GO 0,70 dan GO-BERSYARAT 0,55 tidak diubah):

 1. KRITERIA BARU `splicing_competence`: prediksi efisiensi sirkularisasi dari
    data sumber Shen dkk. 2025 (Fig. 2b-g dan 6d). v1 tidak punya kriteria ini,
    sehingga panjang ekson, lengan, dan spacer tidak berpengaruh pada skor.
 2. `ires_integrity` diukur pada 37 C (suhu translasi). Pada 55 C IRES hanya ada
    sekilas selama reaksi sirkularisasi, lalu dimurnikan dan dilipat ulang; 55 C
    tetap dilaporkan sebagai informasi.
 3. `junction_accessibility`: aksesibilitas situs splice dinilai RELATIF terhadap
    konstruk minimal intron+ekson (pelindung intrinsik intron), bukan terhadap
    nilai mutlak. Pada v1 komponen ini ~0 untuk semua konstruk (IGS secara
    alami berpasangan di dalam struktur intron), sehingga tidak informatif.
 4. `circle_stability` dinormalkan terhadap nol acak: 25 pengacakan dinukleotida
    lingkaran yang sama (z-skor MFE/nt dan keragaman ensemble/nt), bukan ambang
    mutlak (v1: keragaman <= 0,12/nt tidak pernah tercapai oleh lingkaran ber-IRES).
 5. GERBANG KERAS baru: IRES identik 100% dengan AF169665.1, intron identik 100%
    dengan intron Ana acuan, ORF menerjemahkan protein S persis tanpa stop
    internal, lengan homologi <= 50 nt (Fig. 2f: >=75 nt menurunkan ekspresi
    >90%), dan efisiensi dasar prediksi ekson >= 25%.
"""
from __future__ import annotations

import json
import random

from . import access, fold

WEIGHTS2 = {
    "splicing_competence": 0.25,
    "junction_accessibility": 0.15,
    "ires_integrity": 0.15,
    "dsrna_burden": 0.10,
    "start_codon_accessibility": 0.10,
    "circle_stability": 0.10,
    "sequence_cleanliness": 0.10,
    "motif_burden": 0.05,
}
GO, COND = 0.70, 0.55

_MODEL = None
def model(path="data/real/splice_model_data.json"):
    global _MODEL
    if _MODEL is None:
        _MODEL = json.load(open(path))
    return _MODEL

def _interp(table: dict, x: float) -> float:
    pts = sorted((int(k), v) for k, v in table.items())
    if x <= pts[0][0]: return pts[0][1]
    if x >= pts[-1][0]: return pts[-1][1]
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if x0 <= x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)

def _ramp(x, lo, hi):
    return max(0.0, min(1.0, (x - lo) / (hi - lo))) if hi != lo else float(x >= hi)

def predict_splicing(e5_seq: str, e3_seq: str, spacer_nt: int = 0, arm_nt: int = 0) -> dict:
    """Prediksi efisiensi sirkularisasi 55 C (0-1) dari data Shen dkk. 2025.

    Dasar: e_dasar = f5(a) * f3(b) / e(16,16), dengan f5 dari Fig. 2c (ekson 5'
    divariasikan, 3' = 16 nt) dan f3 dari Fig. 2d (kebalikannya). Untuk ekson
    5' 3 nt, identitas basa N-3 N-2 U dikoreksi dengan Fig. 6d relatif terhadap CUU.
    Peningkat (lengan: Fig. 2e; spacer dalam: Fig. 2g, nilai yang konservatif)
    menurunkan fraksi sisa belum-tersambung; hanya yang terkuat dipakai agar tidak
    dihitung ganda. SEMUA ini PERKIRAAN lintas-konteks (muatan EGFP pada paper).
    """
    M = model()
    a, b = len(e5_seq), len(e3_seq)
    ref = M["fig2b"]["16"] / 100
    f5 = _interp(M["fig2c_5exon"], a) / 100
    f3 = _interp(M["fig2d_3exon"], b) / 100
    e_base = f5 * f3 / ref
    mod = 1.0
    if a == 3:
        s = e5_seq.upper().replace("T", "U")
        d = M["fig6d_5exon"]
        if s in d: mod = d[s] / d["CUU"]
    e_base = min(0.98, e_base * mod)
    # peningkat
    base_arm = M["fig2e_arm"]["0"] / 100
    base_sp = M["fig2g_spacer"]["No spacer"] / 100
    r_arm = 1.0
    if arm_nt > 0:
        e_arm = _interp(M["fig2e_arm"], arm_nt) / 100
        r_arm = max(0.0, (1 - e_arm) / (1 - base_arm))
    r_sp = 1.0
    if spacer_nt > 0:
        inner = {int(k.split(",")[1].split()[0]): v for k, v in M["fig2g_spacer"].items() if k.startswith("Inner")}
        e_sp = _interp(inner, spacer_nt) / 100
        r_sp = max(0.0, (1 - e_sp) / (1 - base_sp))
    R = min(r_arm, r_sp)
    e_tot = 1 - (1 - e_base) * R
    return {"e_base": round(e_base, 4), "f5": round(f5, 4), "f3": round(f3, 4), "exon5_identity_modifier": round(mod, 3),
            "residual_factor_arm": round(r_arm, 3), "residual_factor_spacer": round(r_sp, 3),
            "e_total": round(min(0.99, e_tot), 4)}

def arm_expression_factor(arm_nt: int) -> float:
    M = model()["fig2f_arm_MFI"]
    return _interp(M, arm_nt) / max(v for v in M.values())

# ---- pengacakan dinukleotida (Altschul & Erikson 1985; implementasi Jiang dkk. 2008) ----
def dinuc_shuffle(s: str, rng: random.Random) -> str:
    from collections import defaultdict
    last = s[-1]
    edges = defaultdict(list)
    for a, b in zip(s, s[1:]): edges[a].append(b)
    nodes = list(edges)
    while True:
        lastedge = {}
        for x in nodes:
            if x == last: continue
            lastedge[x] = rng.choice(edges[x])
        # cek pohon berakar di last
        ok = True
        for x in nodes:
            if x == last: continue
            seen, y = set(), x
            while y != last and y not in seen:
                seen.add(y); y = lastedge.get(y, last) if y in lastedge else last
            if y != last: ok = False; break
        if ok: break
    out = {}
    for x in nodes:
        rest = list(edges[x])
        if x != last:
            rest.remove(lastedge[x])
            rng.shuffle(rest); rest.append(lastedge[x])
        else:
            rng.shuffle(rest)
        out[x] = rest
    cur, res = s[0], [s[0]]
    for _ in range(len(s) - 1):
        nxt = out[cur].pop(0); res.append(nxt); cur = nxt
    return "".join(res)

def shuffle_null(circle: str, n=25, seed=11, temperature=37.0) -> dict:
    rng = random.Random(seed); mf, dv = [], []
    for _ in range(n):
        f = fold.fold(dinuc_shuffle(circle, rng), circular=True, temperature_c=temperature)
        mf.append(f.mfe_per_nt); dv.append(f.ensemble_diversity / f.sequence_length)
    mu = lambda v: sum(v) / len(v)
    sd = lambda v: (sum((x - mu(v)) ** 2 for x in v) / (len(v) - 1)) ** 0.5
    return {"n": n, "temperature_c": temperature, "mfe_per_nt_mean": mu(mf), "mfe_per_nt_sd": sd(mf),
            "div_per_nt_mean": mu(dv), "div_per_nt_sd": sd(dv)}

def minimal_splice_baseline(intron: str, e2: str, e1: str, igs: tuple[int, int], temps=(55.0,)) -> dict:
    """Aksesibilitas situs splice pada konstruk minimal T7-GGG? + intron + E2 + E1 (pelindung intrinsik)."""
    pre = intron + e2 + e1
    regs = {}
    n_int = len(intron)
    regs["splice5_site"] = (max(1, len(pre) - 9), len(pre))
    regs["splice3_site"] = (max(1, n_int + 1 - 10), min(len(pre), n_int + 1 + 9))
    regs["intron_3half_IGS"] = (igs[0] + 1, igs[1])
    out = {}
    for t in temps:
        out[f"{t:g}C"] = access.accessibility(pre, regs, ulength=10, window=80, temperature_c=t)
    return out

def score2(m: dict, ctx: dict) -> dict:
    """Skor v2. ctx: e5_seq,e3_seq,spacer_nt,arm_nt, null, baseline, ires_ok, intron_ok, orf_ok(dict), cloning."""
    th = {"dsrna_max": 30, "dsrna_hard": 45, "homopolymer_max": 6}
    temps = m["temperatures_c"]; w = WEIGHTS2
    gates = {}
    pre = m["motifs"]["precursor"]
    n_t2 = pre["restriction"].get("BsaI", 0) + pre["restriction"].get("BsmBI_Esp3I", 0)
    gates["no_internal_typeIIS"] = {"passed": n_t2 == 0 or ctx.get("cloning") == "gibson", "value": n_t2}
    gates["no_premature_T7_terminator"] = {"passed": len(pre["t7_terminator"]) == 0, "value": len(pre["t7_terminator"])}
    oc = m.get("orf_checks", {})
    gates["orf_intact_and_unique_start"] = {"passed": bool(oc.get("orf_starts_with_atg") and oc.get("orf_in_frame_length_ok") and not oc.get("upstream_in_frame_atg_near")), "value": None}
    ds_max = max(m["dsrna"][k]["longest_bp"] for k in m["dsrna"])
    gates["dsrna_below_hard_limit"] = {"passed": ds_max < th["dsrna_hard"], "value": ds_max}
    gates["circle_length_within_limit"] = {"passed": m["assembly"]["circle_length"] <= 3000, "value": m["assembly"]["circle_length"]}
    gates["ires_sequence_identical_AF169665"] = {"passed": bool(ctx["ires_ok"]), "value": ctx["ires_ok"]}
    gates["intron_identical_Ana_reference"] = {"passed": bool(ctx["intron_ok"]), "value": ctx["intron_ok"]}
    gates["orf_translates_S_protein_no_internal_stop"] = {"passed": bool(ctx["orf_ok"]), "value": ctx["orf_ok"]}
    gates["homology_arm_le_50nt"] = {"passed": ctx["arm_nt"] <= 50, "value": ctx["arm_nt"]}
    if ctx.get("splice_override") is not None:
        x = ctx["splice_override"]
        sp = {"e_base": x, "e_total": x, "source": "Fig. 6f Source Data (situs split IRES, tanpa bekas)"}
    else:
        sp = predict_splicing(ctx["e5_seq"], ctx["e3_seq"], ctx["spacer_nt"], ctx["arm_nt"])
    gates["exon_base_efficiency_ge_25pct"] = {"passed": sp["e_base"] >= 0.25, "value": sp["e_base"]}
    crit = {}
    # 1 splicing
    crit["splicing_competence"] = (_ramp(sp["e_total"], 0.20, 0.80), sp)
    # 2 junction
    bsj = [m["junction"][f"circle_{t:g}C"]["mean_unpaired_probability"] for t in temps]
    bsj_sub = _ramp(min(bsj), 0.05, 0.50)
    ratios = {}
    for t in temps:
        key = f"{t:g}C"
        if key in ctx["baseline"]:
            for r in ("splice5_site", "splice3_site", "intron_3half_IGS"):
                a = m["accessibility"][key].get(r if r != "intron_3half_IGS" else "intron_3half_IGS", {}).get("mean")
                b = ctx["baseline"][key].get(r, {}).get("mean")
                if a is not None and b:
                    ratios[f"{r}@{key}"] = round(a / b, 3) if b > 1e-4 else None
    rv = [v for v in ratios.values() if v is not None]
    ctx_sub = _ramp(min(rv), 0.25, 0.90) if rv else 0.5
    crit["junction_accessibility"] = (0.6 * bsj_sub + 0.4 * ctx_sub, {"bsj_worst": round(min(bsj), 4), "splice_site_ratio_vs_minimal": ratios})
    # 3 IRES 37C
    ii = m["ires_integrity"]
    id37 = ii["37C"]["identity"] if "37C" in ii else None
    crit["ires_integrity"] = (_ramp(id37, 0.70, 0.90) if id37 is not None else 0.0, {"identity_37C": id37, "identity_55C_informational": ii.get("55C", {}).get("identity")})
    # 4 dsRNA
    crit["dsrna_burden"] = (1 - _ramp(ds_max, th["dsrna_max"] * 0.6, th["dsrna_hard"]), {"longest_bp": ds_max})
    # 5 start codon
    sc = m["start_codon"]
    sv = [v["plfold"]["mean"] for v in sc.values() if v["plfold"]["mean"] is not None]
    gv = [v["mean_unpaired_global_fold"] for v in sc.values()]
    crit["start_codon_accessibility"] = (0.5 * _ramp(min(sv), 0.01, 0.35) + 0.5 * _ramp(min(gv), 0.05, 0.60), {"plfold": {k: v["plfold"]["mean"] for k, v in sc.items()}, "global": {k: v["mean_unpaired_global_fold"] for k, v in sc.items()}})
    # 6 stability vs null (37 C)
    f37 = m["folds"]["circle_37C"]; N = ctx["null"]
    mfe = f37["mfe_per_nt"]; div = f37["ensemble_diversity"] / m["assembly"]["circle_length"]
    z_m = (mfe - N["mfe_per_nt_mean"]) / N["mfe_per_nt_sd"]
    z_d = (div - N["div_per_nt_mean"]) / N["div_per_nt_sd"]
    crit["circle_stability"] = (0.5 * _ramp(-z_m, 0, 6) + 0.5 * _ramp(-z_d, 0, 6), {"z_mfe_per_nt": round(z_m, 2), "z_diversity_per_nt": round(z_d, 2)})
    # 7 cleanliness
    gcv = len(pre["gc_violations"]); hp = m["homopolymer"]["precursor"][1]; gcc = m["gc"]["circle"]
    band = 1.0 if 0.35 <= gcc <= 0.70 else max(0.0, 1 - (0.35 - gcc if gcc < 0.35 else gcc - 0.70) / 0.10)
    crit["sequence_cleanliness"] = (0.4 * (1.0 if n_t2 == 0 else 0.0) + 0.2 * (1 - _ramp(gcv, 0, 10)) + 0.2 * (1 - _ramp(hp, 6, 10)) + 0.2 * band, {"typeIIS": n_t2, "gc_window_violations": gcv, "homopolymer": hp, "gc": gcc})
    # 8 motifs
    m6a = m["motifs"]["circle"]["m6a_density_per_100nt"]; nav = sum(m["motifs"]["circle"]["avoid_motifs"].values())
    crit["motif_burden"] = (0.6 * (1 - _ramp(m6a, 1.0, 4.0)) + 0.4 * (1 - _ramp(nav / max(1.0, m["assembly"]["circle_length"] / 1000), 0, 10)), {"m6a_per_100nt": m6a})
    comp = sum(w[k] * crit[k][0] for k in w) / sum(w.values())
    failed = [k for k, v in gates.items() if not v["passed"]]
    verdict = "NO-GO" if failed else ("GO" if comp >= GO else ("GO-BERSYARAT" if comp >= COND else "NO-GO"))
    return {"composite": round(comp, 4), "verdict": verdict, "gates_failed": failed,
            "criteria": {k: {"weight": w[k], "subscore": round(v[0], 4), "raw": v[1]} for k, v in crit.items()},
            "gates": gates}
