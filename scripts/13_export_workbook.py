#!/usr/bin/env python3
"""Export ONE multi-sheet workbook (xlsx) with every table of a run.

Why: figures are summaries; a reviewer must be able to open the raw numbers
(every variant, every test, every primer, every motif) in Excel or Google Sheets.
Nothing is computed here: each sheet is a flattened copy of a file the other
scripts already wrote, so the workbook can never disagree with the JSON reports.
Sheets whose source file is missing are skipped and listed on the 'Ringkasan'
sheet (never filled with invented values).

Writes: <results>/tables/circhbsag_hasil_lengkap.xlsx (+ .params.json)
        <results>/tables/csv/<sheet>.csv       (same content, for Google Sheets import)
Licence of writer library: openpyxl (MIT).
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import pandas as pd
import yaml
from _common import ROOT, out_dir, parser, sidecar

from circhbsag import content


def flat(prefix: str, x, out: dict) -> dict:
    if isinstance(x, dict):
        for k, v in x.items():
            flat(f"{prefix}.{k}" if prefix else str(k), v, out)
    elif isinstance(x, list) and not (x and isinstance(x[0], (dict, list))):
        out[prefix] = ", ".join(map(str, x))
    elif isinstance(x, list):
        for i, v in enumerate(x):
            flat(f"{prefix}[{i}]", v, out)
    else:
        out[prefix] = x
    return out


def jload(p: Path):
    return json.loads(p.read_text()) if p.exists() else None


def find(res: Path, *names: str) -> Path:
    for base in (res, res / "reports", res / "tables", res / "demo"):
        for n in names:
            if (base / n).exists():
                return base / n
    return res / names[0]


def main() -> int:
    a = parser("Export full results workbook", orf=False).parse_args()
    res = Path(a.results)
    sheets: dict[str, pd.DataFrame] = {}
    missing: list[str] = []

    g0 = jload(find(res, "gate0.json"))
    if g0:
        r = g0["results"]
        dec, g = r["decision"], r["gate0"]
        sheets["Ringkasan"] = pd.DataFrame([
            ("Label data", g0.get("data_label")),
            ("Gate 0", g["decision"]), ("Ringkasan hitung", json.dumps(g["counts"])),
            ("Keputusan komposit", dec["verdict"]), ("Skor komposit", dec["composite_score"]),
            ("Seed", g0.get("seed")), ("Dibuat (UTC)", g0.get("generated_utc")),
            ("Versi circhbsag", g0.get("circhbsag_version")),
            ("Catatan", content.PROJECTION_LABEL),
        ], columns=["Item", "Nilai"])
        sheets["Gate0_uji"] = pd.DataFrame([{
            "id": t["id"], "uji": t["name"], "verdict": t["verdict"], "ambang": t["threshold"],
            "label_ambang": t["threshold_label"], "alasan_ambang": t["threshold_rationale"],
            "catatan": t.get("note", ""),
        } for t in g["tests"]])
        sheets["Gate0_metrik"] = pd.DataFrame(
            [(t["id"], k, v) for t in g["tests"] for k, v in flat("", t["metrics"], {}).items()],
            columns=["uji", "metrik", "nilai"])
        sheets["Skor_kriteria"] = pd.DataFrame([{
            "kriteria": c["name"], "bobot": c["weight"], "subskor": c["subscore"],
            "terbobot": c["weighted"], "alasan": c["rationale"]} for c in dec["criteria"]])
        sheets["Skor_kriteria_raw"] = pd.DataFrame(
            [(c["name"], k, v) for c in dec["criteria"] for k, v in flat("", c["raw"], {}).items()],
            columns=["kriteria", "metrik", "nilai"])
        sheets["Gerbang_keras"] = pd.DataFrame([{
            "gerbang": k, "lulus": v["passed"], "nilai": json.dumps(v["value"], ensure_ascii=False),
            "batas": json.dumps(v["limit"], ensure_ascii=False), "alasan": v["why"],
        } for k, v in dec["hard_gates"].items()])
        sheets["Ambang"] = pd.DataFrame(sorted(dec["all_thresholds"].items()), columns=["ambang", "nilai"])
        sheets["Koordinat"] = pd.DataFrame(
            [(k, v[0], v[1], "prekursor") for k, v in g["construct"]["precursor_coords"].items()]
            + [(k, v[0], v[1], "lingkaran") for k, v in g["construct"]["circle_coords"].items()],
            columns=["bagian", "mulai", "akhir", "cakupan"])
    else:
        missing.append("gate0.json")

    its = jload(find(res, "gate0_iterations.json"))
    if its:
        sheets["Iterasi_redesain"] = pd.DataFrame([
            {**{k: v for k, v in i.items() if k != "verdicts"},
             **{f"uji.{k}": v for k, v in i["verdicts"].items()}} for i in its])
    else:
        missing.append("gate0_iterations.json")

    opt = jload(find(res, "optimization.json"))
    if opt:
        rec = opt["results"]["optimization"]["all_records"]
        pareto = set(opt["results"]["optimization"]["pareto_trial_numbers"])
        sheets["Varian_Optuna"] = pd.DataFrame([{
            "trial": x["trial"], "pareto": x["trial"] in pareto, "skor_komposit": x["composite_score"],
            "verdict": x["verdict"], "gerbang_gagal": ", ".join(x["hard_gates_failed"]),
            "circle_nt": x["circle_length"],
            **{f"p.{k}": v for k, v in x["params"].items()},
            **{f"obj.{k}": v for k, v in x["objectives"].items()},
            "spacer5": x["spacer5"], "spacer3": x["spacer3"], "orf_dna": x["orf_dna"],
        } for x in rec]).sort_values("skor_komposit", ascending=False)
    else:
        missing.append("optimization.json")

    cod = jload(find(res, "codon.json"))
    if cod:
        c = cod["results"]["codon"]
        sheets["Kodon"] = pd.DataFrame(
            [("optimasi", k, v) for k, v in flat("", c["cost"]["metrics"], {}).items()]
            + [("naif (baseline)", k, v) for k, v in flat("", c["baseline_naive"], {}).items()]
            + [("statistik", k, v) for k, v in flat("", c["stats"], {}).items()]
            + [("asal protein", k, v) for k, v in flat("", c["protein_provenance"], {}).items()]
            + [("parameter", k, c[k]) for k in ("iterations", "seed", "length_nt")]
            + [("referensi CAI", "tabel", "data/reference/codon_usage_human.tsv (Homo sapiens)")],
            columns=["kelompok", "metrik", "nilai"])
        sheets["Kodon_ORF"] = pd.DataFrame(
            [(i // 3 + 1, c["orf_dna"][i:i + 3]) for i in range(0, len(c["orf_dna"]), 3)],
            columns=["no_kodon", "kodon"])
    else:
        missing.append("codon.json")

    pr = jload(find(res, "primers.json"))
    if pr:
        p = pr["results"]["primers"]
        ps = []
        for k in ("divergent", "bsj"):
            for x in p[k]["pairs"]:
                ps.append({"jenis": k, **x})
        sheets["Primer_pasangan"] = pd.DataFrame(ps)
        sheets["Primer_BSJ_tunggal"] = pd.DataFrame(p["bsj"]["junction_spanning_primers"])
        sheets["Probe_RNaseH"] = pd.DataFrame(p["rnase_h"]["probes"])
        sheets["Primer_IVT_reverse"] = pd.DataFrame(list(flat("", p["ivt_reverse"], {}).items()),
                                                    columns=["item", "nilai"])
    else:
        missing.append("primers.json")

    bsj = jload(find(res, "bsj_simulation.json"))
    if bsj:
        s = bsj["results"]["bsj_simulation"]
        sheets["BSJ_ringkas"] = pd.DataFrame(list(flat("", {k: v for k, v in s.items()
                                                            if k != "support"}, {}).items())
                                             + list(flat("support", {k: v for k, v in s["support"].items()
                                                                     if k != "per_read"}, {}).items()),
                                             columns=["item", "nilai"])
        sheets["BSJ_per_bacaan"] = pd.DataFrame(s["support"]["per_read"])
    else:
        missing.append("bsj_simulation.json")

    for name, fn in [("Lipatan", "fold_analysis.tsv"), ("m6A_RRACH", "m6a_rrach_sites.tsv"),
                     ("Fitur_konstruk", "construct_features.tsv")]:
        pth = find(res, fn)
        if pth.exists():
            sheets[name] = pd.read_csv(pth, sep="\t")
        else:
            missing.append(fn)
    for name, fn in [("IRES_integritas", "ires_integrity.json"), ("dsRNA", "dsrna_stems.json"),
                     ("Epitop_HLA", "epitope_hla_coverage.json"), ("Off_target", "offtarget.json")]:
        d = jload(find(res, fn))
        if d is None:
            missing.append(fn)
            continue
        sheets[name] = pd.DataFrame(list(flat("", d, {}).items()), columns=["item", "nilai"])

    fa = find(res, "construct.fa")
    if fa.exists():
        seqs, name = {}, None
        for ln in fa.read_text().splitlines():
            if ln.startswith(">"):
                name = ln[1:].strip()
                seqs[name] = ""
            elif name:
                seqs[name] += ln.strip()
        sheets["Sekuens"] = pd.DataFrame([(k, len(v), v) for k, v in seqs.items()],
                                         columns=["nama", "panjang_nt", "sekuens"])

    try:
        c = content.load(ROOT / "content", design_lock=find(res, "design_lock.json"))
        sheets["Register_asumsi"] = pd.DataFrame(content.assumptions_rows(c),
                                                 columns=["ID", "Asumsi", "Alasan", "Dampak bila salah", "Verifikasi"])
        sheets["Pertanyaan_terbuka"] = pd.DataFrame(content.open_questions_rows(c),
                                                    columns=["ID", "Pertanyaan", "Dijawab oleh", "Memblokir"])
        sheets["Software_versi"] = pd.DataFrame(content.software_rows(c),
                                                columns=["Alat", "Peran", "Versi", "Lisensi", "URL", "Status"])
    except Exception as e:  # content layer unavailable: say so, do not invent
        missing.append(f"content/ ({e.__class__.__name__})")


    # ---- bobot, justifikasi, sensitivitas, peta alur, status sumber ------------------------
    if g0:
        src = (ROOT / "src" / "circhbsag" / "score.py").read_text()
        lab = {}
        for m_ in re.finditer(r'"(\w+)":\s*[-\d.]+,\s*#\s*(\[[A-Z]+\])?\s*(.*)', src):
            lab[m_.group(1)] = (m_.group(2) or "[PERKIRAAN]", m_.group(3).strip())
        wsrc = {}
        blk = src[src.index("DEFAULT_WEIGHTS: dict"):src.index("DEFAULT_THRESHOLDS: dict")]
        cur = ""
        for ln in blk.splitlines():
            t = ln.strip()
            if t.startswith("#"):
                cur = t.lstrip("# ").strip()
            elif t.startswith('"'):
                wsrc[t.split('"')[1]] = cur
        sheets["Bobot_dan_justifikasi"] = pd.DataFrame([{
            "kriteria": c["name"], "bobot": c["weight"], "label_bobot": "[PERKIRAAN] bobot = keputusan desain pipeline, BUKAN nilai literatur",
            "alasan_mekanistik": wsrc.get(c["name"], ""), "rumus_dan_normalisasi": c["rationale"],
            "subskor_run_ini": c["subscore"], "kontribusi": c["weighted"],
        } for c in dec["criteria"]] + [{
            "kriteria": "TOTAL", "bobot": round(sum(c["weight"] for c in dec["criteria"]), 4),
            "label_bobot": "", "alasan_mekanistik": "komposit = jumlah bobot x subskor",
            "rumus_dan_normalisasi": f"GO >= {dec['thresholds']['go']}; GO-BERSYARAT >= {dec['thresholds']['conditional']}; "
                                     "satu gerbang keras gagal = NO-GO apa pun skornya",
            "subskor_run_ini": "", "kontribusi": dec["composite_score"]}])
        sheets["Ambang_berlabel"] = pd.DataFrame([
            (k, v, lab.get(k, ("", ""))[0], lab.get(k, ("", ""))[1]) for k, v in sorted(dec["all_thresholds"].items())],
            columns=["ambang", "nilai", "label", "keterangan"])
        sheets["Ambang_berlabel"]["label"] = sheets["Ambang_berlabel"]["label"].replace("", "[PERKIRAAN]")
        t0 = jload(find(res, "tahap0.json"))
        if t0 and "sensitivity" in t0["results"]:
            sv = t0["results"]["sensitivity"]
            sheets["Sensitivitas_bobot"] = pd.DataFrame(
                [(k, json.dumps(v) if isinstance(v, dict) else v) for k, v in sv.items()],
                columns=["item", "nilai"])
        sheets["Sumber_dan_status"] = pd.DataFrame([
            ("Label data run", g0.get("data_label")),
            ("Bagian sintetis (placeholder)", ", ".join(g0["provenance"].get("synthetic_parts", [])) or "tidak ada"),
            ("Boleh dikutip sebagai hasil ilmiah", g0["provenance"].get("results_may_be_reported_as_scientific")),
            ("Rubrik di berkas ini", "v1 (7 kriteria), dihitung oleh src/circhbsag/score.py"),
            ("Rubrik yang dipakai angka di TA v3/v4", "v2 (8 kriteria), lihat sheet Rubrik_v2_TA (snapshot; belum di repo)"),
            ("Kontrol positif/negatif", "tidak ada di rubrik v1 repo; kontrol PC1-3 dan NC1-5 ada pada rubrik v2 (TA)"),
        ], columns=["item", "nilai"])
    pm = yaml.safe_load((ROOT / "content" / "pipeline_map.yaml").read_text())
    sheets["Peta_alur"] = pd.DataFrame(pm["steps"])
    v2 = pm["rubric_v2_snapshot"]
    sheets["Rubrik_v2_TA"] = pd.DataFrame(
        [("SUMBER", v2["source"], "SNAPSHOT, bukan keluaran run ini"), ("GO >=", v2["go"], ""), ("GO-BERSYARAT >=", v2["conditional"], "")]
        + [(f"bobot: {k}", v, "") for k, v in v2["weights"].items()]
        + [("gerbang keras tambahan", g, "") for g in v2["hard_gates_added_vs_v1"]]
        + [("alasan perubahan dari v1", v2["why_changed"], "")] + [("kontrol", c, "") for c in v2["controls"]],
        columns=["item", "nilai", "catatan"])

    if missing:
        extra = pd.DataFrame([("Sumber tidak ditemukan (sheet dilewati)", ", ".join(missing))], columns=["Item", "Nilai"])
        sheets["Ringkasan"] = pd.concat([sheets.get("Ringkasan", pd.DataFrame(columns=["Item", "Nilai"])), extra],
                                        ignore_index=True)
    if not sheets:
        raise SystemExit(f"Tidak ada hasil di {res}: jalankan scripts/04 dan seterusnya dulu.")

    d = out_dir(a, "tables")
    out = d / "circhbsag_hasil_lengkap.xlsx"
    order = ["Ringkasan"] + [k for k in sheets if k != "Ringkasan"]
    with pd.ExcelWriter(out, engine="openpyxl") as w:
        for k in order:
            df = sheets[k]
            df.to_excel(w, sheet_name=k[:31], index=False)
            ws = w.sheets[k[:31]]
            ws.freeze_panes = "A2"
            for col in ws.columns:
                width = max(len(str(c.value)) if c.value is not None else 0 for c in col[:200])
                ws.column_dimensions[col[0].column_letter].width = min(max(10, width + 2), 60)
    cdir = d / "csv"
    cdir.mkdir(exist_ok=True)
    for k in order:
        sheets[k].to_csv(cdir / f"{k}.csv", index=False, quoting=csv.QUOTE_MINIMAL)
    print(out, f"({len(order)} sheet)")
    sidecar(out, a, {"sheets": order, "missing": missing})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
