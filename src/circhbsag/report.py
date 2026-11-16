"""Penulisan laporan Tahap 0/Tahap 1 dalam JSON, Markdown, dan HTML mandiri.

Setiap laporan memuat blok PROVENANS di bagian paling atas: status tiap bagian
sekuens (disediakan pengguna atau sintetis-uji), versi alat, seed, dan apakah
hasil boleh dikutip secara ilmiah. Ini mencegah angka uji perangkat lunak
terbawa ke dokumen tugas akhir tanpa label.
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

from . import __version__, licenses


def _json_default(o):
    if isinstance(o, tuple):
        return list(o)
    return str(o)


def envelope(payload: dict, *, provenance: dict | None = None, seed: int | None = None) -> dict:
    """Bungkus hasil dengan metadata reprodusibilitas."""
    prov = provenance or {}
    synthetic = prov.get("synthetic_parts") or []
    return {
        "report_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "circhbsag_version": __version__,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "seed": seed,
        "provenance": prov,
        "data_label": (
            "[HASIL IN SILICO] — sekuens masukan disediakan pengguna"
            if not synthetic
            else "[UJI PERANGKAT LUNAK] — memuat sekuens SINTETIS; "
            "angka di bawah TIDAK boleh dikutip sebagai hasil ilmiah"
        ),
        "tool_licenses": licenses(),
        "results": payload,
    }


def write_json(data: dict, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, default=_json_default, ensure_ascii=False))
    return path


# --------------------------------------------------------------------------
# Markdown
# --------------------------------------------------------------------------


def markdown(report: dict) -> str:
    """Laporan Markdown ringkas namun lengkap untuk dilampirkan ke dokumen TA."""
    r = report["results"]
    dec = r.get("decision", {})
    asm = r.get("measurements", {}).get("assembly", {})
    L: list[str] = []
    a = L.append

    a("# Laporan Tahap 0 — Kelayakan Sirkularisasi circRNA-HBsAg (PIE)")
    a("")
    a(f"**Label data:** {report['data_label']}")
    a("")
    a(f"- Dihasilkan (UTC): `{report['generated_utc']}`")
    a(f"- circhbsag `{report['circhbsag_version']}`, Python `{report['python']}`")
    a(f"- Seed: `{report.get('seed')}`")
    a(f"- ViennaRNA: `{report['tool_licenses']['ViennaRNA']['version']}` "
      f"({report['tool_licenses']['ViennaRNA']['license']})")
    a("")

    prov = report.get("provenance", {}).get("parts", {})
    if prov:
        a("## Provenans sekuens masukan")
        a("")
        a("| Bagian | Panjang (nt) | Status | Aksesi/sumber |")
        a("|---|---:|---|---|")
        for k, v in prov.items():
            src = v.get("accession") or v.get("source") or v.get("file") or "inline"
            a(f"| `{k}` | {v.get('length','?')} | {v.get('provenance_status','?')} | {src} |")
        a("")

    if dec:
        a("## Keputusan")
        a("")
        a(f"### {dec['verdict']}  (skor komposit {dec['composite_score']:.3f})")
        a("")
        a(f"Ambang: GO >= {dec['thresholds']['go']}, "
          f"GO-BERSYARAT >= {dec['thresholds']['conditional']}, selain itu NO-GO.")
        a("")
        a("#### Gerbang keras")
        a("")
        a("| Gerbang | Lulus | Nilai | Batas | Mengapa |")
        a("|---|:--:|---|---|---|")
        for k, v in dec["hard_gates"].items():
            a(f"| `{k}` | {'YA' if v['passed'] else '**TIDAK**'} | `{v['value']}` | `{v['limit']}` | {v['why']} |")
        a("")
        a("#### Kriteria berbobot")
        a("")
        a("| Kriteria | Bobot | Subskor | Kontribusi |")
        a("|---|---:|---:|---:|")
        for c in dec["criteria"]:
            a(f"| `{c['name']}` | {c['weight']:.2f} | {c['subscore']:.3f} | {c['weighted']:.3f} |")
        a("")
        a("Rincian penalaran tiap kriteria:")
        a("")
        for c in dec["criteria"]:
            a(f"- **{c['name']}** (subskor {c['subscore']:.3f}): {c['rationale']}")
        a("")

    m = r.get("measurements", {})
    if m:
        a("## Pengukuran kunci")
        a("")
        a(f"- Panjang prekursor: **{asm.get('precursor_length')} nt**; "
          f"lingkaran: **{asm.get('circle_length')} nt**; BSJ pada posisi {asm.get('bsj_position')}")
        a("")
        a("### Lipatan (RNAfold, `--circ` untuk lingkaran)")
        a("")
        a("| Molekul | Suhu | MFE (kcal/mol) | MFE/nt | Energi ensemble | Keragaman ensemble | Freq. MFE |")
        a("|---|---:|---:|---:|---:|---:|---:|")
        for k, v in m.get("folds", {}).items():
            a(f"| `{k}` | {v['temperature_c']:g} C | {v['mfe_kcal']:.2f} | {v['mfe_per_nt']:.4f} | "
              f"{v['ensemble_free_energy_kcal']:.2f} | {v['ensemble_diversity']:.2f} | "
              f"{v['frequency_mfe_in_ensemble']:.4g} |")
        a("")
        if m.get("junction"):
            a("### Sambungan balik (BSJ)")
            a("")
            a("| Lipatan | P(tak-berpasangan) rerata | P(berpasangan) rerata | Struktur jendela |")
            a("|---|---:|---:|---|")
            for k, v in m["junction"].items():
                a(f"| `{k}` | {v['mean_unpaired_probability']:.4f} | "
                  f"{v['mean_paired_probability']:.4f} | `{v['structure_window']}` |")
            a("")
        if m.get("accessibility"):
            a("### Aksesibilitas RNAplfold pada prekursor")
            a("")
            a("| Suhu | Wilayah | Rerata | Minimum |")
            a("|---|---|---:|---:|")
            for t, regions in m["accessibility"].items():
                for name, v in regions.items():
                    if v.get("mean") is None:
                        continue
                    a(f"| {t} | `{name}` | {v['mean']:.5f} | {v['min']:.5f} |")
            a("")
        if m.get("ires_integrity"):
            a("### Integritas struktur IRES")
            a("")
            a("| Suhu | Identitas status lipatan | MFE IRES sendiri (kcal/mol) |")
            a("|---|---:|---:|")
            for t, v in m["ires_integrity"].items():
                a(f"| {t} | {v['identity']:.3f} | {v['ires_alone_mfe_kcal']:.2f} |")
            a("")
        if m.get("dsrna"):
            a("### Beban dsRNA (heliks beruntun terpanjang)")
            a("")
            a("| Lipatan | Terpanjang (bp) | Jumlah heliks >= ambang |")
            a("|---|---:|---:|")
            for k, v in m["dsrna"].items():
                a(f"| `{k}` | {v['longest_bp']} | {len(v.get('violations', []))} |")
            a("")
        if m.get("motifs"):
            a("### Pemindaian motif")
            a("")
            for which in ("precursor", "circle"):
                mm = m["motifs"].get(which, {})
                if not mm:
                    continue
                a(f"**{which}** (panjang {mm['length']} nt): "
                  f"situs Tipe IIS = {mm['restriction']}, "
                  f"terminator T7 = {len(mm['t7_terminator'])}, "
                  f"m6A RRACH = {mm['m6a_sites']} ({mm['m6a_density_per_100nt']}/100 nt), "
                  f"pelanggaran jendela GC = {len(mm['gc_violations'])}")
                a("")

    if r.get("sensitivity"):
        s = r["sensitivity"]
        a("## Sensitivitas terhadap bobot")
        a("")
        a(f"{s['n']} set bobot dengan jitter +-{int(s['jitter_fraction']*100)}%: "
          f"skor komposit p05 = {s['composite_p05']}, median = {s['composite_median']}, "
          f"p95 = {s['composite_p95']}. Distribusi keputusan: {s['verdict_counts']}.")
        a("")

    if r.get("codon"):
        c = r["codon"]
        a("## Optimasi kodon ORF")
        a("")
        met = c.get("cost", {}).get("metrics", {})
        a(f"- Panjang ORF: {c.get('length_nt')} nt; iterasi {c.get('iterations')}; seed {c.get('seed')}")
        a(f"- Hasil: CAI = {met.get('cai')}, GC = {met.get('gc')}, "
          f"situs m6A = {met.get('m6a_sites')}, situs Tipe IIS = {met.get('restriction_sites')}, "
          f"motif dihindari = {met.get('avoid_motifs')}")
        bl = c.get("baseline_naive", {})
        if bl:
            a(f"- Baseline naif (kodon paling disukai saja): CAI = {bl.get('cai')}, "
              f"GC = {bl.get('gc')}, situs m6A = {bl.get('m6a_sites')}")
        a("")

    if r.get("primers"):
        p = r["primers"]
        a("## Oligo yang dirancang")
        a("")
        for key, title in (("divergent", "Primer divergen"), ("bsj", "Primer BSJ")):
            if key not in p:
                continue
            a(f"### {title}")
            a("")
            a("| # | Forward | Tm | Reverse | Tm | Ukuran produk |")
            a("|---:|---|---:|---|---:|---:|")
            for pair in p[key].get("pairs", [])[:5]:
                a(f"| {pair['rank']} | `{pair['left_seq']}` | {pair['left_tm']} | "
                  f"`{pair['right_seq']}` | {pair['right_tm']} | {pair['product_size']} |")
            a("")
        if "rnase_h" in p:
            a("### Probe RNase H")
            a("")
            a("| Posisi | Probe DNA (5'->3') | Tm | GC | Aksesibilitas |")
            a("|---:|---|---:|---:|---:|")
            for pr in p["rnase_h"].get("probes", []):
                a(f"| {pr['start_1based']}-{pr['end_1based']} | `{pr['probe_dna_5to3']}` | "
                  f"{pr['tm']} | {pr['gc']} | {pr['accessibility_plfold']:.4g} |")
            a("")
        if "ivt_reverse" in p:
            iv = p["ivt_reverse"]
            a(f"### Primer reverse IVT\n\n`{iv['with_modifications']}` — {iv['modification']}. "
              f"Tm {iv['tm']} C. Tujuan: {iv['purpose']}.")
            a("")

    if r.get("bsj_simulation"):
        b = r["bsj_simulation"]
        s = b["support"]
        a("## Filter statistik BSJ (simulasi bacaan)")
        a("")
        a(f"Parameter: {b['parameters']}")
        a("")
        a(f"- Kandidat BSJ: {s['n_candidates']}; lolos ambang posterior "
          f"{s['posterior_threshold']}: {s['n_passing']}")
        a(f"- Posterior gabungan sirkular: {s['joint_posterior_circular']:.6g} "
          f"(log10 BF gabungan {s['log10_joint_bayes_factor']})")
        a(f"- Terhadap asal bacaan simulasi: {s['confusion_on_simulated_origin']}")
        a("")
        a(s["note"])
        a("")

    if r.get("offtarget"):
        a("## Homologi off-target")
        a("")
        for k, v in r["offtarget"].items():
            if v.get("status") == "OK":
                n = v.get("n_hits", v.get("n_reference_sequences_with_hits"))
                a(f"- `{k}`: metode {v.get('method')}, status OK, jumlah hit = {n}")
            else:
                a(f"- `{k}`: **{v.get('status')}** — {v.get('reason')}")
        a("")

    a("## Catatan keterbatasan")
    a("")
    a("- Semua nilai di atas adalah PREDIKSI termodinamika dan heuristik sekuens, "
      "bukan pengukuran laboratorium. Tidak satu pun menggantikan uji RNase R, "
      "RT-PCR divergen dengan sekuensing Sanger, atau uji RNase H.")
    a("- Model energi ViennaRNA tidak memodelkan pseudoknot, sehingga domain IRES "
      "dan struktur tersier intron grup I tidak terwakili sepenuhnya; metrik "
      "integritas IRES di sini adalah proksi sekunder, bukan prediksi fungsi.")
    a("- Skor komposit adalah alat triase untuk mengurutkan kandidat sebelum "
      "sintesis, dengan bobot yang ditetapkan pipeline ini ([PERKIRAAN]) dan "
      "dapat diubah; lihat bagian sensitivitas.")
    return "\n".join(L)


def write_markdown(report: dict, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(markdown(report))
    return path


# --------------------------------------------------------------------------
# HTML mandiri
# --------------------------------------------------------------------------

_HTML = """<!doctype html>
<html lang="id"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Laporan Tahap 0 circRNA-HBsAg</title>
<style>
 :root{{--bg:#ffffff;--fg:#17201a;--mut:#5c6b61;--line:#dfe6e0;--card:#f6f8f6;
        --go:#1b6e3c;--cond:#8a6400;--no:#9b2226;--accent:#2a6b5f}}
 @media (prefers-color-scheme:dark){{:root:not([data-theme="light"]){{
   --bg:#11150f;--fg:#e8ede7;--mut:#9fb0a4;--line:#2b332c;--card:#181d18;
   --go:#6fd49a;--cond:#e5bd5a;--no:#f08a84;--accent:#7fd6c4}}}}
 :root[data-theme="dark"]{{--bg:#11150f;--fg:#e8ede7;--mut:#9fb0a4;--line:#2b332c;
   --card:#181d18;--go:#6fd49a;--cond:#e5bd5a;--no:#f08a84;--accent:#7fd6c4}}
 *{{box-sizing:border-box}}
 body{{margin:0;background:var(--bg);color:var(--fg);
   font:16px/1.65 ui-sans-serif,system-ui,-apple-system,"Segoe UI",sans-serif}}
 .wrap{{max-width:60rem;margin:0 auto;padding:2.5rem 1rem 5rem}}
 h1{{font-size:1.6rem;line-height:1.25;margin:0 0 .4rem}}
 h2{{font-size:1.15rem;margin:2.2rem 0 .6rem;padding-bottom:.3rem;border-bottom:1px solid var(--line)}}
 .label{{display:inline-block;padding:.2rem .55rem;border-radius:.4rem;
   background:var(--card);border:1px solid var(--line);font-size:.8rem;color:var(--mut)}}
 .verdict{{font-size:2rem;font-weight:700;margin:.3rem 0}}
 .GO{{color:var(--go)}} .COND{{color:var(--cond)}} .NO{{color:var(--no)}}
 table{{width:100%;border-collapse:collapse;margin:.6rem 0;font-size:.9rem}}
 th,td{{text-align:left;padding:.45rem .5rem;border-bottom:1px solid var(--line);
   vertical-align:top}}
 th{{color:var(--mut);font-weight:600;font-size:.8rem;text-transform:uppercase;
   letter-spacing:.03em}}
 td.num,th.num{{text-align:right;font-variant-numeric:tabular-nums}}
 code{{font:.85em ui-monospace,SFMono-Regular,Menlo,monospace;background:var(--card);
   padding:.1rem .3rem;border-radius:.25rem;word-break:break-all}}
 .bar{{height:.5rem;background:var(--card);border-radius:.25rem;overflow:hidden}}
 .bar>i{{display:block;height:100%;background:var(--accent)}}
 .meta{{color:var(--mut);font-size:.85rem}}
 .note{{background:var(--card);border-left:3px solid var(--accent);padding:.7rem .9rem;
   border-radius:.3rem;margin:1rem 0;font-size:.9rem}}
 .struct{{font:.75rem/1.4 ui-monospace,monospace;word-break:break-all;color:var(--mut)}}
 @media (max-width:40rem){{.wrap{{padding:1.5rem 1rem 4rem}} table{{font-size:.82rem}}}}
</style></head><body><div class="wrap">{body}</div></body></html>
"""


def html(report: dict) -> str:
    """Laporan HTML mandiri (tanpa dependensi eksternal) untuk diunduh."""
    import html as _h

    r = report["results"]
    dec = r.get("decision", {})
    m = r.get("measurements", {})
    cls = {"GO": "GO", "GO-BERSYARAT": "COND", "NO-GO": "NO"}.get(dec.get("verdict", ""), "")
    B: list[str] = []
    a = B.append
    a("<h1>Laporan Tahap 0 — kelayakan sirkularisasi circRNA-HBsAg (PIE)</h1>")
    a(f'<p><span class="label">{_h.escape(report["data_label"])}</span></p>')
    a(f'<p class="meta">Dihasilkan {report["generated_utc"]} UTC · circhbsag '
      f'{report["circhbsag_version"]} · ViennaRNA {report["tool_licenses"]["ViennaRNA"]["version"]} '
      f'· seed {report.get("seed")}</p>')
    if dec:
        a(f'<div class="verdict {cls}">{_h.escape(dec["verdict"])}</div>')
        a(f'<p class="meta">Skor komposit <b>{dec["composite_score"]:.3f}</b> '
          f'(GO &ge; {dec["thresholds"]["go"]}, bersyarat &ge; {dec["thresholds"]["conditional"]})</p>')
        a("<h2>Gerbang keras</h2><table><tr><th>Gerbang</th><th>Lulus</th><th>Nilai</th>"
          "<th>Batas</th><th>Alasan</th></tr>")
        for k, v in dec["hard_gates"].items():
            a(f"<tr><td><code>{k}</code></td><td>{'YA' if v['passed'] else '<b>TIDAK</b>'}</td>"
              f"<td><code>{_h.escape(str(v['value']))}</code></td>"
              f"<td><code>{_h.escape(str(v['limit']))}</code></td>"
              f"<td class='meta'>{_h.escape(v['why'])}</td></tr>")
        a("</table>")
        a("<h2>Kriteria berbobot</h2><table><tr><th>Kriteria</th><th class='num'>Bobot</th>"
          "<th class='num'>Subskor</th><th>Visual</th><th class='num'>Kontribusi</th></tr>")
        for c in dec["criteria"]:
            pct = int(round(c["subscore"] * 100))
            a(f"<tr><td><code>{c['name']}</code><div class='meta'>{_h.escape(c['rationale'])}</div></td>"
              f"<td class='num'>{c['weight']:.2f}</td><td class='num'>{c['subscore']:.3f}</td>"
              f"<td><div class='bar'><i style='width:{pct}%'></i></div></td>"
              f"<td class='num'>{c['weighted']:.3f}</td></tr>")
        a("</table>")
    if m.get("folds"):
        a("<h2>Lipatan RNA</h2><table><tr><th>Molekul</th><th class='num'>Suhu</th>"
          "<th class='num'>MFE</th><th class='num'>MFE/nt</th><th class='num'>Ensemble</th>"
          "<th class='num'>Keragaman</th></tr>")
        for k, v in m["folds"].items():
            a(f"<tr><td><code>{k}</code></td><td class='num'>{v['temperature_c']:g}&nbsp;C</td>"
              f"<td class='num'>{v['mfe_kcal']:.2f}</td><td class='num'>{v['mfe_per_nt']:.4f}</td>"
              f"<td class='num'>{v['ensemble_free_energy_kcal']:.2f}</td>"
              f"<td class='num'>{v['ensemble_diversity']:.2f}</td></tr>")
        a("</table>")
    if m.get("junction"):
        a("<h2>Sambungan balik</h2><table><tr><th>Lipatan</th>"
          "<th class='num'>P(tak-berpasangan)</th><th>Struktur jendela</th></tr>")
        for k, v in m["junction"].items():
            a(f"<tr><td><code>{k}</code></td><td class='num'>{v['mean_unpaired_probability']:.4f}</td>"
              f"<td class='struct'>{_h.escape(v['structure_window'])}</td></tr>")
        a("</table>")
    if m.get("accessibility"):
        a("<h2>Aksesibilitas RNAplfold</h2><table><tr><th>Suhu</th><th>Wilayah</th>"
          "<th class='num'>Rerata</th><th class='num'>Min</th></tr>")
        for t, regions in m["accessibility"].items():
            for name, v in regions.items():
                if v.get("mean") is None:
                    continue
                a(f"<tr><td>{t}</td><td><code>{name}</code></td>"
                  f"<td class='num'>{v['mean']:.5f}</td><td class='num'>{v['min']:.5f}</td></tr>")
        a("</table>")
    if m.get("motifs"):
        a("<h2>Motif</h2><table><tr><th>Molekul</th><th>Tipe IIS</th><th class='num'>Terminator T7</th>"
          "<th class='num'>m6A RRACH</th><th class='num'>/100 nt</th>"
          "<th class='num'>Pelanggaran GC</th></tr>")
        for which, mm in m["motifs"].items():
            rs = ", ".join(f"{k}={v}" for k, v in mm["restriction"].items() if v) or "tidak ada"
            a(f"<tr><td>{which}</td><td>{rs}</td><td class='num'>{len(mm['t7_terminator'])}</td>"
              f"<td class='num'>{mm['m6a_sites']}</td>"
              f"<td class='num'>{mm['m6a_density_per_100nt']}</td>"
              f"<td class='num'>{len(mm['gc_violations'])}</td></tr>")
        a("</table>")
    a('<div class="note">Semua nilai adalah prediksi termodinamika dan heuristik sekuens, '
      "bukan pengukuran laboratorium. Model energi ViennaRNA tidak memodelkan pseudoknot, "
      "sehingga integritas IRES dan struktur tersier intron grup I hanya terwakili sebagai "
      "proksi sekunder.</div>")
    return _HTML.format(body="\n".join(B))


def write_html(report: dict, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html(report))
    return path


def write_all(report: dict, out_dir: str | Path, stem: str = "tahap0") -> dict[str, Path]:
    out_dir = Path(out_dir)
    return {
        "json": write_json(report, out_dir / f"{stem}.json"),
        "markdown": write_markdown(report, out_dir / f"{stem}.md"),
        "html": write_html(report, out_dir / f"{stem}.html"),
    }
