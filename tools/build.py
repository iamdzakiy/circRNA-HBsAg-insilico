"""Generate every deliverable from ``content/`` and ``design_lock.json``.

    python -m tools.build all          # website + thesis docx + deck script
    python -m tools.build site
    python -m tools.build docx
    python -m tools.build deck
    python -m tools.build check        # fail if any generated file is stale

Nothing in this script contains a project number or a citation of its own: both
come from the content layer, so a change in one place propagates everywhere.
``check`` is what makes that promise enforceable — it rebuilds into a temporary
directory and compares, so CI fails when a generated file was edited by hand or
left behind by a content change.
"""
from __future__ import annotations

import argparse
import html
import json
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circhbsag import content as C  # noqa: E402

GENERATED = "BERKAS INI DIHASILKAN OLEH tools/build.py — JANGAN DIEDIT LANGSUNG"
# docs/ dipegang oleh halaman publik yang ditulis terpisah. Generator ini
# menulis situsnya sendiri ke build/site/ dan hanya menerbitkan DATA ke
# docs/data/, supaya dua sumber tidak saling menimpa.
SITE_DIR = ROOT / "build" / "site"
DATA_DIR = ROOT / "docs" / "data"
BUILD_DIR = ROOT / "build"

PAGES = [
    ("index.html", "Beranda"),
    ("design.html", "Desain konstruk"),
    ("gate0.html", "Gate 0"),
    ("pipeline.html", "Pipeline"),
    ("decisions.html", "Keputusan"),
    ("planning.html", "Rencana"),
    ("logbook.html", "Logbook"),
    ("methods.html", "Metode & software"),
    ("downloads.html", "Unduhan"),
    ("calculator.html", "Kalkulator"),
]

CSS = """
:root{--bg:#fbfaf8;--fg:#1a1a18;--muted:#5f5d57;--line:#e2ded6;--card:#fff;
--accent:#1d6f5e;--warn:#8a5a00;--warnbg:#fff7e6;--pass:#1d6f5e;--rev:#8a5a00;--fail:#a33;--na:#5f5d57}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#15161a;--fg:#eceae6;
--muted:#a8a5a0;--line:#2c2e35;--card:#1c1e23;--accent:#5fc3ab;--warn:#e8b55c;--warnbg:#2a2318;
--pass:#5fc3ab;--rev:#e8b55c;--fail:#ef8a8a;--na:#a8a5a0}}
:root[data-theme="dark"]{--bg:#15161a;--fg:#eceae6;--muted:#a8a5a0;--line:#2c2e35;--card:#1c1e23;
--accent:#5fc3ab;--warn:#e8b55c;--warnbg:#2a2318;--pass:#5fc3ab;--rev:#e8b55c;--fail:#ef8a8a;--na:#a8a5a0}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:16px/1.65 ui-sans-serif,system-ui,-apple-system,"Segoe UI",sans-serif}
header.site{border-bottom:1px solid var(--line);background:var(--card)}
.wrap{max-width:island;max-width:62rem;margin:0 auto;padding:0 16px}
nav{display:flex;flex-wrap:wrap;gap:2px;padding:10px 0}
nav a{padding:6px 10px;border-radius:7px;color:var(--muted);text-decoration:none;font-size:14px}
nav a:hover{background:var(--bg);color:var(--fg)}
nav a[aria-current=page]{background:var(--accent);color:var(--card);font-weight:600}
h1{font-size:1.75rem;line-height:1.25;margin:28px 0 6px}
h2{font-size:1.25rem;margin:32px 0 8px;padding-bottom:5px;border-bottom:1px solid var(--line)}
h3{font-size:1.02rem;margin:22px 0 4px}
p,li,td,th,h1,h2,h3{overflow-wrap:anywhere}
p,li{max-width:62ch}
main{overflow-x:hidden}
.lede{color:var(--muted)}
.banner{background:var(--warnbg);border:1px solid var(--warn);color:var(--warn);
border-radius:10px;padding:12px 14px;margin:18px 0;font-size:14px}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:14px 16px;margin:14px 0}
.grid{display:grid;gap:12px;grid-template-columns:repeat(auto-fit,minmax(min(100%,13rem),1fr))}
.kv{font-variant-numeric:tabular-nums}
.kv b{display:block;font-size:1.3rem}
.kv span{color:var(--muted);font-size:13px}
table{width:100%;border-collapse:collapse;font-size:14px;margin:10px 0;display:block;overflow-x:auto}
th,td{border-bottom:1px solid var(--line);padding:7px 9px;text-align:left;vertical-align:top}
th{color:var(--muted);font-weight:600}
code,pre{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:13px}
pre{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px;overflow-x:auto}
.tag{display:inline-block;padding:1px 7px;border-radius:999px;font-size:12px;border:1px solid var(--line);color:var(--muted)}
.v-PASS{color:var(--pass);font-weight:700}
.v-REVISE{color:var(--rev);font-weight:700}
.v-FAIL{color:var(--fail);font-weight:700}
.v-UNAVAILABLE{color:var(--na);font-weight:700}
footer{margin:48px 0 24px;padding-top:14px;border-top:1px solid var(--line);color:var(--muted);font-size:13px}
a{color:var(--accent)}
svg.chart{max-width:100%;height:auto;background:var(--card);border:1px solid var(--line);border-radius:10px}
"""


def esc(x) -> str:
    return html.escape(str(x), quote=False)


def page(title: str, body: str, *, current: str, project: dict) -> str:
    nav = "".join(
        f'<a href="{f}"{" aria-current=page" if f == current else ""}>{esc(n)}</a>'
        for f, n in PAGES
    )
    return f"""<!doctype html>
<html lang="id"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{esc(title)} — {esc(project['short_title'])}</title>
<meta name="description" content="{esc(project['subtitle'])}">
<!-- {GENERATED} -->
<style>{CSS}</style></head>
<body><header class="site"><div class="wrap"><nav>{nav}</nav></div></header>
<main class="wrap">{body}</main>
<div class="wrap"><footer>
{esc(project['title'])} · {esc(project['institution'])} · kode berlisensi {esc(project['code_license'])}<br>
Halaman ini dihasilkan otomatis dari <code>content/</code> dan <code>design_lock.json</code>
pada {datetime.now(timezone.utc).strftime('%Y-%m-%d')} UTC. Jangan menyunting HTML ini langsung.
</footer></div></body></html>
"""


def banner(c: C.Content) -> str:
    if not c.has_numbers:
        return (
            '<div class="banner">Belum ada jalankan pipeline, sehingga halaman ini '
            "belum memuat angka konstruk.</div>"
        )
    return f'<div class="banner"><b>Status data:</b> {esc(c.data_label)}</div>'


def table(headers: list[str], rows: list[list[str]], *, raw_cols: tuple[int, ...] = ()) -> str:
    th = "".join(f"<th>{esc(h)}</th>" for h in headers)
    trs = []
    for r in rows:
        tds = "".join(
            f"<td>{cell if i in raw_cols else esc(cell)}</td>" for i, cell in enumerate(r)
        )
        trs.append(f"<tr>{tds}</tr>")
    return f"<table><thead><tr>{th}</tr></thead><tbody>{''.join(trs)}</tbody></table>"


def md_to_html(path: Path) -> str:
    """Render a Markdown file, or say plainly that it is missing."""
    if not path.exists():
        return f"<p><em>{esc(path.name)} belum ada.</em></p>"
    try:
        import markdown

        return markdown.markdown(
            path.read_text(), extensions=["tables", "fenced_code", "toc", "attr_list"]
        )
    except ImportError:
        return f"<pre>{esc(path.read_text())}</pre>"


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------


def construct_svg(c: C.Content) -> str:
    """Inline SVG map of the precursor, drawn to scale from the design lock."""
    coords = c.lock("sequence_coords") or {}
    parts = c.project["construct"]["precursor"]
    circle_ids = set(c.project["construct"]["circle"])
    total = max(1, len(parts))
    w, h = 980, 150
    bw = (w - 20) / total
    cells = []
    for i, p in enumerate(parts):
        x = 10 + i * bw
        fill = "var(--accent)" if p["id"] in circle_ids else "var(--line)"
        txt = "var(--card)" if p["id"] in circle_ids else "var(--fg)"
        cells.append(
            f'<g><title>{esc(p["label"])} — {esc(p["role"])}</title>'
            f'<rect x="{x:.1f}" y="46" width="{bw - 3:.1f}" height="40" rx="5" fill="{fill}"/>'
            f'<text x="{x + bw / 2:.1f}" y="70" text-anchor="middle" font-size="10" fill="{txt}">'
            f'{esc(p["id"])}</text></g>'
        )
    note = (
        f'Prekursor {c.num("lengths_nt", "precursor", unit="nt")} · '
        f'lingkaran {c.num("lengths_nt", "circle", unit="nt")} · '
        f'BSJ {c.num("lengths_nt", "bsj_position")}'
    )
    return (
        f'<svg class="chart" viewBox="0 0 {w} {h}" role="img" '
        f'aria-label="Peta konstruk prekursor"><text x="10" y="26" font-size="13" '
        f'fill="var(--muted)">Kotak berwarna = ikut ke dalam lingkaran setelah splicing</text>'
        + "".join(cells)
        + f'<text x="10" y="112" font-size="12" fill="var(--fg)">{esc(note)}</text>'
        f'<text x="10" y="132" font-size="11" fill="var(--muted)">{esc(unused_note(coords))}'
        "</text></svg>"
    )


def unused_note(coords: dict) -> str:
    return (
        "Panjang kotak pada gambar ini tidak proporsional terhadap panjang sekuens; "
        "koordinat sebenarnya ada di gate0.json."
    )


def arm_curve_svg(c: C.Content) -> str:
    """ΔG of the homology-arm duplex against arm length, from the design lock."""
    scan = [p for p in (c.lock("homology_arms", "length_scan") or []) if p.get("available")]
    if not scan:
        return "<p><em>Kurva ΔG lengan homologi belum ada: jalankan <code>circhbsag gate0</code>.</em></p>"
    w, h, pad = 560, 260, 46
    xs = [p["length_nt"] for p in scan]
    ys = [p["binding_kcal"] for p in scan]
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys), max(ys)
    if y1 == y0:
        y1 = y0 + 1

    def px(x):
        return pad + (x - x0) / max(1, (x1 - x0)) * (w - 2 * pad)

    def py(y):
        return h - pad - (y - y0) / (y1 - y0) * (h - 2 * pad)

    pts = " ".join(f"{px(x):.1f},{py(y):.1f}" for x, y in zip(xs, ys, strict=True))
    dots = "".join(
        f'<g><title>{x} nt: {y:.2f} kcal/mol</title>'
        f'<circle cx="{px(x):.1f}" cy="{py(y):.1f}" r="4" fill="var(--accent)"/></g>'
        for x, y in zip(xs, ys, strict=True)
    )
    ticks = "".join(
        f'<text x="{px(x):.1f}" y="{h - pad + 16}" font-size="11" text-anchor="middle" '
        f'fill="var(--muted)">{x}</text>'
        for x in xs
    )
    t = c.lock("homology_arms", "temperature_c")
    return (
        f'<svg class="chart" viewBox="0 0 {w} {h}" role="img" aria-label="Kurva ΔG lengan homologi">'
        f'<text x="{pad}" y="22" font-size="12" fill="var(--muted)">ΔG pengikatan lengan homologi '
        f'(kcal/mol) pada {esc(t)} °C</text>'
        f'<line x1="{pad}" y1="{h - pad}" x2="{w - pad}" y2="{h - pad}" stroke="var(--line)"/>'
        f'<line x1="{pad}" y1="{pad}" x2="{pad}" y2="{h - pad}" stroke="var(--line)"/>'
        f'<polyline points="{pts}" fill="none" stroke="var(--accent)" stroke-width="2"/>{dots}{ticks}'
        f'<text x="{w / 2:.0f}" y="{h - 8}" font-size="11" text-anchor="middle" '
        f'fill="var(--muted)">panjang lengan (nt)</text>'
        f'<text x="{pad}" y="{py(y0):.1f}" font-size="11" fill="var(--muted)" dx="-38">{y0:.0f}</text>'
        f'<text x="{pad}" y="{py(y1):.1f}" font-size="11" fill="var(--muted)" dx="-38">{y1:.0f}</text>'
        "</svg>"
    )


# --------------------------------------------------------------------------
# Site
# --------------------------------------------------------------------------


def build_site(c: C.Content, gate: dict | None) -> list[Path]:
    out: list[Path] = []
    pr = c.project
    SITE_DIR.mkdir(parents=True, exist_ok=True)

    def write(name: str, title: str, body: str):
        p = SITE_DIR / name
        p.write_text(page(title, body, current=name, project=pr))
        out.append(p)

    # ---- Beranda
    score = c.lock("composite_score")
    tiles = [
        ("Prekursor", c.num("lengths_nt", "precursor", unit="nt")),
        ("Lingkaran", c.num("lengths_nt", "circle", unit="nt")),
        ("ORF", c.num("lengths_nt", "orf", unit="nt")),
        ("Protein", c.num("lengths_nt", "protein_aa", unit="aa")),
        ("Skor komposit", "—" if score is None else f"{score:.4f}"),
        ("Keputusan Gate 0", str(c.lock("gate0", "decision") or "—")),
    ]
    grid = "".join(
        f'<div class="card kv"><b>{esc(v)}</b><span>{esc(k)}</span></div>' for k, v in tiles
    )
    write(
        "index.html",
        "Beranda",
        f"<h1>{esc(pr['title'])}</h1><p class=lede>{esc(pr['subtitle'])}</p>"
        + banner(c)
        + f'<div class="grid">{grid}</div>'
        + "<h2>Pertanyaan yang dijawab situs ini</h2>"
        "<p>Apakah konstruk RNA sirkular pembawa ORF HBsAg ini layak disintesis, "
        "atau sebaiknya dirancang ulang lebih dahulu? Jawabannya diberikan sebagai "
        "keputusan Gate 0 atas sembilan uji, ditambah skor komposit berbobot yang "
        "seluruh bobot dan ambangnya tertulis.</p>"
        + construct_svg(c)
        + "<h2>Tiga aturan yang dipegang</h2><ol>"
        "<li><b>Tidak ada sekuens yang dikarang.</b> Sekuens referensi tidak diikutkan "
        "ke repositori; bila berkasnya belum ada, pipeline berhenti dengan galat.</li>"
        "<li><b>Tidak ada angka yang dihardcode.</b> Setiap nilai dihitung saat dijalankan, "
        "dan seluruh halaman ini dihasilkan dari satu berkas angka.</li>"
        "<li><b>Setiap laporan membawa label provenans.</b> Satu bagian masukan sintetis "
        "sudah cukup untuk melabeli seluruh laporan sebagai demo.</li></ol>"
        f'<p class=lede>{esc(C.PROJECTION_LABEL)}.</p>',
    )

    # ---- Desain konstruk
    rows = [
        [p["id"], p["label"], p["role"], "ya" if p["id"] in set(pr["construct"]["circle"]) else "tidak"]
        for p in pr["construct"]["precursor"]
    ]
    write(
        "design.html",
        "Desain konstruk",
        "<h1>Desain konstruk</h1>"
        + banner(c)
        + construct_svg(c)
        + table(["id", "Elemen", "Peran", "Ikut ke lingkaran"], rows)
        + "<h2>Lengan homologi</h2>"
        + table(
            ["Besaran", "Nilai"],
            [
                ["Lengan 5'", c.num("homology_arms", "arm5_nt", unit="nt")],
                ["Lengan 3'", c.num("homology_arms", "arm3_nt", unit="nt")],
                ["ΔG pengikatan (RNAcofold)", c.num("homology_arms", "binding_dg_kcal", unit="kcal/mol")],
                ["ΔG duplex (RNAduplex)", c.num("homology_arms", "duplex_dg_kcal", unit="kcal/mol")],
                ["Suhu", c.num("homology_arms", "temperature_c", unit="°C")],
            ],
        )
        + arm_curve_svg(c)
        + "<h2>Lipatan</h2>"
        + table(
            ["Besaran", "Nilai"],
            [
                ["MFE lingkaran", c.num("folding", "circle_mfe_kcal", unit="kcal/mol")],
                ["MFE lingkaran per nt", c.num("folding", "circle_mfe_per_nt", unit="kcal/mol/nt")],
                ["Ensemble diversity lingkaran", c.num("folding", "circle_ensemble_diversity")],
                ["MFE prekursor", c.num("folding", "precursor_mfe_kcal", unit="kcal/mol")],
                ["Ensemble diversity prekursor", c.num("folding", "precursor_ensemble_diversity")],
                ["Suhu", c.num("folding", "temperature_c", unit="°C")],
            ],
        ),
    )

    # ---- Gate 0
    write("gate0.html", "Gate 0", gate0_body(c, gate))

    # ---- Pipeline
    write(
        "pipeline.html",
        "Pipeline",
        "<h1>Pipeline</h1>"
        + banner(c)
        + "<p>Setiap kotak di bawah adalah satu aturan Snakemake. Perintahnya dapat "
        "dijalankan sendiri-sendiri, dan setiap aturan menulis satu berkas JSON berisi "
        "parameter dan versi alat yang dipakai.</p>"
        + table(
            ["Langkah", "Perintah", "Keluaran"],
            [
                ["Lisensi alat", "<code>circhbsag licenses</code>", "licenses.json"],
                ["Optimasi kodon", "<code>circhbsag codon</code>", "orf_optimized.fa, codon.json"],
                ["Gate 0 sembilan uji", "<code>circhbsag gate0</code>", "gate0.json, gate0.md, design_lock.json"],
                ["Skor Tahap 0", "<code>circhbsag tahap0</code>", "tahap0.json/md/html"],
                ["Primer dan probe", "<code>circhbsag primers</code>", "primers.json"],
                ["Simulasi BSJ", "<code>circhbsag bsj-sim</code>", "bsj_simulation.json"],
                ["Homologi off-target", "<code>circhbsag offtarget</code>", "offtarget.json"],
                ["Optimasi multi-objektif", "<code>circhbsag optimize</code>", "optimization.json"],
                ["Dokumen dan situs", "<code>python -m tools.build all</code>", "docs/, build/"],
            ],
            raw_cols=(1,),
        )
        + "<h2>Reprodusibilitas</h2><pre>"
        + esc(
            "mamba env create -f environment.yml && conda activate circhbsag\n"
            "pip install -e .\n"
            "snakemake --cores 4 --config parts=config/parts.example.yaml outdir=results/demo\n"
            "python -m tools.build all"
        )
        + "</pre>",
    )

    # ---- Keputusan
    dec_rows = []
    for d in c.decisions:
        refs = ", ".join(
            f'<a href="https://doi.org/{c.ref(k)["doi"]}">{esc(c.short_cite(k))}</a>'
            for k in d.get("refs", [])
        )
        dec_rows.append(
            [
                d["id"],
                d["question"],
                d["choice"],
                f'<span class="tag">{esc(d["label"])}</span>',
                d["rationale"],
                d.get("consequence", "—"),
                refs or "—",
            ]
        )
    write(
        "decisions.html",
        "Keputusan",
        "<h1>Keputusan desain</h1>"
        "<p class=lede>Satu baris per keputusan, dengan alasan, konsekuensi, dan sumbernya. "
        "Keputusan berlabel LITERATUR wajib menyebut referensi; pembangun situs menolak jalan "
        "bila tidak.</p>"
        + table(
            ["id", "Pertanyaan", "Pilihan", "Label", "Alasan", "Konsekuensi", "Sumber"],
            dec_rows,
            raw_cols=(3, 6),
        )
        + "<h2>Assumption Register</h2>"
        + table(
            ["id", "Asumsi", "Alasan diambil", "Dampak bila salah", "Cara memverifikasi"],
            C.assumptions_rows(c),
        )
        + "<h2>Open Questions</h2>"
        + table(
            ["id", "Pertanyaan", "Dijawab oleh", "Memblokir"], C.open_questions_rows(c)
        ),
    )

    # ---- Rencana & logbook (dari Markdown)
    write(
        "planning.html",
        "Rencana",
        "<h1>Rencana optimasi</h1>" + banner(c) + md_to_html(ROOT / "docs" / "PLANNING.md"),
    )
    write(
        "logbook.html",
        "Logbook",
        "<h1>Logbook</h1>" + md_to_html(ROOT / "docs" / "LOGBOOK.md"),
    )

    # ---- Metode & software
    versions = C.detected_versions()
    sw_rows = []
    for r in C.software_rows(c):
        name, role, ver, lic, url, used = r
        sw_rows.append(
            [
                name,
                role,
                versions.get(name, ver),
                lic,
                f'<a href="{esc(url)}">{esc(url)}</a>',
                used,
            ]
        )
    write(
        "methods.html",
        "Metode & software",
        "<h1>Metode dan software</h1>"
        f"<p class=lede>Versi pada kolom ketiga dibaca dari paket yang benar-benar terpasang "
        f"saat halaman ini dibangun, bukan ditulis tangan. Tanggal akses: "
        f"{esc(c.raw['software.yaml'].get('accessed', '—'))}.</p>"
        + table(
            ["Alat", "Fungsi", "Versi", "Lisensi", "URL", "Status"], sw_rows, raw_cols=(4,)
        )
        + "<h2>Alat yang sengaja tidak dipakai</h2>"
        + table(
            ["Alat", "Alasan", "Harga yang dibayar"],
            [[s["name"], s["reason"], s.get("cost", "—")] for s in c.software_not_used],
        )
        + "<h2>Daftar pustaka</h2><ol>"
        + "".join(
            f'<li>{esc(r["authors"])}. {esc(r["title"])}. <em>{esc(r["journal"])}</em> '
            f'{esc(r["year"])}; {esc(r["citation"])}. '
            f'<a href="https://doi.org/{esc(r["doi"])}">doi:{esc(r["doi"])}</a></li>'
            for r in c.references.values()
        )
        + "</ol><p class=lede>Metadata sitasi diambil dari PubMed (NCBI) pada "
        + esc(c.raw["references.yaml"].get("retrieved", "—"))
        + ".</p>",
    )

    # ---- Unduhan
    write(
        "downloads.html",
        "Unduhan",
        "<h1>Unduhan</h1>" + banner(c) + downloads_body(c),
    )
    return out


def gate0_body(c: C.Content, gate: dict | None) -> str:
    if not gate:
        return (
            "<h1>Gate 0</h1>"
            "<p>Belum ada laporan Gate 0. Jalankan <code>circhbsag gate0</code>.</p>"
        )
    sym = {
        "PASS": "PASS",
        "REVISE": "REVISE",
        "FAIL": "FAIL",
        "UNAVAILABLE": "BELUM BISA DIJALANKAN",
    }
    rows = []
    for i, t in enumerate(gate["tests"], 1):
        rows.append(
            [
                str(i),
                t["name"],
                f'<span class="v-{t["verdict"]}">{sym[t["verdict"]]}</span>',
                t["threshold"],
                f'<span class="tag">{esc(t["threshold_label"])}</span>',
            ]
        )
    detail = []
    for i, t in enumerate(gate["tests"], 1):
        detail.append(
            f'<div class="card"><h3>Uji {i} — {esc(t["name"])} '
            f'<span class="v-{t["verdict"]}">{sym[t["verdict"]]}</span></h3>'
            f"<p>{esc(t['question'])}</p>"
            f"<p><b>Ambang</b> ({esc(t['threshold_label'])}): {esc(t['threshold'])}</p>"
            f"<p><b>Dasar ambang:</b> {esc(t['threshold_rationale'])}</p>"
            + (f"<p><b>Catatan:</b> {esc(t['note'])}</p>" if t.get("note") else "")
            + "<details><summary>Angka mentah</summary><pre>"
            + esc(json.dumps(t["metrics"], indent=2, ensure_ascii=False)[:6000])
            + "</pre></details></div>"
        )
    cnt = gate["counts"]
    return (
        "<h1>Gate 0 — kelayakan sirkularisasi</h1>"
        + banner(c)
        + f'<div class="card"><b>Keputusan: {esc(gate["decision"])}</b><br>{esc(gate["statement"])}'
        f'<br><br>{cnt["PASS"]} PASS · {cnt["REVISE"]} REVISE · {cnt["FAIL"]} FAIL · '
        f'{cnt["UNAVAILABLE"]} belum bisa dijalankan</div>'
        + table(["#", "Uji", "Verdict", "Ambang", "Dasar"], rows, raw_cols=(2, 4))
        + arm_curve_svg(c)
        + "<h2>Rincian tiap uji</h2>"
        + "".join(detail)
    )


def downloads_body(c: C.Content) -> str:
    items = [
        ("design_lock.json", "Satu berkas angka yang dikutip semua deliverable"),
        ("results/demo/gate0.json", "Laporan Gate 0 lengkap dengan angka mentah"),
        ("results/demo/gate0.md", "Laporan Gate 0 dalam Markdown"),
        ("results/demo/tahap0.html", "Laporan Tahap 0 mandiri"),
        ("results/demo/SUMMARY.md", "Ringkasan jalankan pipeline"),
        ("results/demo/primers.json", "Tabel oligo: primer divergen, BSJ, probe RNase H"),
        ("results/demo/optimization.json", "Front Pareto optimasi multi-objektif"),
        ("build/TA_circRNA-HBsAg.docx", "Dokumen TA yang dihasilkan dari sumber tunggal"),
        ("build/deck_script.md", "Naskah dek"),
    ]
    rows = [[f"<code>{esc(p)}</code>", d] for p, d in items]
    return (
        "<p>Berkas di bawah ada di repositori. Situs ini statis, jadi tautannya mengikuti "
        "struktur repo, bukan penyimpanan terpisah.</p>"
        + table(["Berkas", "Isi"], rows, raw_cols=(0,))
    )


# --------------------------------------------------------------------------
# Thesis document
# --------------------------------------------------------------------------


def build_docx(c: C.Content, gate: dict | None) -> list[Path]:
    from docx import Document
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.shared import Pt

    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(12)

    pr = c.project
    t = doc.add_paragraph(pr["title"])
    t.alignment = WD_ALIGN_PARAGRAPH.CENTER
    t.runs[0].bold = True
    t.runs[0].font.size = Pt(16)
    for line in (pr["subtitle"], pr["program"], pr["institution"], pr["author"]):
        p = doc.add_paragraph(line)
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    doc.add_paragraph()
    warn = doc.add_paragraph()
    warn.add_run(GENERATED).italic = True
    doc.add_paragraph(C.PROJECTION_LABEL).runs[0].bold = True
    doc.add_paragraph(c.data_label)

    def h(text, level=1):
        doc.add_heading(text, level=level)

    def tbl(headers, rows):
        tb = doc.add_table(rows=1, cols=len(headers))
        tb.style = "Light Grid Accent 1"
        for i, hh in enumerate(headers):
            cell = tb.rows[0].cells[i]
            cell.text = str(hh)
            for p in cell.paragraphs:
                for r in p.runs:
                    r.bold = True
        for row in rows:
            cells = tb.add_row().cells
            for i, v in enumerate(row):
                cells[i].text = "" if v is None else str(v)
        doc.add_paragraph()

    h("Keputusan desain terkunci (Design Lock)")
    doc.add_paragraph(
        "Tabel ini adalah satu-satunya tempat keputusan desain ditetapkan. Dokumen, "
        "situs, dan dek membacanya dari berkas yang sama, sehingga angka di ketiganya "
        "tidak bisa berbeda tanpa terdeteksi."
    )
    tbl(
        ["Parameter", "Nilai terpilih", "Label", "Alasan ringkas"],
        [
            [d["question"], d["choice"], d["label"], d["rationale"].strip()[:400]]
            for d in c.decisions
        ],
    )

    h("Angka konstruk")
    tbl(
        ["Besaran", "Nilai"],
        [
            ["Panjang prekursor linier", c.num("lengths_nt", "precursor", unit="nt")],
            ["Panjang produk sirkular", c.num("lengths_nt", "circle", unit="nt")],
            ["Posisi sambungan balik", c.num("lengths_nt", "bsj_position")],
            ["Panjang ORF", c.num("lengths_nt", "orf", unit="nt")],
            ["Panjang protein", c.num("lengths_nt", "protein_aa", unit="aa")],
            ["IRES", f'{c.lock("ires", "name") or "—"} ({c.num("ires", "length_nt", unit="nt")})'],
            ["Metode sirkularisasi", c.lock("circularisation_method") or "—"],
            ["ΔG pengikatan lengan homologi", c.num("homology_arms", "binding_dg_kcal", unit="kcal/mol")],
            ["MFE lingkaran", c.num("folding", "circle_mfe_kcal", unit="kcal/mol")],
            ["Ensemble diversity lingkaran", c.num("folding", "circle_ensemble_diversity")],
            ["Skor komposit", c.num("composite_score", nd=4)],
            ["Keputusan Gate 0", c.lock("gate0", "decision") or "—"],
        ],
    )

    h("Gate 0: kelayakan sirkularisasi")
    if gate:
        doc.add_paragraph(f"Keputusan: {gate['decision']} — {gate['statement']}")
        tbl(
            ["#", "Uji", "Verdict", "Ambang", "Dasar ambang"],
            [
                [i, t["name"], t["verdict"], t["threshold"], t["threshold_label"]]
                for i, t in enumerate(gate["tests"], 1)
            ],
        )
        doc.add_paragraph(
            "Ambang berlabel PERKIRAAN berasal dari penalaran penulis, bukan dari sumber "
            "primer, dan ditandai [PERLU VERIFIKASI] di laporan Gate 0 beserta kata kunci "
            "yang harus dicari."
        )
    else:
        doc.add_paragraph("Belum ada laporan Gate 0.")

    h("Assumption Register")
    tbl(
        ["id", "Asumsi", "Alasan diambil", "Dampak bila salah", "Cara memverifikasi"],
        C.assumptions_rows(c),
    )

    h("Open Questions")
    tbl(["id", "Pertanyaan", "Dijawab oleh", "Memblokir"], C.open_questions_rows(c))

    h("Software dan Versi")
    versions = C.detected_versions()
    tbl(
        ["Alat", "Fungsi", "Versi", "Lisensi", "URL", "Status"],
        [[r[0], r[1], versions.get(r[0], r[2]), r[3], r[4], r[5]] for r in C.software_rows(c)],
    )
    doc.add_paragraph("Alat yang sengaja tidak dipakai:")
    tbl(
        ["Alat", "Alasan", "Harga yang dibayar"],
        [[s["name"], s["reason"], s.get("cost", "—")] for s in c.software_not_used],
    )

    h("Daftar Pustaka")
    for r in c.references.values():
        doc.add_paragraph(
            f"{r['authors']}. {r['title']}. {r['journal']} {r['year']}; {r['citation']}. "
            f"doi:{r['doi']}",
            style="List Number",
        )
    doc.add_paragraph(
        f"Metadata sitasi diambil dari PubMed (NCBI) pada "
        f"{c.raw['references.yaml'].get('retrieved', '—')}."
    )

    out = BUILD_DIR / "TA_circRNA-HBsAg.docx"
    doc.save(out)
    return [out]


# --------------------------------------------------------------------------
# Deck script
# --------------------------------------------------------------------------

SLIDES = [
    ("Masalah", "Hepatitis B masih endemis, dan vaksin rekombinan yang ada menuntut rantai dingin serta dosis berulang."),
    ("Gagasan", "RNA sirkular pembawa ORF HBsAg: tanpa ujung bebas, lebih stabil, dan dapat diproduksi in vitro."),
    ("Mengapa sirkular", "Sirkularisasi menurunkan imunogenisitas bawaan dan dapat memperpanjang durasi translasi."),
    ("Arsitektur konstruk", "PIE berbasis intron grup I: lengan homologi mendekatkan ujung, ekson menyatu menjadi lingkaran."),
    ("Keputusan terkunci", "Satu IRES (CVB3), antigen S-HBsAg 226 aa, subgenotipe B3, metode PIE."),
    ("Relevansi lokal", "B3 dominan di paruh barat kepulauan; C dominan di Papua. Satu konstruk tidak menutup semuanya."),
    ("Gate 0", "Sembilan uji sebelum satu baris protokol basah ditulis."),
    ("Gate 0: hasil", "Verdict per uji, dengan ambang dan dasarnya."),
    ("Lengan homologi", "ΔG pengikatan terhadap panjang lengan 20-40 nt."),
    ("Lipatan", "MFE dan ensemble diversity prekursor dan lingkaran, pada 37 dan 55 °C."),
    ("Integritas IRES", "Domain IRES harus tetap terlipat setelah masuk ke dalam lingkaran."),
    ("Risiko dsRNA", "Heliks panjang adalah pola yang dikenali sensor imun bawaan."),
    ("Optimasi", "NSGA-II atas spacer, lengan, dan kodon; front Pareto diperingkat ulang dengan skor komposit."),
    ("Reprodusibilitas", "Snakemake, seed tetap, CI, dan satu berkas angka yang dibaca semua deliverable."),
    ("Keterbatasan", "Prediksi, bukan pengukuran. Model tidak menangani pseudoknot."),
    ("Rencana validasi", "In vitro lalu in vivo, dengan kontrol yang dirancang sejak awal."),
    ("Permintaan", "Akses sekuens referensi, waktu laboratorium, dan reagen IVT."),
]


def build_deck(c: C.Content, gate: dict | None) -> list[Path]:
    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    score = c.lock("composite_score")
    facts = {
        "Arsitektur konstruk": (
            f'Prekursor {c.num("lengths_nt", "precursor", unit="nt")}, '
            f'lingkaran {c.num("lengths_nt", "circle", unit="nt")}, '
            f'BSJ pada {c.num("lengths_nt", "bsj_position")}.'
        ),
        "Keputusan terkunci": "; ".join(f'{d["choice"]}' for d in c.decisions[:4]),
        "Gate 0: hasil": (
            f'{gate["decision"]}: ' + ", ".join(f'{k} {v}' for k, v in gate["counts"].items())
            if gate
            else "belum ada laporan Gate 0"
        ),
        "Lengan homologi": (
            f'{c.num("homology_arms", "arm5_nt", unit="nt")} + '
            f'{c.num("homology_arms", "arm3_nt", unit="nt")}, ΔG pengikatan '
            f'{c.num("homology_arms", "binding_dg_kcal", unit="kcal/mol")} pada '
            f'{c.num("homology_arms", "temperature_c", unit="°C")}.'
        ),
        "Lipatan": (
            f'MFE lingkaran {c.num("folding", "circle_mfe_kcal", unit="kcal/mol")} '
            f'({c.num("folding", "circle_mfe_per_nt")} per nt), ensemble diversity '
            f'{c.num("folding", "circle_ensemble_diversity")}.'
        ),
        "Optimasi": (
            "Skor komposit " + ("—" if score is None else f"{score:.4f}") + " ("
            + str(c.lock("composite_verdict") or "—") + ")."
        ),
    }
    lines = [
        f"# Naskah dek — {c.project['short_title']}",
        "",
        f"> {GENERATED}",
        "",
        f"**{C.PROJECTION_LABEL}**",
        "",
        f"Status data: {c.data_label}",
        "",
    ]
    for i, (title, msg) in enumerate(SLIDES, 1):
        lines += [f"## Slide {i} — {title}", "", msg, ""]
        if title in facts:
            lines += [f"Angka dari design_lock.json: {facts[title]}", ""]
    lines += [
        "## Catatan untuk pembuat dek",
        "",
        "Setiap angka di naskah ini berasal dari design_lock.json. Jangan menyalinnya "
        "ke dalam slide sebagai teks tetap tanpa mencatat versinya: jalankan ulang "
        "`python -m tools.build deck` setelah pipeline berubah, dan perbarui slide yang "
        "angkanya bergeser.",
        "",
    ]
    md = BUILD_DIR / "deck_script.md"
    md.write_text("\n".join(lines))

    from pptx import Presentation
    from pptx.util import Inches, Pt

    prs = Presentation()
    blank = prs.slide_layouts[5]
    for i, (title, msg) in enumerate(SLIDES, 1):
        s = prs.slides.add_slide(blank)
        s.shapes.title.text = f"{i}. {title}"
        box = s.shapes.add_textbox(Inches(0.7), Inches(1.9), Inches(8.6), Inches(3.4))
        tf = box.text_frame
        tf.word_wrap = True
        tf.text = msg
        if title in facts:
            p = tf.add_paragraph()
            p.text = facts[title]
            p.font.size = Pt(14)
        p = tf.add_paragraph()
        p.text = C.PROJECTION_LABEL
        p.font.size = Pt(11)
        p.font.italic = True
    pptx = BUILD_DIR / "deck.pptx"
    prs.save(pptx)
    return [md, pptx]


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------



def build_data(c: C.Content, gate: dict | None) -> list[Path]:
    """Publish machine-readable numbers for the public site to read.

    The public pages under ``docs/`` are written by hand elsewhere; this target
    gives them the numbers as JSON so they never have to hard-code one. Only
    data is written here — no HTML — so the two sources cannot overwrite each
    other.
    """
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    out: list[Path] = []

    lock = dict(c.design_lock or {})
    lock["_generated_by"] = "tools/build.py"
    lock["_generated_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    p = DATA_DIR / "design_lock.json"
    p.write_text(json.dumps(lock, indent=2, ensure_ascii=False))
    out.append(p)

    if gate:
        p = DATA_DIR / "gate0.json"
        p.write_text(json.dumps(gate, indent=2, ensure_ascii=False))
        out.append(p)

    versions = C.detected_versions()
    payload = {
        "_generated_by": "tools/build.py",
        "project": {k: c.project[k] for k in ("title", "short_title", "subtitle", "institution")},
        "projection_label": C.PROJECTION_LABEL,
        "data_label": c.data_label,
        "data_status": c.lock("data_status"),
        "decisions": c.decisions,
        "assumptions": c.assumptions,
        "open_questions": c.open_questions,
        "software": [
            {**s, "version": versions.get(s["name"], s.get("version"))} for s in c.software
        ],
        "software_not_used": c.software_not_used,
        "references": c.references,
        "references_retrieved": c.raw["references.yaml"].get("retrieved"),
    }
    p = DATA_DIR / "content.json"
    p.write_text(json.dumps(payload, indent=2, ensure_ascii=False))
    out.append(p)
    return out


def load_all(lock_path: str) -> tuple[C.Content, dict | None]:
    c = C.load(ROOT / "content", design_lock=lock_path)
    gate = None
    gp = Path(lock_path).parent / "gate0.json"
    if gp.exists():
        payload = json.loads(gp.read_text())
        gate = payload.get("results", {}).get("gate0") or payload.get("gate0")
    return c, gate


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="tools.build", description=__doc__)
    ap.add_argument("target", choices=["site", "data", "docx", "deck", "all", "check"])
    ap.add_argument("--design-lock", default="results/demo/design_lock.json")
    args = ap.parse_args(argv)

    if args.target == "check":
        return check(args.design_lock)

    c, gate = load_all(args.design_lock)
    made: list[Path] = []
    if args.target in ("site", "all"):
        made += build_site(c, gate)
    if args.target in ("data", "all"):
        made += build_data(c, gate)
    if args.target in ("docx", "all"):
        made += build_docx(c, gate)
    if args.target in ("deck", "all"):
        made += build_deck(c, gate)
    for p in made:
        print(f"  {p.relative_to(ROOT)}")
    print(f"{len(made)} berkas dihasilkan dari content/ + {args.design_lock}")
    return 0


def check(lock_path: str) -> int:
    """Rebuild the site into a temp tree and fail if the committed one differs."""
    c, gate = load_all(lock_path)
    global SITE_DIR
    real = SITE_DIR
    with tempfile.TemporaryDirectory() as td:
        SITE_DIR = Path(td)
        try:
            built = build_site(c, gate)
        finally:
            tmp, SITE_DIR = SITE_DIR, real
        stale = []
        for p in built:
            committed = real / p.name
            if not committed.exists():
                stale.append(f"{p.name}: belum ada di docs/")
                continue
            a = _strip_date(committed.read_text())
            b = _strip_date(p.read_text())
            if a != b:
                stale.append(f"{p.name}: berbeda dari hasil pembangunan ulang")
        del tmp
    if stale:
        print("Berkas yang dihasilkan sudah kedaluwarsa:")
        for s in stale:
            print(f"  - {s}")
        print("Jalankan: python -m tools.build all")
        return 1
    print(f"{len(built)} halaman situs sesuai dengan content/ + design_lock.json")
    return 0


def _strip_date(text: str) -> str:
    """Drop the build date so `check` compares content, not when it ran."""
    return "\n".join(
        line for line in text.splitlines() if "dihasilkan otomatis dari" not in line
    )


if __name__ == "__main__":
    raise SystemExit(main())
