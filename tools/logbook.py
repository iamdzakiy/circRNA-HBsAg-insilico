"""Generate docs/logbook.html and docs/data/timeline.json from content/timeline.yaml.

Writes ONLY those two paths (plus nothing else in docs/), so it cannot clobber the
hand-maintained pages. ``--check`` exits 1 when the committed files are stale.

    python -m tools.logbook          # regenerate
    python -m tools.logbook --check  # CI: fail if stale
"""
from __future__ import annotations

import argparse
import html
import json
import sys
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import yaml

ROOT = Path(__file__).resolve().parents[1]
WIB = ZoneInfo("Asia/Jakarta")
MONTHS = ["Januari", "Februari", "Maret", "April", "Mei", "Juni", "Juli", "Agustus",
          "September", "Oktober", "November", "Desember"]
KIND = {"catatan": "catatan penulis", "thread": "catatan kerja proyek", "file": "dokumen proyek",
        "commit": "riwayat git"}
STATUS_CLS = {"BERJALAN": "info"}
PHASE_COLORS = {"p0": "#9aa59d", "p1": "#6b8f7a", "p2": "#15705e", "p3": "#8a6400", "p4": "#7a3b8f"}


def esc(s) -> str:
    return html.escape(str(s), quote=True)


def parse_when(e: dict) -> datetime:
    w = str(e["when"])
    if len(w) == 10:
        return datetime.fromisoformat(w + "T00:00:00+00:00")
    return datetime.fromisoformat(w.replace("Z", "+00:00"))


def fmt_date(dt: datetime, precision: str | None) -> str:
    if precision == "month":
        return f"{MONTHS[dt.month - 1]} {dt.year}"
    d = dt if precision == "date" else dt.astimezone(WIB)
    return f"{d.day} {MONTHS[d.month - 1]} {d.year}"


def fmt_time(dt: datetime, precision: str | None) -> str:
    if precision in ("date", "month"):
        return "tanpa jam"
    return dt.astimezone(WIB).strftime("%H.%M") + " WIB"


def day_key(dt: datetime, precision: str | None) -> str:
    d = dt if precision in ("date", "month") else dt.astimezone(WIB)
    return d.strftime("%Y-%m-%d")


def load() -> dict:
    data = yaml.safe_load((ROOT / "content" / "timeline.yaml").read_text())
    phases = {p["id"] for p in data["phases"]}
    ids = set()
    for e in data["entries"]:
        if e["id"] in ids:
            raise SystemExit(f"id ganda: {e['id']}")
        ids.add(e["id"])
        if e["phase"] not in phases:
            raise SystemExit(f"{e['id']}: fase tidak dikenal {e['phase']}")
        for img in [e.get("image")] if e.get("image") else []:
            if not (ROOT / "docs" / "logbook" / "img" / img).exists():
                raise SystemExit(f"{e['id']}: gambar tidak ada: {img}")
        if not e.get("kind"):
            raise SystemExit(f"{e['id']}: entri bukti harus punya `kind` (sumber)")
    for g in data["gallery"]:
        if not (ROOT / "docs" / "logbook" / "img" / g["file"]).exists():
            raise SystemExit(f"galeri: gambar tidak ada: {g['file']}")
    data["entries"].sort(key=parse_when)
    return data


def scheme_svg(data: dict) -> str:
    ents = [e for e in data["entries"] if not e.get("period")]
    t0 = min(parse_when(e) for e in ents)
    t1 = max(parse_when(e) for e in ents)
    span = (t1 - t0).total_seconds() or 1
    W, H, PAD = 880, 190, 24
    x = lambda dt: PAD + (dt - t0).total_seconds() / span * (W - 2 * PAD - 70) + 70  # noqa: E731
    out = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-labelledby="schT schD">',
           '<title id="schT">Skema fase proyek</title>',
           '<desc id="schD">Garis waktu dari 26 Oktober sampai 20 November 2026 dengan titik '
           'untuk tiap catatan, diwarnai menurut fase. Mei sampai 25 Oktober adalah masa kajian awal '
           'yang diringkas per bulan.</desc>']
    out.append(f'<rect x="{PAD}" y="62" width="58" height="26" rx="4" fill="none" '
               f'stroke="{PHASE_COLORS["p0"]}"/>')
    out.append(f'<text x="{PAD + 29}" y="79" text-anchor="middle" font-size="10" fill="var(--fg-soft)">Mei–25 Okt</text>')
    out.append(f'<text x="{PAD + 29}" y="102" text-anchor="middle" font-size="9" fill="var(--mut)">kajian awal</text>')
    out.append(f'<line x1="{PAD + 70}" y1="75" x2="{W - PAD}" y2="75" stroke="var(--line)" stroke-width="2"/>')
    # day ticks
    seen = set()
    for e in ents:
        dt = parse_when(e)
        key = day_key(dt, e.get("precision"))
        if key in seen:
            continue
        seen.add(key)
    ticks = sorted(seen)
    for k in ticks:
        dt = datetime.fromisoformat(k + "T12:00:00+07:00")
        px = x(dt.astimezone(ZoneInfo("UTC")))
        px = min(max(px, PAD + 70), W - PAD)
        out.append(f'<line x1="{px:.1f}" y1="70" x2="{px:.1f}" y2="80" stroke="var(--mut)"/>')
        lab = dt.strftime("%d/%m")
        out.append(f'<text x="{px:.1f}" y="96" text-anchor="middle" font-size="9" fill="var(--mut)">{lab}</text>')
    lane = {"p1": 40, "p2": 56, "p3": 112, "p4": 128}
    for e in ents:
        dt = parse_when(e)
        px = min(max(x(dt), PAD + 70), W - PAD)
        py = lane.get(e["phase"], 75)
        col = PHASE_COLORS[e["phase"]]
        out.append(f'<a href="#{esc(e["id"])}"><circle cx="{px:.1f}" cy="{py}" r="4.5" fill="{col}">'
                   f'<title>{esc(e["title"])}</title></circle></a>')
    # legend
    lx = PAD
    for p in data["phases"]:
        out.append(f'<circle cx="{lx + 5}" cy="{H - 14}" r="4.5" fill="{PHASE_COLORS[p["id"]]}"/>')
        out.append(f'<text x="{lx + 14}" y="{H - 10}" font-size="10" fill="var(--fg-soft)">{esc(p["name"])}</text>')
        lx += 14 + 6.3 * len(p["name"]) + 16
    out.append("</svg>")
    return "\n".join(out)


def entry_html(e: dict, phase_names: dict) -> str:
    dt = parse_when(e)
    prec = e.get("precision")
    cls = STATUS_CLS.get(e.get("status", ""), "")
    parts = [f'<article class="entry {cls}" id="{esc(e["id"])}" data-phase="{e["phase"]}">',
             f'<header><span class="when">{esc(fmt_date(dt, prec))} <small>{esc(fmt_time(dt, prec))}</small></span>',
             f'<span class="chip ph-{e["phase"]}">{esc(phase_names[e["phase"]])}</span>']
    if e.get("status"):
        parts.append(f'<span class="chip st st-{cls or "x"}">{esc(e["status"].replace("-", " "))}</span>')
    parts.append(f'</header><h3>{esc(e["title"])}</h3>')
    if e.get("summary"):
        parts.append(f'<p>{esc(e["summary"])}</p>')
    if e.get("decisions"):
        parts.append('<p class="lbl">Keputusan</p><ul>' + "".join(f"<li>{esc(d)}</li>" for d in e["decisions"]) + "</ul>")
    if e.get("result"):
        parts.append('<dl class="res">' + "".join(f"<div><dt>{esc(k)}</dt><dd>{esc(v)}</dd></div>"
                                                   for k, v in e["result"].items()) + "</dl>")
    if e.get("outputs"):
        parts.append('<p class="lbl">Keluaran</p><ul class="files">' + "".join(
            f"<li><code>{esc(o)}</code></li>" for o in e["outputs"]) + "</ul>")
    if e.get("note"):
        parts.append(f'<p class="note-in">{esc(e["note"])}</p>')
    if e.get("image"):
        parts.append(f'<figure><img src="logbook/img/{esc(e["image"])}" alt="{esc(e.get("caption", e["title"]))}" '
                     f'loading="lazy"><figcaption>{esc(e.get("caption", ""))}</figcaption></figure>')
    foot = []
    if e.get("kind"):
        foot.append(f"Sumber: {esc(KIND[e['kind']])}")
    if e.get("commit"):
        foot.append(f"commit <code>{esc(e['commit'])}</code>")
    if e.get("link"):
        foot.append(f'<a href="{esc(e["link"])}" rel="noopener">{esc(e["link"])}</a>')
    if foot:
        parts.append('<p class="src">' + " · ".join(foot) + "</p>")
    parts.append("</article>")
    return "\n".join(parts)


CSS = """
.lb-hero{padding:2.2rem 0 1rem}
.lb-hero h1{margin:.2rem 0 .4rem}
.scheme{background:var(--surface);border:1px solid var(--line);border-radius:var(--r-lg);padding:1rem;margin:1.2rem 0}
.scheme svg{width:100%;height:auto;display:block;font-family:var(--sans)}
.filters{display:flex;flex-wrap:wrap;gap:.4rem;margin:1rem 0}
.filters button{border:1px solid var(--line);background:var(--surface);color:var(--fg);border-radius:999px;padding:.3rem .8rem;font:inherit;font-size:.85rem;cursor:pointer}
.filters button[aria-pressed=true]{background:var(--accent);color:var(--accent-ink);border-color:var(--accent)}
.day{margin:1.6rem 0 .4rem;font-size:.8rem;letter-spacing:.07em;text-transform:uppercase;color:var(--mut);border-bottom:1px solid var(--line);padding-bottom:.3rem}
.entry{background:var(--surface);border:1px solid var(--line);border-left:4px solid var(--accent);border-radius:var(--r-md);padding:.9rem 1.1rem;margin:.7rem 0;scroll-margin-top:5rem}
.entry.warn{border-left-color:var(--cond);background:var(--cond-soft)}
.entry header{display:flex;flex-wrap:wrap;gap:.5rem;align-items:center;font-size:.8rem;color:var(--mut)}
.entry .when small{color:var(--mut);margin-left:.3rem}
.entry h3{margin:.4rem 0 .3rem;font-size:1.05rem}
.entry p{margin:.3rem 0;max-width:var(--measure)}
.chip{border-radius:999px;padding:.05rem .6rem;font-size:.72rem;border:1px solid var(--line);background:var(--card);color:var(--fg-soft)}
.chip.st-warn{background:var(--cond-soft);color:var(--cond);border-color:var(--cond)}
.chip.st-info{background:var(--accent-soft);color:var(--accent);border-color:var(--accent)}
.lbl{font-size:.75rem;text-transform:uppercase;letter-spacing:.06em;color:var(--mut);margin:.6rem 0 .1rem!important}
.entry ul{margin:.1rem 0 .3rem 1.1rem;padding:0}
.entry code,.src code{overflow-wrap:anywhere;white-space:normal;word-break:break-word}
.res{display:flex;flex-wrap:wrap;gap:.5rem 1.2rem;margin:.5rem 0}
.res div{background:var(--card);border-radius:var(--r-sm);padding:.3rem .7rem}
.res dt{font-size:.7rem;color:var(--mut)}.res dd{margin:0;font-weight:650}
.note-in{font-size:.88rem;color:var(--cond)}
.entry figure{margin:.7rem 0 0}.entry img{max-width:100%;height:auto;border:1px solid var(--line);border-radius:var(--r-sm)}
figcaption{font-size:.8rem;color:var(--mut);margin-top:.3rem}
.src{font-size:.78rem;color:var(--mut);margin-top:.5rem!important}
.gal{display:grid;grid-template-columns:repeat(auto-fill,minmax(18rem,1fr));gap:1rem}
.gal figure{margin:0;background:var(--surface);border:1px solid var(--line);border-radius:var(--r-md);padding:.7rem}
.gal img{width:100%;height:auto;border-radius:var(--r-sm)}
.gal h3{font-size:.95rem;margin:.1rem 0 .4rem}
.tpl{background:var(--cond-soft);border:1px dashed var(--cond);border-radius:var(--r-md);padding:1rem}
.tpl pre{overflow:auto;font-size:.8rem}
.hidden{display:none!important}
@media (max-width:40rem){.entry{padding:.8rem}}
"""

JS = """
(function(){var b=document.querySelectorAll('.filters button'),e=document.querySelectorAll('.entry'),d=document.querySelectorAll('.day');
function apply(p){b.forEach(function(x){x.setAttribute('aria-pressed',String(x.dataset.p===p))});
e.forEach(function(x){x.classList.toggle('hidden',p!=='all'&&x.dataset.phase!==p)});
d.forEach(function(h){var n=h.nextElementSibling,any=false;while(n&&n.classList.contains('entry')){if(!n.classList.contains('hidden'))any=true;n=n.nextElementSibling}h.classList.toggle('hidden',!any)})}
b.forEach(function(x){x.addEventListener('click',function(){apply(x.dataset.p)})});
var t=document.getElementById('theme');if(t)t.addEventListener('click',function(){var r=document.documentElement,c=r.getAttribute('data-theme')||(matchMedia('(prefers-color-scheme: dark)').matches?'dark':'light'),n=c==='dark'?'light':'dark';r.setAttribute('data-theme',n);try{localStorage.setItem('theme',n)}catch(_){}});})();
"""


def render(data: dict) -> str:
    phases = {p["id"]: p["name"] for p in data["phases"]}
    ents = data["entries"]
    real = [e for e in ents if not e.get("period")]
    first = min(parse_when(e) for e in real)
    last = max(parse_when(e) for e in real)
    days: dict[str, list] = {}
    for e in ents:
        days.setdefault(day_key(parse_when(e), e.get("precision")), []).append(e)
    body = []
    for k in sorted(days):
        dt = datetime.fromisoformat(k + "T00:00:00+00:00")
        body.append(f'<h2 class="day">{dt.day} {MONTHS[dt.month - 1]} {dt.year}</h2>')
        body.extend(entry_html(e, phases) for e in days[k])
    filt = '<button type="button" data-p="all" aria-pressed="true">Semua</button>' + "".join(
        f'<button type="button" data-p="{p["id"]}" aria-pressed="false">{esc(p["name"])}</button>' for p in data["phases"])
    kn = "".join(
        f'<tr><td>{esc(datetime.fromisoformat(r["when"]).strftime("%d/%m/%Y"))}</td><td>{esc(r["what"])}</td>'
        f'<td class="num">{esc(r["value"])}</td><td>{esc(r["verdict"])}</td><td>{esc(r["status"])}</td></tr>'
        for r in data["key_numbers"])
    gal = "".join(
        f'<figure><h3>{esc(g["title"])}</h3><img src="logbook/img/{esc(g["file"])}" alt="{esc(g["title"])}" loading="lazy">'
        f'<figcaption>{esc(g["caption"])}</figcaption></figure>' for g in data["gallery"])
    generated = (f"{len(ents)} catatan, {fmt_date(first, None)} – {fmt_date(last, None)} "
                 f"(didahului kajian awal sejak Mei 2026)")
    return f"""<!doctype html>
<html lang="id">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Logbook proyek circRNA-HBsAg</title>
<meta name="description" content="Linimasa bertanggal proses, keputusan, hasil, dan gambar proyek circRNA-HBsAg.">
<meta name="color-scheme" content="light dark">
<link rel="stylesheet" href="site.css">
<style>{CSS}</style>
<script>(function(){{try{{var t=localStorage.getItem("theme");if(t==="dark"||t==="light")document.documentElement.setAttribute("data-theme",t)}}catch(e){{}}}})();</script>
</head>
<body>
<!-- DIHASILKAN OLEH tools/logbook.py dari content/timeline.yaml. Jangan disunting tangan. -->
<div class="topbar"><div class="wrap topbar-in">
<a class="brand" href="./">circRNA&#8209;HBsAg <span class="thin">in silico</span></a>
<nav class="navlinks" aria-label="Bagian halaman">
<a href="#skema">Skema</a><a href="#linimasa">Linimasa</a><a href="#angka">Angka kunci</a><a href="#galeri">Galeri</a></nav>
<div class="topbar-actions"><a class="btn-link" href="./">Beranda</a><a class="btn-link" href="DESIGN.html">Justifikasi desain</a>
<button class="icon-btn" id="theme" type="button" aria-label="Ganti tema terang atau gelap" title="Ganti tema">◐</button></div>
</div></div>
<main class="wrap" id="isi">
<header class="lb-hero"><span class="eyebrow">Logbook dan dokumentasi proses</span>
<h1>Logbook proyek circRNA&#8209;HBsAg</h1>
<p class="sub">{esc(generated)}. Semua waktu dalam WIB (UTC+7). Angka hasil in silico tidak boleh dibaca
sebagai data eksperimen; setiap entri menyebut sumbernya. Hal yang belum diverifikasi ditandai
<strong>[PERLU VERIFIKASI]</strong>, bukan dikarang.</p></header>

<section id="skema"><h2>Skema fase</h2>
<div class="scheme">{scheme_svg(data)}</div>
<p>Titik = satu catatan; klik untuk melompat ke entrinya. Fase 2 (pipeline) dan fase 3 (dokumen) berjalan paralel pada 16–18 Nov.</p></section>

<section id="linimasa"><h2>Linimasa</h2>
<div class="filters" role="group" aria-label="Saring menurut fase">{filt}</div>
{''.join(body)}
</section>

<section id="angka"><h2>Angka kunci menurut tanggal</h2>
<p>Skor komposit bukan ukuran yang sama sepanjang waktu: rubrik diganti antara 19 dan 20 Nov (lihat catatan 20 Nov), dan baris placeholder tidak sebanding dengan baris sekuens nyata.</p>
<div class="table-scroll"><table><thead><tr><th>Tanggal</th><th>Apa</th><th class="num">Skor</th><th>Keputusan</th><th>Status</th></tr></thead><tbody>{kn}</tbody></table></div></section>

<section id="galeri"><h2>Galeri gambar dan tangkapan layar</h2>
<p>Gambar hasil in silico berasal dari analisis proyek ini; tidak ada gambar jurnal berhak cipta di halaman ini.</p>
<div class="gal">{gal}</div></section>

<section id="cara"><h2>Cara halaman ini dibuat</h2>
<p>Linimasa ditulis dari satu berkas sumber, <code>content/timeline.yaml</code>, lalu dihasilkan oleh <code>tools/logbook.py</code>.
Setiap entri menyebut sumbernya dan generator menolak entri tanpa sumber. Catatan lengkap pipeline ada di
<a href="https://github.com/iamdzakiy/circRNA-HBsAg-insilico/blob/main/docs/LOGBOOK.md">docs/LOGBOOK.md</a>.</p></section>
</main>
<footer class="site"><div class="wrap"><p>Logbook proyek circRNA&#8209;HBsAg · kode MIT · <a href="./">Beranda</a></p></div></footer>
<script>{JS}</script>
</body></html>
"""


def timeline_json(data: dict) -> str:
    out = []
    for e in data["entries"]:
        dt = parse_when(e)
        out.append({**{k: v for k, v in e.items() if k not in ("when", "kind")}, "source": KIND[e["kind"]],
                    "when_utc": str(e["when"]), "date_wib": dt.astimezone(WIB).strftime("%Y-%m-%d")
                    if e.get("precision") not in ("date", "month") else str(e["when"])[:10]})
    return json.dumps({"phases": data["phases"], "entries": out, "key_numbers": data["key_numbers"]},
                      indent=2, ensure_ascii=False)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    data = load()
    targets = {ROOT / "docs" / "logbook.html": render(data),
               ROOT / "docs" / "data" / "timeline.json": timeline_json(data) + "\n"}
    stale = [p for p, t in targets.items() if not p.exists() or p.read_text() != t]
    if a.check:
        for p in stale:
            print(f"kedaluwarsa: {p.relative_to(ROOT)}")
        return 1 if stale else 0
    for p, t in targets.items():
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(t)
        print(p.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
