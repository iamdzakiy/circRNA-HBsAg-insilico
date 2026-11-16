"""Antarmuka web Tahap 0: masukkan sekuens, jalankan skor GO/NO-GO, unduh laporan.

Jalankan:
    streamlit run app/streamlit_app.py

Lisensi Streamlit: Apache-2.0.

Catatan provenans: apa pun yang Anda tempel di sini diperlakukan sebagai
masukan yang Anda sediakan sendiri. Bila Anda memakai tombol "Muat data contoh",
seluruh keluaran dilabeli [UJI PERANGKAT LUNAK] karena data contoh adalah
sekuens acak, bukan sekuens biologis.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import streamlit as st

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from circhbsag import (  # noqa: E402
    access,
    bsj_stats,
    codon,
    construct,
    fold,
    motifs,
    primers,
    report,
    score,
    seqio,
)

EXAMPLE_DIR = ROOT / "data" / "example"

PART_FIELDS = [
    ("intron_3half", "Separuh intron 3' (PIE)", True, "intron_3half.fa"),
    ("E2", "Eksón E2 (sisi 3' sambungan)", True, "E2.fa"),
    ("spacer5", "Spacer antara E2 dan IRES", False, None),
    ("IRES", "IRES (satu saja, mis. CVB3)", True, "ires_demo.fa"),
    ("ORF", "ORF HBsAg (DNA, kelipatan 3)", True, None),
    ("spacer3", "Spacer antara ORF dan E1", False, None),
    ("E1", "Eksón E1 (sisi 5' sambungan)", True, "E1.fa"),
    ("intron_5half", "Separuh intron 5' (PIE)", True, "intron_5half.fa"),
    ("HA5", "Lengan homologi 5'", False, "HA5.fa"),
    ("HA3", "Lengan homologi 3'", False, "HA3.fa"),
]

st.set_page_config(page_title="Tahap 0 circRNA-HBsAg", page_icon="🧬", layout="wide")


# --------------------------------------------------------------------------
# Status sesi
# --------------------------------------------------------------------------

def _init_state() -> None:
    for key, _label, _req, _fname in PART_FIELDS:
        st.session_state.setdefault(f"seq_{key}", "")
    st.session_state.setdefault("synthetic", False)
    st.session_state.setdefault("results", None)


_init_state()


def load_example() -> None:
    """Isi formulir dengan data contoh SINTETIS dan tandai sesi sebagai uji."""
    for key, _label, _req, fname in PART_FIELDS:
        if fname and (EXAMPLE_DIR / fname).exists():
            st.session_state[f"seq_{key}"] = (EXAMPLE_DIR / fname).read_text()
    st.session_state["seq_spacer5"] = "ACTAACTACAATCAACTACAATCAACTAACAATC"
    st.session_state["seq_spacer3"] = "AAACAAAGAAACAAAGAAACAAAGAAACAAAG"
    # ORF demo dari protein acak 150 aa pada data contoh.
    prot_file = EXAMPLE_DIR / "demo_antigen_protein.fa"
    if prot_file.exists():
        protein = "".join(prot_file.read_text().splitlines()[1:]).strip()
        st.session_state["seq_ORF"] = codon.optimize(protein, iterations=1500, seed=1)["orf_dna"]
    st.session_state["synthetic"] = True
    st.session_state["results"] = None


def clear_all() -> None:
    for key, _label, _req, _fname in PART_FIELDS:
        st.session_state[f"seq_{key}"] = ""
    st.session_state["synthetic"] = False
    st.session_state["results"] = None


# --------------------------------------------------------------------------
# Sidebar
# --------------------------------------------------------------------------

with st.sidebar:
    st.header("Pengaturan")
    temps_choice = st.multiselect(
        "Suhu lipatan (°C)",
        options=[25.0, 30.0, 37.0, 42.0, 50.0, 55.0, 60.0],
        default=[37.0, 55.0],
        help="37 °C = fisiologis. 55 °C = kondisi splicing PIE in vitro. "
             "Setiap kriteria memakai nilai TERBURUK dari suhu yang dipilih.",
    )
    seed = st.number_input("Seed", value=20260101, step=1,
                           help="Seed tetap membuat hasil dapat diulang.")
    st.divider()
    st.subheader("Ambang keputusan")
    th_go = st.slider("Ambang GO", 0.5, 0.95, float(score.DEFAULT_THRESHOLDS["go"]), 0.01)
    th_cond = st.slider("Ambang GO-BERSYARAT", 0.3, 0.9,
                        float(score.DEFAULT_THRESHOLDS["conditional"]), 0.01)
    st.caption("Kedua ambang ini [PERKIRAAN], ditetapkan pipeline ini dan bukan "
               "dari pustaka. Laporan selalu menyertakan uji sensitivitas bobot.")
    st.divider()
    st.subheader("Bobot kriteria")
    weights: dict[str, float] = {}
    for name, default in score.DEFAULT_WEIGHTS.items():
        weights[name] = st.slider(name, 0.0, 0.5, float(default), 0.01)
    total_w = sum(weights.values())
    st.caption(f"Jumlah bobot = {total_w:.2f} (skor dinormalkan, jadi tidak harus 1,00)")
    st.divider()
    st.subheader("Lebih lanjut")
    run_primers = st.checkbox("Rancang primer dan probe RNase H", value=True)
    run_bsj = st.checkbox("Simulasi filter statistik BSJ", value=False,
                          help="Menambah ~10-30 detik.")
    n_sens = st.slider("Jumlah undian uji sensitivitas", 0, 400, 100, 20)


# --------------------------------------------------------------------------
# Header dan masukan
# --------------------------------------------------------------------------

st.title("Tahap 0 — kelayakan sirkularisasi circRNA-HBsAg")
st.markdown(
    "Masukkan sekuens setiap bagian konstruk **permuted intron–exon (PIE)**, "
    "lalu jalankan keputusan GO / NO-GO. Semua angka dihitung saat ini juga oleh "
    f"ViennaRNA `{fold.VIENNA_VERSION}` dan primer3 — tidak ada nilai yang dihardcode."
)
st.code(
    "prekursor:  HA5 | intron 3' | E2 | spacer5 | IRES | ORF | spacer3 | E1 | intron 5' | HA3\n"
    "lingkaran:              E2 | spacer5 | IRES | ORF | spacer3 | E1\n"
    "                        ^-------------------------------------^ sambungan balik (BSJ)",
    language=None,
)

c1, c2, _ = st.columns([1, 1, 4])
c1.button("Muat data contoh (sintetis)", on_click=load_example, use_container_width=True)
c2.button("Kosongkan", on_click=clear_all, use_container_width=True)

if st.session_state["synthetic"]:
    st.warning(
        "**Sesi ini memakai data contoh SINTETIS.** Sekuens contoh dibangkitkan "
        "acak dengan seed tetap dan bukan sekuens biologis: tidak ada nukleotida "
        "dari CVB3, dari intron Anabaena, maupun dari HBV. Seluruh keluaran "
        "dilabeli [UJI PERANGKAT LUNAK] dan tidak boleh dikutip sebagai hasil ilmiah.",
        icon="⚠️",
    )

st.subheader("Bagian konstruk")
st.caption("Tempel FASTA atau sekuens mentah. Spasi, nomor baris, dan baris header FASTA diabaikan. "
           "U otomatis dibaca sebagai T. Bagian bertanda * wajib.")

cols = st.columns(2)
for i, (key, label, required, _fname) in enumerate(PART_FIELDS):
    with cols[i % 2]:
        st.text_area(
            f"{label}{' *' if required else ''}",
            key=f"seq_{key}",
            height=90,
            placeholder="ACGT… atau >header\nACGT…",
        )

st.divider()
run = st.button("Jalankan Tahap 0", type="primary", use_container_width=True)


# --------------------------------------------------------------------------
# Jalankan
# --------------------------------------------------------------------------

def _parse(text: str) -> str:
    """Terima FASTA atau sekuens mentah; kembalikan DNA bersih ('' bila kosong)."""
    lines = [ln for ln in text.splitlines() if not ln.strip().startswith(">")]
    joined = "".join(lines).strip()
    return seqio.clean(joined) if joined else ""


if run:
    try:
        parts: dict[str, construct.Part] = {}
        missing: list[str] = []
        for key, label, required, _fname in PART_FIELDS:
            s = _parse(st.session_state[f"seq_{key}"])
            if not s:
                if required:
                    missing.append(label)
                continue
            features = {}
            if key in ("intron_5half", "intron_3half") and len(s) >= 6:
                # Anotasi IGS bawaan: 6 nt pertama (separuh 5') atau terakhir (separuh 3').
                features = {"IGS": (0, 6)} if key == "intron_5half" else {"IGS": (len(s) - 6, len(s))}
            parts[key] = construct.Part(name=key, seq=s, features=features)
        if missing:
            st.error("Bagian wajib belum diisi: " + ", ".join(missing))
            st.stop()

        parts.setdefault("T7_promoter", construct.Part("T7_promoter", motifs.T7_PROMOTER))
        parts.setdefault("gg_leader", construct.Part("gg_leader", "GGG"))

        orf = parts["ORF"].seq
        if len(orf) % 3:
            st.error(f"Panjang ORF {len(orf)} nt bukan kelipatan 3.")
            st.stop()
        if not orf.startswith("ATG"):
            st.warning("ORF tidak dimulai dengan ATG — gerbang keras akan gagal.", icon="⚠️")

        temps = tuple(sorted(temps_choice)) or (37.0,)
        with st.spinner(f"Melipat dan mengukur pada {len(temps)} suhu…"):
            asm = construct.build(parts)
            m = score.measure(asm, temperatures=temps)
            dec = score.score(m, weights=weights,
                              thresholds={"go": th_go, "conditional": th_cond})
            payload = {
                "measurements": m,
                "decision": dec.to_dict(),
                "permuted_template": construct.permuted_dna_template(asm),
                "guide_sequence_pairs": construct.guide_sequence_pairs(asm),
            }
            if n_sens:
                payload["sensitivity"] = score.sensitivity(m, n=int(n_sens))
            if run_primers:
                payload["primers"] = primers.design_all(asm)
            if run_bsj:
                payload["bsj_simulation"] = bsj_stats.run_simulation(
                    asm.circle, asm.precursor, n_circular=800, n_linear=800, seed=int(seed)
                )
            prov = {
                "parts": {
                    k: {
                        "length": len(v.seq),
                        "sha256": seqio.Record(k, v.seq).sha256,
                        "provenance_status": "SINTETIS-UJI" if st.session_state["synthetic"]
                        else "DISEDIAKAN-PENGGUNA",
                        "source": "data contoh sintetis" if st.session_state["synthetic"]
                        else "ditempel di antarmuka web",
                    }
                    for k, v in asm.parts.items()
                },
                "synthetic_parts": list(asm.parts) if st.session_state["synthetic"] else [],
                "provenance_policy": "web",
            }
            st.session_state["results"] = {
                "asm": asm,
                "report": report.envelope(payload, provenance=prov, seed=int(seed)),
            }
    except seqio.SequenceError as exc:
        st.error(f"Sekuens tidak valid: {exc}")
        st.stop()
    except Exception as exc:  # pragma: no cover - jaring pengaman antarmuka
        st.exception(exc)
        st.stop()


# --------------------------------------------------------------------------
# Hasil
# --------------------------------------------------------------------------

res = st.session_state.get("results")
if not res:
    st.info("Isi bagian konstruk lalu tekan **Jalankan Tahap 0**. "
            "Tidak punya sekuens? Tekan **Muat data contoh (sintetis)** untuk mencoba antarmuka.")
    st.stop()

asm = res["asm"]
rep = res["report"]
r = rep["results"]
dec = r["decision"]
m = r["measurements"]

st.divider()
st.header("Hasil")
st.caption(rep["data_label"])

colour = {"GO": "green", "GO-BERSYARAT": "orange", "NO-GO": "red"}[dec["verdict"]]
st.markdown(f"## :{colour}[{dec['verdict']}]  ·  skor komposit {dec['composite_score']:.3f}")

k1, k2, k3, k4 = st.columns(4)
k1.metric("Panjang lingkaran", f"{m['assembly']['circle_length']} nt")
k2.metric("Panjang prekursor", f"{m['assembly']['precursor_length']} nt")
t0 = f"{min(m['temperatures_c']):g}C"
k3.metric(f"MFE lingkaran @{t0}", f"{m['folds'][f'circle_{t0}']['mfe_kcal']:.1f} kcal/mol")
k4.metric("dsRNA terpanjang", f"{max(v['longest_bp'] for v in m['dsrna'].values())} bp")

failed = dec["hard_gates_failed"]
if failed:
    st.error("**Gerbang keras gagal:** " + ", ".join(f"`{f}`" for f in failed)
             + ". Satu kegagalan gerbang keras memaksa NO-GO terlepas dari skor.")

tabs = st.tabs([
    "Keputusan", "Lipatan & sambungan", "Aksesibilitas", "IRES & dsRNA",
    "Motif & m6A", "Primer & probe", "Unduh",
])

# --- Keputusan
with tabs[0]:
    st.subheader("Gerbang keras")
    st.dataframe(
        [
            {
                "gerbang": k,
                "lulus": "YA" if v["passed"] else "TIDAK",
                "nilai": str(v["value"]),
                "batas": str(v["limit"]),
                "alasan": v["why"],
            }
            for k, v in dec["hard_gates"].items()
        ],
        use_container_width=True, hide_index=True,
    )
    st.subheader("Kriteria berbobot")
    for c in dec["criteria"]:
        st.markdown(f"**`{c['name']}`** — subskor {c['subscore']:.3f} "
                    f"× bobot {c['weight']:.2f} = {c['weighted']:.3f}")
        st.progress(min(1.0, max(0.0, c["subscore"])))
        st.caption(c["rationale"])
        with st.expander("Nilai mentah"):
            st.json(c["raw"], expanded=False)
    if "sensitivity" in r:
        s = r["sensitivity"]
        st.subheader("Sensitivitas terhadap bobot")
        st.markdown(
            f"{s['n']} set bobot dengan jitter ±{int(s['jitter_fraction'] * 100)}%: "
            f"p05 **{s['composite_p05']}**, median **{s['composite_median']}**, "
            f"p95 **{s['composite_p95']}**."
        )
        st.json(s["verdict_counts"])
        st.caption("Bila distribusi keputusan terbelah, keputusan bergantung pada "
                   "pilihan bobot dan harus dinyatakan demikian.")

# --- Lipatan
with tabs[1]:
    st.subheader("Lipatan RNA")
    st.dataframe(
        [
            {
                "molekul": k, "suhu (°C)": v["temperature_c"], "sirkular": v["circular"],
                "MFE (kcal/mol)": v["mfe_kcal"], "MFE/nt": v["mfe_per_nt"],
                "energi ensemble": v["ensemble_free_energy_kcal"],
                "keragaman ensemble": v["ensemble_diversity"],
                "freq. MFE": v["frequency_mfe_in_ensemble"],
            }
            for k, v in m["folds"].items()
        ],
        use_container_width=True, hide_index=True,
    )
    st.subheader("Sambungan balik (BSJ)")
    for k, v in m["junction"].items():
        st.markdown(f"**{k}** — P(tak-berpasangan) rerata **{v['mean_unpaired_probability']:.4f}** "
                    f"pada jendela ±{v['window_span_each_side']} nt")
        st.code(v["structure_window"], language=None)
    st.subheader("Struktur lingkaran (titik-kurung)")
    tkey = st.selectbox("Pilih lipatan", list(m["folds"]), key="struct_pick")
    st.text_area("Struktur", m["folds"][tkey]["mfe_structure"], height=160)
    st.caption("Posisi 1 lingkaran = basa pertama E2; BSJ ada di antara basa terakhir "
               f"({m['assembly']['bsj_position']}) dan basa 1.")

# --- Aksesibilitas
with tabs[2]:
    st.subheader("Aksesibilitas RNAplfold pada prekursor")
    st.caption("Nilai = P(u nukleotida beruntun tak-berpasangan). Situs splice dan "
               "IGS (P1/P10) yang tertutup struktur menandakan risiko sirkularisasi.")
    rows = []
    for t, regions in m["accessibility"].items():
        for name, v in regions.items():
            if v.get("mean") is None:
                continue
            rows.append({"suhu": t, "wilayah": name, "awal": v["start"], "akhir": v["end"],
                         "rerata": v["mean"], "min": v["min"], "maks": v.get("max")})
    st.dataframe(rows, use_container_width=True, hide_index=True)
    st.subheader("Aksesibilitas start codon")
    st.json(m["start_codon"])

# --- IRES & dsRNA
with tabs[3]:
    st.subheader("Integritas struktur IRES")
    if m["ires_integrity"]:
        st.dataframe(
            [{"suhu": k, "identitas status lipatan": v["identity"],
              "MFE IRES sendiri": v["ires_alone_mfe_kcal"], "panjang IRES": v["ires_length"]}
             for k, v in m["ires_integrity"].items()],
            use_container_width=True, hide_index=True,
        )
        st.caption("Identitas = fraksi posisi IRES yang mempertahankan status "
                   "berpasangan/tak-berpasangan seperti saat IRES dilipat sendiri. "
                   "Ini proksi struktural, BUKAN prediksi fungsi translasi: model "
                   "energi ViennaRNA tidak memodelkan pseudoknot.")
    st.subheader("Beban dsRNA")
    st.dataframe(
        [{"lipatan": k, "heliks terpanjang (bp)": v["longest_bp"],
          "jumlah heliks": v["n_helices"], "pelanggaran ≥ ambang": len(v.get("violations", []))}
         for k, v in m["dsrna"].items()],
        use_container_width=True, hide_index=True,
    )
    worst = max(m["dsrna"], key=lambda k: m["dsrna"][k]["longest_bp"])
    st.json(m["dsrna"][worst].get("top_helices", [])[:5])

# --- Motif
with tabs[4]:
    for which in ("precursor", "circle"):
        mm = m["motifs"][which]
        st.subheader(f"{which} — {mm['length']} nt")
        a, b, c, d = st.columns(4)
        a.metric("Situs Tipe IIS", sum(mm["restriction"].values()))
        b.metric("Terminator T7", len(mm["t7_terminator"]))
        c.metric("Situs m6A RRACH", mm["m6a_sites"])
        d.metric("m6A per 100 nt", mm["m6a_density_per_100nt"])
        st.write("Situs restriksi per enzim:", mm["restriction"])
        st.write("Motif yang dihindari:", mm["avoid_motifs"])
        if mm["gc_violations"]:
            st.write(f"Jendela GC di luar 30–70%: {len(mm['gc_violations'])}")
            st.json(mm["gc_violations"][:10])
        st.divider()
    st.caption("Keputusan m6A: bawaan pipeline MENURUNKAN jumlah RRACH. Bila Anda "
               "memang berencana menyisipkan m6A ~1% sebagai fitur, bobot kriteria "
               "`motif_burden` harus diubah, bukan diminimalkan.")

# --- Primer
with tabs[5]:
    if "primers" not in r:
        st.info("Aktifkan **Rancang primer dan probe RNase H** di sidebar, lalu jalankan ulang.")
    else:
        p = r["primers"]
        st.subheader("Primer divergen")
        st.caption("Dirancang pada titik antipodal BSJ: pada templat LINEAR kedua primer "
                   "menjauh sehingga tidak ada produk; hanya lingkaran memberi amplikon.")
        st.dataframe(p["divergent"]["pairs"], use_container_width=True, hide_index=True)
        st.subheader("Primer menumpang BSJ")
        st.dataframe(p["bsj"]["junction_spanning_primers"], use_container_width=True,
                     hide_index=True)
        st.caption(p["bsj"]["note"])
        st.subheader("Probe RNase H")
        st.dataframe(p["rnase_h"]["probes"], use_container_width=True, hide_index=True)
        st.caption(p["rnase_h"]["interpretation"])
        st.subheader("Primer reverse untuk templat IVT")
        st.json(p["ivt_reverse"])
    if "bsj_simulation" in r:
        st.subheader("Filter statistik BSJ (simulasi bacaan)")
        st.json(r["bsj_simulation"]["support"])

# --- Unduh
with tabs[6]:
    st.subheader("Unduh laporan")
    st.caption("Setiap berkas memuat blok provenans, versi alat, dan seed, sehingga "
               "labelnya ikut terbawa ke dokumen Anda.")
    html_bytes = report.html(rep).encode()
    json_bytes = json.dumps(rep, indent=2, default=str, ensure_ascii=False).encode()
    md_bytes = report.markdown(rep).encode()
    d1, d2, d3 = st.columns(3)
    d1.download_button("Laporan HTML", html_bytes, "tahap0.html", "text/html",
                       use_container_width=True)
    d2.download_button("Laporan Markdown", md_bytes, "tahap0.md", "text/markdown",
                       use_container_width=True)
    d3.download_button("Data JSON lengkap", json_bytes, "tahap0.json", "application/json",
                       use_container_width=True)
    st.divider()
    st.subheader("Sekuens hasil rakitan")
    e1, e2 = st.columns(2)
    e1.download_button("Lingkaran (FASTA)",
                       f">circRNA len={len(asm.circle)}\n{asm.circle}\n".encode(),
                       "circle.fa", "text/plain", use_container_width=True)
    e2.download_button("Prekursor / templat (FASTA)",
                       f">precursor len={len(asm.precursor)}\n{asm.precursor}\n".encode(),
                       "precursor.fa", "text/plain", use_container_width=True)
    with st.expander("Koordinat bagian pada lingkaran (1-based, inklusif)"):
        st.json(m["assembly"]["circle_coords"])
    with st.expander("Templat DNA terpermutasi — pemeriksaan T7"):
        st.json(r["permuted_template"])

st.divider()
st.caption(
    "Semua nilai di halaman ini adalah **prediksi** termodinamika dan heuristik "
    "sekuens, bukan pengukuran laboratorium, dan tidak menggantikan uji RNase R, "
    "RT-PCR divergen dengan sekuensing Sanger, atau uji RNase H. "
    f"ViennaRNA {fold.VIENNA_VERSION} · primer3 ({primers.PRIMER3_LICENSE}) · "
    f"circhbsag {rep['circhbsag_version']} (MIT). "
    f"Heliks dsRNA diuraikan oleh `{access.__name__}`."
)
