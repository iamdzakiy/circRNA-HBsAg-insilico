"""Tes antarmuka Streamlit dengan AppTest (tanpa peramban).

Tes ini menjalankan skrip aplikasi sungguhan, menekan tombol, dan memeriksa
hasil yang muncul — bukan hanya memeriksa bahwa server menyala.
"""
from __future__ import annotations

from pathlib import Path

import pytest

APP = Path(__file__).resolve().parents[1] / "app" / "streamlit_app.py"
pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402

TIMEOUT = 600


@pytest.fixture(scope="module")
def app():
    at = AppTest.from_file(str(APP), default_timeout=TIMEOUT)
    at.run()
    return at


def test_app_starts_without_exception(app):
    assert not app.exception
    assert any("Tahap 0" in t.value for t in app.title)


def test_initial_state_asks_for_input(app):
    # Sebelum dijalankan, aplikasi memberi petunjuk, bukan hasil.
    assert any("Jalankan Tahap 0" in i.value for i in app.info)


def test_load_example_warns_about_synthetic_data():
    at = AppTest.from_file(str(APP), default_timeout=TIMEOUT)
    at.run()
    at.button[0].click().run()          # "Muat data contoh (sintetis)"
    assert not at.exception
    assert at.session_state["synthetic"] is True
    joined = " ".join(w.value for w in at.warning)
    assert "SINTETIS" in joined
    assert "UJI PERANGKAT LUNAK" in joined
    # Semua bagian wajib terisi.
    for key in ("intron_3half", "E2", "IRES", "ORF", "E1", "intron_5half"):
        assert at.session_state[f"seq_{key}"].strip()


def test_run_produces_verdict_and_labels_results_as_software_test():
    at = AppTest.from_file(str(APP), default_timeout=TIMEOUT)
    at.run()
    at.button[0].click().run()                    # muat contoh
    # Kurangi beban: satu suhu, tanpa primer/BSJ/sensitivitas.
    at.sidebar.multiselect[0].set_value([37.0]).run()
    for cb in at.sidebar.checkbox:
        cb.set_value(False)
    # Slider terakhir di sidebar adalah "Jumlah undian uji sensitivitas".
    assert at.sidebar.slider[-1].label == "Jumlah undian uji sensitivitas"
    at.sidebar.slider[-1].set_value(0).run()
    at.button[2].click().run()                    # "Jalankan Tahap 0"
    assert not at.exception

    res = at.session_state["results"]
    assert res is not None
    rep = res["report"]
    assert "[UJI PERANGKAT LUNAK]" in rep["data_label"]
    dec = rep["results"]["decision"]
    assert dec["verdict"] in {"GO", "GO-BERSYARAT", "NO-GO"}
    assert 0.0 <= dec["composite_score"] <= 1.0
    # Keputusan tampil di halaman.
    assert any(dec["verdict"] in md.value for md in at.markdown)
    # Provenans tiap bagian tercatat sebagai sintetis.
    prov = rep["provenance"]["parts"]
    assert prov and all(v["provenance_status"] == "SINTETIS-UJI" for v in prov.values())


def test_missing_required_part_gives_error():
    at = AppTest.from_file(str(APP), default_timeout=TIMEOUT)
    at.run()
    at.session_state["seq_ORF"] = "ATGGCCTGA"      # hanya ORF
    at.run()
    at.button[2].click().run()
    assert not at.exception
    assert at.error, "aplikasi harus menolak masukan yang tidak lengkap"
    assert "belum diisi" in " ".join(e.value for e in at.error)


def test_bad_orf_length_is_rejected():
    at = AppTest.from_file(str(APP), default_timeout=TIMEOUT)
    at.run()
    at.button[0].click().run()
    at.session_state["seq_ORF"] = "ATGGCCTG"       # bukan kelipatan 3
    at.run()
    at.button[2].click().run()
    assert not at.exception
    assert any("kelipatan 3" in e.value for e in at.error)


def test_invalid_characters_are_rejected():
    at = AppTest.from_file(str(APP), default_timeout=TIMEOUT)
    at.run()
    at.button[0].click().run()
    at.session_state["seq_IRES"] = "ACGTNNNNACGT"
    at.run()
    at.button[2].click().run()
    assert not at.exception
    assert any("tidak valid" in e.value or "non-DNA" in e.value for e in at.error)
