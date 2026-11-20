"""Logbook generator: the timeline must be validated and its outputs fresh."""
from __future__ import annotations

import json

import pytest
import yaml

from tools import logbook as L


def test_outputs_are_fresh():
    assert L.main(["--check"]) == 0


def test_timeline_json_has_entries_and_placeholders():
    d = json.loads((L.ROOT / "docs" / "data" / "timeline.json").read_text())
    ents = d["entries"]
    assert len(ents) > 10
    assert any(e.get("status") == "PERLU-DIISI" for e in ents)


@pytest.mark.parametrize("mutate,msg", [
    (lambda d: d["entries"].append(dict(d["entries"][0])), "id ganda"),
    (lambda d: d["entries"][0].update(phase="zz"), "fase tidak dikenal"),
    (lambda d: d["entries"][0].update(image="tidak-ada.png"), "gambar tidak ada"),
])
def test_load_rejects_bad_entries(tmp_path, monkeypatch, mutate, msg):
    data = yaml.safe_load((L.ROOT / "content" / "timeline.yaml").read_text())
    ev = next(e for e in data["entries"] if e.get("status") != "PERLU-DIISI")
    data["entries"].remove(ev)
    data["entries"].insert(0, ev)
    mutate(data)
    (tmp_path / "content").mkdir()
    (tmp_path / "content" / "timeline.yaml").write_text(yaml.safe_dump(data, allow_unicode=True))
    (tmp_path / "docs" / "logbook" / "img").mkdir(parents=True)
    for g in data["gallery"]:
        (tmp_path / "docs" / "logbook" / "img" / g["file"]).write_bytes(b"x")
    for e in data["entries"]:
        if e.get("image") and e["image"] != "tidak-ada.png":
            (tmp_path / "docs" / "logbook" / "img" / e["image"]).write_bytes(b"x")
    monkeypatch.setattr(L, "ROOT", tmp_path)
    with pytest.raises(SystemExit, match=msg):
        L.load()
