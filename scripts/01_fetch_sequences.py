#!/usr/bin/env python3
"""Fetch reference sequences listed in data/raw/accessions.tsv and record provenance.

Reads accession numbers YOU provide (never guessed), downloads them from NCBI
E-utilities (public, free), and writes data/raw/provenance.tsv with accession,
access date, SHA-256 and status. If the network or an accession is unavailable the
row is marked UNAVAILABLE; nothing is substituted.
"""
from __future__ import annotations

import csv
import hashlib
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

from _common import ROOT, parser, sidecar

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fetch(acc: str, db: str) -> str:
    q = urllib.parse.urlencode({"db": db, "id": acc, "rettype": "fasta", "retmode": "text"})
    with urllib.request.urlopen(f"{EFETCH}?{q}", timeout=30) as r:  # noqa: S310 - fixed https host
        return r.read().decode()


def main() -> int:
    p = parser("Fetch sequences by accession and write provenance.tsv", orf=False)
    p.add_argument("--accessions", default=str(ROOT / "data/raw/accessions.tsv"))
    p.add_argument("--out", default=str(ROOT / "data/raw"))
    a = p.parse_args()
    acc_file, out = Path(a.accessions), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(acc_file.open(), delimiter="\t"))
    prov = []
    for r in rows:
        acc = (r.get("accession") or "").strip()
        rec = {"part": r["part"], "accession": acc or "[PERLU VERIFIKASI]", "db": r.get("db", "nuccore"),
               "access_date": date.today().isoformat(), "sha256": "", "status": ""}
        if not acc:
            rec["status"] = "NO-ACCESSION (isi data/raw/accessions.tsv)"
        else:
            try:
                text = fetch(acc, rec["db"])
                (out / f"{r['part']}.fa").write_text(text)
                rec["sha256"] = hashlib.sha256(text.encode()).hexdigest()
                rec["status"] = "OK"
            except Exception as exc:  # network blocked, bad accession, ...
                rec["status"] = f"UNAVAILABLE: {type(exc).__name__}: {exc}"
        prov.append(rec)
        print(f"{rec['part']:16s} {rec['status']}")
    tsv = out / "provenance.tsv"
    with tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(prov[0]), delimiter="\t")
        w.writeheader()
        w.writerows(prov)
    sidecar(tsv, a)
    return 0 if all(r["status"] == "OK" for r in prov) else 1


if __name__ == "__main__":
    raise SystemExit(main())
