#!/usr/bin/env python3
"""Assemble the precursor and circle; write FASTA, feature table and GenBank."""
from __future__ import annotations

import csv

from _common import load_parts, out_dir, parser, sidecar

from circhbsag import construct


def main() -> int:
    a = parser("Assemble construct").parse_args()
    parts, meta = load_parts(a)
    asm = construct.build(parts)
    t = out_dir(a, "tables")
    label = "DEMO-PLACEHOLDER" if meta.get("synthetic_parts") else "USER-SEQUENCES"
    fa = t / "construct.fa"
    fa.write_text(
        f">precursor len={len(asm.precursor)} status={label}\n{asm.precursor}\n"
        f">circle len={len(asm.circle)} status={label}\n{asm.circle}\n"
    )
    with (t / "construct_features.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["region", "start_1based", "end_1based", "length_nt"])
        for name, (s, e) in asm.precursor_coords.items():
            w.writerow([name, s, e, e - s + 1])
    try:
        from Bio.Seq import Seq
        from Bio.SeqFeature import FeatureLocation, SeqFeature
        from Bio.SeqRecord import SeqRecord

        rec = SeqRecord(Seq(asm.precursor), id="precursor", name="circHBsAg",
                        description=f"circRNA-HBsAg PIE precursor ({label})",
                        annotations={"molecule_type": "DNA"})
        for name, (s, e) in asm.precursor_coords.items():
            rec.features.append(SeqFeature(FeatureLocation(s - 1, e), type="misc_feature",
                                           qualifiers={"label": [name]}))
        from Bio import SeqIO

        SeqIO.write(rec, t / "precursor.gb", "genbank")
    except Exception as exc:  # pragma: no cover
        print(f"GenBank skipped: {exc}")
    print(f"precursor {len(asm.precursor)} nt, circle {len(asm.circle)} nt -> {t}")
    sidecar(fa, a, {"data_status": label})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
