# Changelog

## Unreleased
- Gate 0 (nine tests) with Infernal/RF00028 positive-controlled Uji 2, design_lock.json.
- Redesign loop (one parameter per iteration) with iteration table; Uji 3 resolved with 20 nt arms, Uji 9 with spacer GC 0.55 (demo data).
- Shared hygiene constants between codon optimiser and Gate 0.
- scripts/01-12 CLIs with parameter/version sidecars; Snakemake `gate0` rule.
- Content single-source-of-truth layer (content/*.yaml) + generators (site, docx, deck).

## Tidak dirilis — Logbook web
- `content/timeline.yaml` + `tools/logbook.py` → `docs/logbook.html`, `docs/data/timeline.json` (linimasa, skema fase, galeri). Juli–25 Okt ditandai `[PERLU DIISI JEK]` (tidak ada sumber). Pemeriksaan kedaluwarsa di CI.
- `scripts/13_export_workbook.py` + rule `workbook`: satu xlsx multi-sheet (≈29 sheet) dan CSV per sheet dari setiap run.
- Sekuens nyata HBsAg D00329.1 (CDS 157..837, AAB59969.1) di `data/reference/`; terjemahan diverifikasi.
