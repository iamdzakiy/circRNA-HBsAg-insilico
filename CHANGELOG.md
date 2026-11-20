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
- Sekuens nyata HBsAg protein AAB59969.1 (CDS 157..837; genom induk U95551.1, subtipe ayw) di `data/reference/`; terjemahan diverifikasi.
- `analisis_nyata/`: hasil V1–V18, kontrol, pemindaian 45 desain, mutagenesis, primer, GenBank untuk protein AAB59969.1 (genom U95551.1) dengan rubrik v2; `src/circhbsag/rubric2.py`.
- `construct.py`: tata letak CIRC (satu intron utuh) bila `intron_5half` tidak ada. `primers.py`: pasangan primer divergen sejati (komplemen terbalik) dan produk konvergen pada linear dicatat terbuka.
