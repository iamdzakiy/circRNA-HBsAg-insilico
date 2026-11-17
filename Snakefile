# ============================================================================
# Alur kerja reprodusibel Tahap 0 + Tahap 1 — circRNA-HBsAg (PIE, IRES CVB3)
# ============================================================================
# Jalankan:
#   snakemake --cores 4                        # konfigurasi sungguhan (config/parts.yaml)
#   snakemake --cores 4 --config parts=config/parts.example.yaml outdir=results/example
#   snakemake --cores 1 --forcerun all         # paksa ulang semuanya
#   snakemake --dag | dot -Tsvg > docs/dag.svg # diagram alur
#
# Semua langkah memakai seed tetap (config `seed`) sehingga keluaran identik
# antar-jalankan pada versi alat yang sama. Versi alat dicatat di setiap
# berkas laporan JSON.
#
# Lisensi Snakemake: MIT.
# ============================================================================

import sys
from pathlib import Path

configfile: "config/config.yaml"

PARTS = config.get("parts", "config/parts.yaml")
OUT = Path(config.get("outdir", "results"))
SEED = int(config.get("seed", 20260101))
TEMPS = config.get("temperatures", "37,55")
SENS = int(config.get("sensitivity_draws", 200))
CODON_ITER = int(config.get("codon_iterations", 20000))
OPT_TRIALS = int(config.get("optuna_trials", 40))
OPT_CODON_ITER = int(config.get("optuna_codon_iterations", 2000))
BSJ = config.get("bsj_simulation", {})
BLAST_DB = config.get("blast_db") or None
REF_FASTA = config.get("offtarget_reference_fasta") or None

# `src` ditambahkan ke PYTHONPATH agar paket dapat dipakai tanpa instalasi.
PY = f"PYTHONPATH={Path('src').resolve()} python3 -m circhbsag.cli"


rule all:
    """Target lengkap: ORF teroptimasi, keputusan Tahap 0, oligo, BSJ, off-target, ringkasan."""
    input:
        OUT / "orf_optimized.fa",
        OUT / "tahap0.json",
        OUT / "tahap0.md",
        OUT / "tahap0.html",
        OUT / "primers.json",
        OUT / "bsj_simulation.json",
        OUT / "offtarget_report.json",
        OUT / "licenses.json",
        OUT / "SUMMARY.md",


rule licenses:
    """Catat versi dan lisensi semua alat yang dipakai pada jalankan ini."""
    output:
        OUT / "licenses.json",
    shell:
        "{PY} licenses > {output}"


rule codon_optimize:
    """Tahap 1: optimasi kodon ORF HBsAg (CAI, GC, motif, m6A)."""
    input:
        parts=PARTS,
    output:
        fa=OUT / "orf_optimized.fa",
        json=OUT / "codon.json",
    params:
        iters=CODON_ITER,
    shell:
        "{PY} codon --parts {input.parts} --out {OUT} "
        "--iterations {params.iters} --seed {SEED}"


rule tahap0:
    """Tahap 0: pengukuran lipatan/aksesibilitas/motif + keputusan GO/NO-GO.

    Dijalankan pada konstruk yang memuat ORF hasil optimasi kodon. `|| true`
    tidak dipakai: bila keputusan NO-GO, CLI keluar dengan kode 2 dan Snakemake
    menghentikan alur, sehingga keputusan NO-GO tidak pernah terlewat diam-diam.
    Untuk tetap melanjutkan ke desain oligo pada kasus NO-GO, jalankan target
    masing-masing secara eksplisit.
    """
    input:
        parts=PARTS,
        orf=OUT / "orf_optimized.fa",
    output:
        json=OUT / "tahap0.json",
        md=OUT / "tahap0.md",
        html=OUT / "tahap0.html",
    params:
        sens=SENS,
    shell:
        "{PY} tahap0 --parts {input.parts} --orf-fasta {input.orf} --out {OUT} "
        "--temperatures {TEMPS} --sensitivity {params.sens} --seed {SEED}"


rule gate0:
    """Gate 0: sembilan uji + design_lock.json + Gate 0 Report.

    Bukan bagian dari `rule all`: CLI keluar dengan kode 2 pada FAIL-PINDAH-JALUR,
    dan keputusan itu harus menghentikan alur, bukan ditelan. Jalankan eksplisit:
        snakemake --cores 2 gate0
    Prasyarat `cmsearch` (Infernal) dan data/reference/RF00028.cm untuk Uji 2;
    tanpanya Uji 2 dilaporkan UNAVAILABLE, bukan lulus.
    """
    input:
        parts=PARTS,
        orf=OUT / "orf_optimized.fa",
    output:
        gate=OUT / "reports" / "gate0.json",
        lock=OUT / "reports" / "design_lock.json",
    shell:
        "{PY} gate0 --parts {input.parts} --orf-fasta {input.orf} --out {OUT}/reports "
        "--seed {SEED} --cm data/reference/RF00028.cm --temperatures {TEMPS}"


rule primers:
    """Tahap 1: primer divergen, primer BSJ, probe RNase H, primer reverse IVT."""
    input:
        parts=PARTS,
        orf=OUT / "orf_optimized.fa",
    output:
        OUT / "primers.json",
    shell:
        "{PY} primers --parts {input.parts} --orf-fasta {input.orf} "
        "--out {OUT} --seed {SEED}"


rule bsj_simulation:
    """Tahap 1: simulasi bacaan + filter statistik BSJ (posterior Bayes)."""
    input:
        parts=PARTS,
        orf=OUT / "orf_optimized.fa",
    output:
        OUT / "bsj_simulation.json",
    params:
        nc=int(BSJ.get("n_circular", 1500)),
        nl=int(BSJ.get("n_linear", 1500)),
        rl=int(BSJ.get("read_len", 150)),
        anchor=int(BSJ.get("anchor", 20)),
    shell:
        "{PY} bsj-sim --parts {input.parts} --orf-fasta {input.orf} --out {OUT} "
        "--n-circular {params.nc} --n-linear {params.nl} "
        "--read-len {params.rl} --anchor {params.anchor} --seed {SEED}"


rule offtarget:
    """Validasi: homologi off-target ke transkriptom/genom manusia.

    Memakai BLAST+ bila `blast_db` diisi di config dan `blastn` ada di PATH;
    bila tidak, memakai filter k-mer terhadap FASTA lokal; bila keduanya tidak
    tersedia, menulis status SKIPPED yang eksplisit (tidak mengarang hasil).
    """
    input:
        parts=PARTS,
        orf=OUT / "orf_optimized.fa",
    output:
        OUT / "offtarget_report.json",
    params:
        db=f"--db {BLAST_DB}" if BLAST_DB else "",
        ref=f"--reference {REF_FASTA}" if REF_FASTA else "",
    shell:
        "{PY} offtarget --parts {input.parts} --orf-fasta {input.orf} "
        "--out {OUT} {params.db} {params.ref} --seed {SEED}"


rule optimize:
    """Siklus optimasi multi-objektif (Optuna NSGA-II) atas spacer/HA/kodon.

    Bukan bagian dari `rule all` karena memakan waktu; jalankan eksplisit:
        snakemake --cores 4 results/optimization.json
    atau setelah keputusan NO-GO/GO-BERSYARAT pada Tahap 0.
    """
    input:
        parts=PARTS,
    output:
        OUT / "optimization.json",
    params:
        trials=OPT_TRIALS,
        citer=OPT_CODON_ITER,
    shell:
        "{PY} optimize --parts {input.parts} --out {OUT} "
        "--n-trials {params.trials} --codon-iterations {params.citer} --seed {SEED}"


rule summary:
    """Ringkasan satu halaman dari semua berkas JSON hasil."""
    input:
        tahap0=OUT / "tahap0.json",
        codon=OUT / "codon.json",
        primers=OUT / "primers.json",
        bsj=OUT / "bsj_simulation.json",
        offtarget=OUT / "offtarget_report.json",
        licenses=OUT / "licenses.json",
    output:
        OUT / "SUMMARY.md",
    run:
        import json

        def load(p):
            return json.loads(Path(p).read_text())

        t0 = load(input.tahap0)
        cd = load(input.codon)
        pr = load(input.primers)
        bs = load(input.bsj)
        dec = t0["results"]["decision"]
        asm = t0["results"]["measurements"]["assembly"]
        met = cd["results"]["codon"]["cost"]["metrics"]
        sup = bs["results"]["bsj_simulation"]["support"]
        lines = [
            "# Ringkasan jalankan pipeline circRNA-HBsAg",
            "",
            f"**Label data:** {t0['data_label']}",
            "",
            f"- Dihasilkan (UTC): {t0['generated_utc']}",
            f"- Seed: {t0.get('seed')}",
            f"- ViennaRNA: {t0['tool_licenses']['ViennaRNA']['version']}",
            "",
            "## Keputusan Tahap 0",
            "",
            f"**{dec['verdict']}** — skor komposit {dec['composite_score']:.4f} "
            f"(ambang GO {dec['thresholds']['go']}, bersyarat {dec['thresholds']['conditional']})",
            "",
            f"Gerbang keras yang gagal: "
            f"{dec['hard_gates_failed'] if dec['hard_gates_failed'] else 'tidak ada'}",
            "",
            "| Kriteria | Bobot | Subskor |",
            "|---|---:|---:|",
        ]
        for c in dec["criteria"]:
            lines.append(f"| {c['name']} | {c['weight']:.2f} | {c['subscore']:.3f} |")
        lines += [
            "",
            "## Konstruk",
            "",
            f"- Prekursor {asm['precursor_length']} nt, lingkaran {asm['circle_length']} nt, "
            f"BSJ pada {asm['bsj_position']}",
            f"- ORF teroptimasi: CAI {met['cai']}, GC {met['gc']}, "
            f"m6A {met['m6a_sites']} situs, Tipe IIS {met['restriction_sites']}",
            "",
            "## Oligo",
            "",
            f"- Primer divergen: {len(pr['results']['primers']['divergent']['pairs'])} pasangan",
            f"- Primer BSJ: {len(pr['results']['primers']['bsj']['pairs'])} pasangan",
            f"- Probe RNase H: {len(pr['results']['primers']['rnase_h']['probes'])}",
            "",
            "## Filter statistik BSJ (simulasi)",
            "",
            f"- Kandidat {sup['n_candidates']}, lolos {sup['n_passing']}, "
            f"posterior gabungan {sup['joint_posterior_circular']:.6g}",
            f"- {sup['confusion_on_simulated_origin']}",
            "",
            "Laporan lengkap: `tahap0.md`, `tahap0.html`, dan berkas JSON di direktori ini.",
        ]
        Path(output[0]).write_text("\n".join(lines) + "\n")
