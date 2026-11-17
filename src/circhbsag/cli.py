"""Antarmuka baris perintah: `python -m circhbsag.cli <subcommand>`.

Subcommand
----------
  licenses       cetak daftar alat dan lisensinya
  tahap0         jalankan pengukuran Tahap 0 + skor GO/NO-GO + laporan
  gate0          sembilan uji Gate 0 + design_lock.json (satu sumber angka)
  codon          optimasi kodon ORF dari FASTA protein
  primers        desain primer divergen/BSJ, probe RNase H, primer IVT
  bsj-sim        simulasi bacaan + filter statistik BSJ (posterior Bayes)
  offtarget      pemeriksaan homologi off-target (BLAST+ atau k-mer)
  optimize       optimasi multi-objektif Optuna atas spacer/HA/kodon
  codon-table    hitung tabel penggunaan kodon dari FASTA CDS inang
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _load(cfg_path: str, orf_fasta: str | None = None):
    """Muat bagian konstruk; `orf_fasta` menimpa bagian ORF bila diberikan.

    Alur Snakemake memakai ini untuk menyuntikkan ORF hasil optimasi kodon ke
    dalam konstruk tanpa mengubah parts.yaml.
    """
    from . import parts as parts_mod
    from .construct import Part
    from .seqio import read_fasta

    parts, meta = parts_mod.load_parts(cfg_path)
    arm_nt = (meta.get("design") or {}).get("arm_nt")
    if arm_nt and "HA5" in parts and "HA3" in parts:
        from . import construct

        a5, a3 = construct.trim_arms(parts["HA5"].seq, parts["HA3"].seq, int(arm_nt))
        parts["HA5"] = Part(name="HA5", seq=a5, role=parts["HA5"].role)
        parts["HA3"] = Part(name="HA3", seq=a3, role=parts["HA3"].role)
    if orf_fasta:
        rec = read_fasta(orf_fasta)[0]
        parts["ORF"] = Part(name="ORF", seq=rec.seq, role="ORF teroptimasi kodon")
        meta.setdefault("parts", {})["ORF"] = {
            "file": orf_fasta,
            "fasta_header": rec.seq_id,
            "sha256": rec.sha256,
            "length": len(rec.seq),
            "provenance_status": "DIHASILKAN-PIPELINE (optimasi kodon)",
            "synthetic": False,
        }
    return parts, meta


def cmd_licenses(args) -> int:
    from . import licenses

    print(json.dumps(licenses(), indent=2, ensure_ascii=False))
    return 0


def cmd_tahap0(args) -> int:
    from . import construct, report, score

    parts, meta = _load(args.parts, getattr(args, 'orf_fasta', None))
    asm = construct.build(parts)
    temps = tuple(float(t) for t in args.temperatures.split(","))
    m = score.measure(asm, temperatures=temps)
    dec = score.score(m)
    payload = {
        "measurements": m,
        "decision": dec.to_dict(),
        "permuted_template": construct.permuted_dna_template(asm),
        "guide_sequence_pairs": construct.guide_sequence_pairs(asm),
    }
    if args.sensitivity:
        payload["sensitivity"] = score.sensitivity(m, n=args.sensitivity)
    rep = report.envelope(payload, provenance=meta, seed=args.seed)
    paths = report.write_all(rep, args.out, stem="tahap0")
    print(f"Keputusan: {dec.verdict}  (skor {dec.composite:.4f})")
    for k, v in paths.items():
        print(f"  {k}: {v}")
    return 0 if dec.verdict != "NO-GO" else 2


def cmd_gate0(args) -> int:
    """Sembilan uji Gate 0, laporan, dan design_lock.json.

    design_lock.json sengaja menjadi satu-satunya tempat angka konstruk
    diterbitkan: dokumen TA, dek, dan situs membacanya, tidak menyalinnya.
    """
    from . import construct, report, score
    from . import gate0 as g0

    parts, meta = _load(args.parts, getattr(args, "orf_fasta", None))
    asm = construct.build(parts)
    temps = tuple(float(t) for t in args.temperatures.split(","))
    m = score.measure(asm, temperatures=temps)
    dec = score.score(m)
    gate = g0.run(
        asm,
        m,
        cm_path=args.cm,
        arm_lengths=tuple(int(x) for x in args.arm_lengths.split(",")),
        splice_temperature_c=float(args.splice_temperature),
    )
    lock = g0.design_lock(gate, m, dec.to_dict(), provenance=meta)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rep = report.envelope({"gate0": gate, "decision": dec.to_dict()}, provenance=meta, seed=args.seed)
    (out / "gate0.json").write_text(json.dumps(rep, indent=2, ensure_ascii=False))
    (out / "gate0.md").write_text(
        f"> {rep['data_label']}\n\n" + g0.markdown(gate)
    )
    lock["report_seed"] = args.seed
    lock["circhbsag_version"] = rep.get("circhbsag_version")
    lock["generated_utc"] = rep.get("generated_utc")
    (out / "design_lock.json").write_text(json.dumps(lock, indent=2, ensure_ascii=False))

    print(f"Gate 0: {gate['decision']} — {gate['statement']}")
    for t_ in gate["tests"]:
        print(f"  {t_['id']:32s} {t_['verdict']}")
    print(f"  design_lock: {out / 'design_lock.json'}")
    return 0 if gate["decision"] in ("PASS", "REVISE", "BELUM-LENGKAP") else 2


def cmd_redesign(args) -> int:
    """Loop redesain Gate 0: satu parameter per iterasi, dicatat dalam tabel versi."""
    from . import redesign

    parts, meta = _load(args.parts, getattr(args, "orf_fasta", None))
    steps = [redesign.trim_arms_step(int(x)) for x in args.arm_nt.split(",") if x]
    steps += [redesign.spacer_gc_step(float(x)) for x in args.spacer_gc.split(",") if x]
    its, _, final = redesign.run_loop(
        parts, steps, temperatures=(37.0,), cm_path=args.cm, arm_lengths=(20, 25, 30, 35, 40)
    )
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "gate0_iterations.json").write_text(json.dumps(its, indent=2, ensure_ascii=False))
    (out / "gate0_iterations.md").write_text(redesign.markdown(its))
    print(redesign.markdown(its))
    return 0


def cmd_codon(args) -> int:
    from . import codon, report
    from . import parts as parts_mod

    protein, pmeta = parts_mod.load_protein(args.parts)
    obj = codon.CodonObjective(m6a_mode=args.m6a_mode)
    res = codon.optimize(protein, objective=obj, iterations=args.iterations, seed=args.seed)
    res["protein_provenance"] = pmeta
    res["stats"] = codon.codon_stats(res["orf_dna"])
    rep = report.envelope({"codon": res}, provenance={"parts": {}, "synthetic_parts":
                          ["protein"] if pmeta.get("synthetic") else []}, seed=args.seed)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    report.write_json(rep, out / "codon.json")
    (out / "orf_optimized.fa").write_text(
        f">ORF_HBsAg_codon_optimized len={res['length_nt']} seed={args.seed}\n{res['orf_dna']}\n"
    )
    print(json.dumps(res["cost"]["metrics"], indent=2))
    print(f"ORF -> {out/'orf_optimized.fa'}")
    return 0


def cmd_primers(args) -> int:
    from . import construct, primers, report

    parts, meta = _load(args.parts, getattr(args, 'orf_fasta', None))
    asm = construct.build(parts)
    res = primers.design_all(asm)
    rep = report.envelope({"primers": res}, provenance=meta, seed=args.seed)
    report.write_json(rep, Path(args.out) / "primers.json")
    print(f"Primer divergen: {len(res['divergent']['pairs'])} pasangan; "
          f"BSJ: {len(res['bsj']['pairs'])} pasangan; "
          f"probe RNase H: {len(res['rnase_h']['probes'])}")
    return 0


def cmd_bsj_sim(args) -> int:
    from . import bsj_stats, construct, report

    parts, meta = _load(args.parts, getattr(args, 'orf_fasta', None))
    asm = construct.build(parts)
    res = bsj_stats.run_simulation(
        asm.circle, asm.precursor, n_circular=args.n_circular, n_linear=args.n_linear,
        read_len=args.read_len, anchor=args.anchor, seed=args.seed,
    )
    rep = report.envelope({"bsj_simulation": res}, provenance=meta, seed=args.seed)
    report.write_json(rep, Path(args.out) / "bsj_simulation.json")
    s = res["support"]
    print(f"Kandidat {s['n_candidates']}, lolos {s['n_passing']}, "
          f"posterior gabungan {s['joint_posterior_circular']:.6g}")
    print(f"Confusion: {s['confusion_on_simulated_origin']}")
    return 0


def cmd_offtarget(args) -> int:
    from . import construct, offtarget, report

    parts, meta = _load(args.parts, getattr(args, 'orf_fasta', None))
    asm = construct.build(parts)
    res = offtarget.check_construct(asm, db=args.db, reference_fasta=args.reference, out_dir=args.out)
    rep = report.envelope({"offtarget": res}, provenance=meta, seed=args.seed)
    report.write_json(rep, Path(args.out) / "offtarget_report.json")
    for k, v in res.items():
        print(f"{k}: {v.get('status')} ({v.get('method')})")
    return 0


def cmd_optimize(args) -> int:
    from . import optimize, report
    from . import parts as parts_mod

    parts, meta = _load(args.parts)
    protein, pmeta = parts_mod.load_protein(args.parts)
    base = {k: v for k, v in parts.items() if k != "ORF"}
    study, records = optimize.run_study(
        base, protein, n_trials=args.n_trials, seed=args.seed,
        codon_iterations=args.codon_iterations,
    )
    ranked = optimize.rank_records(records)
    payload = {
        "optimization": {
            "n_trials": args.n_trials,
            "seed": args.seed,
            "n_pareto": len(study.best_trials),
            "pareto_trial_numbers": [t.number for t in study.best_trials],
            "ranked": ranked[: args.top],
            "all_records": records,
        }
    }
    rep = report.envelope(payload, provenance={**meta, "protein": pmeta}, seed=args.seed)
    report.write_json(rep, Path(args.out) / "optimization.json")
    print(f"{len(records)} trial, {len(study.best_trials)} di front Pareto.")
    for r in ranked[:5]:
        print(f"  trial {r['trial']}: skor {r['composite_score']:.4f} {r['verdict']} "
              f"({r['params']})")
    return 0


def cmd_codon_table(args) -> int:
    """Hitung tabel penggunaan kodon dari FASTA CDS (menggantikan tabel [PERKIRAAN])."""
    from .seqio import CODON_TABLE, read_fasta, to_dna

    counts = {c: 0 for c in CODON_TABLE}
    total = 0
    skipped = 0
    for rec in read_fasta(args.fasta):
        s = to_dna(rec.seq)
        if len(s) % 3:
            skipped += 1
            continue
        for i in range(0, len(s), 3):
            c = s[i : i + 3]
            if c in counts:
                counts[c] += 1
                total += 1
    if not total:
        print("Tidak ada kodon terbaca.", file=sys.stderr)
        return 1
    lines = [
        f"# Tabel penggunaan kodon dihitung dari {args.fasta}",
        f"# {total} kodon dari FASTA; {skipped} rekaman dilewati (panjang bukan kelipatan 3)",
        "# STATUS PROVENANS: [HASIL IN SILICO] dihitung dari data yang disediakan pengguna",
        "codon\taa\tfrequency_per_thousand",
    ]
    for c in sorted(counts):
        lines.append(f"{c}\t{CODON_TABLE[c]}\t{1000.0 * counts[c] / total:.3f}")
    Path(args.out).write_text("\n".join(lines) + "\n")
    print(f"Tabel ditulis ke {args.out} ({total} kodon)")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="circhbsag", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp, *, orf: bool = True):
        sp.add_argument("--parts", default="config/parts.yaml", help="berkas parts.yaml")
        sp.add_argument("--out", default="results", help="direktori keluaran")
        sp.add_argument("--seed", type=int, default=20260101)
        if orf:
            sp.add_argument("--orf-fasta", dest="orf_fasta", default=None,
                            help="FASTA ORF teroptimasi yang menimpa bagian ORF pada parts.yaml")
        return sp

    s = common(sub.add_parser("tahap0", help="pengukuran + keputusan GO/NO-GO"))
    s.add_argument("--temperatures", default="37,55")
    s.add_argument("--sensitivity", type=int, default=200,
                   help="jumlah set bobot untuk uji sensitivitas (0 = lewati)")
    s.set_defaults(func=cmd_tahap0)

    s = common(sub.add_parser("gate0", help="sembilan uji Gate 0 + design_lock.json"))
    s.add_argument("--temperatures", default="37,55")
    s.add_argument("--arm-lengths", dest="arm_lengths", default="20,25,30,35,40",
                   help="panjang lengan homologi yang dipindai untuk kurva dG")
    s.add_argument("--splice-temperature", dest="splice_temperature", default="55",
                   help="suhu reaksi splicing untuk uji pasangan lengan homologi")
    s.add_argument("--cm", default=None,
                   help="berkas model kovarians Rfam RF00028 (.cm) untuk cmsearch")
    s.set_defaults(func=cmd_gate0)

    s = common(sub.add_parser("redesign", help="loop redesain Gate 0 (satu parameter per iterasi)"))
    s.add_argument("--arm-nt", dest="arm_nt", default="20", help="panjang lengan, dipisah koma")
    s.add_argument("--spacer-gc", dest="spacer_gc", default="0.45,0.55")
    s.add_argument("--cm", default=None)
    s.set_defaults(func=cmd_redesign)

    s = common(sub.add_parser("codon", help="optimasi kodon ORF"), orf=False)
    s.add_argument("--iterations", type=int, default=20000)
    s.add_argument("--m6a-mode", dest="m6a_mode", choices=["minimize", "ignore"], default="minimize")
    s.set_defaults(func=cmd_codon)

    s = common(sub.add_parser("primers", help="desain primer dan probe"))
    s.set_defaults(func=cmd_primers)

    s = common(sub.add_parser("bsj-sim", help="simulasi bacaan + posterior BSJ"))
    s.add_argument("--n-circular", dest="n_circular", type=int, default=1500)
    s.add_argument("--n-linear", dest="n_linear", type=int, default=1500)
    s.add_argument("--read-len", dest="read_len", type=int, default=150)
    s.add_argument("--anchor", type=int, default=20)
    s.set_defaults(func=cmd_bsj_sim)

    s = common(sub.add_parser("offtarget", help="homologi off-target"))
    s.add_argument("--db", default=None, help="nama basis data BLAST")
    s.add_argument("--reference", default=None, help="FASTA referensi lokal untuk filter k-mer")
    s.set_defaults(func=cmd_offtarget)

    s = common(sub.add_parser("optimize", help="optimasi multi-objektif Optuna"))
    s.add_argument("--n-trials", dest="n_trials", type=int, default=40)
    s.add_argument("--codon-iterations", dest="codon_iterations", type=int, default=2000)
    s.add_argument("--top", type=int, default=10)
    s.set_defaults(func=cmd_optimize)

    s = sub.add_parser("licenses", help="daftar alat dan lisensi")
    s.set_defaults(func=cmd_licenses)

    s = sub.add_parser("codon-table", help="hitung tabel kodon dari FASTA CDS")
    s.add_argument("--fasta", required=True)
    s.add_argument("--out", default="data/reference/codon_usage_from_cds.tsv")
    s.set_defaults(func=cmd_codon_table)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
