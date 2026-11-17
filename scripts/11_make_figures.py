#!/usr/bin/env python3
"""Figures: dG vs homology-arm length, and GC profile of the transcribed precursor."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from _common import load_parts, out_dir, parser, sidecar  # noqa: E402
from _gate_view import test_metrics  # noqa: E402

from circhbsag import construct, gate0, motifs, seqio  # noqa: E402
from circhbsag import parts as parts_mod  # noqa: E402

LABEL = "Proyeksi berbasis literatur — bukan data eksperimen"


def main() -> int:
    a = parser("Figures").parse_args()
    fd = out_dir(a, "figures")
    t3 = test_metrics(a.results, "uji3")["metrics"]
    # Scan the UNtrimmed arms so the curve shows why the chosen length was chosen.
    raw, meta = parts_mod.load_parts(a.parts)
    label = LABEL + (" | [DEMO — SEKUENS PLACEHOLDER]" if meta.get("synthetic_parts") else "")
    scan = [r for r in gate0.arm_length_scan(raw["HA5"].seq, raw["HA3"].seq,
                                             lengths=(10, 15, 20, 25, 30, 35, 40),
                                             temperature_c=t3["temperature_c"]) if r.get("available")]
    fig, ax = plt.subplots(figsize=(5.5, 3.6), dpi=200)
    ax.plot([r["length_nt"] for r in scan], [r["binding_kcal"] for r in scan], "o-", color="#1f4e79")
    ax.axhspan(-25, -8, color="#2e7d32", alpha=0.15, label="band target (−8…−25)")
    ax.set_xlabel("Panjang lengan homologi (nt)")
    ax.set_ylabel(f"ΔG dupleks lengan (kcal/mol, {t3['temperature_c']:.0f} °C)")
    ax.legend(frameon=False, fontsize=8)
    fig.text(0.01, 0.01, label, fontsize=6, color="#555")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    f1 = fd / "dG_vs_arm_length.png"
    fig.savefig(f1)
    plt.close(fig)

    parts, _ = load_parts(a)
    asm = construct.build(parts)
    prof = seqio.gc_windows(asm.precursor, window=50, step=10)
    fig, ax = plt.subplots(figsize=(6.5, 3.2), dpi=200)
    ax.plot([s for s, _ in prof], [g for _, g in prof], color="#8e3b00")
    ax.axhspan(motifs.GC_WINDOW_LO, motifs.GC_WINDOW_HI, color="#2e7d32", alpha=0.12)
    ax.set_xlabel("Posisi pada prekursor (nt)")
    ax.set_ylabel("GC jendela 50 nt")
    fig.text(0.01, 0.01, label, fontsize=6, color="#555")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    f2 = fd / "gc_profile.png"
    fig.savefig(f2)
    plt.close(fig)
    print(f1, f2, sep="\n")
    sidecar(f1, a, {"label": label})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
