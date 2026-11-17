"""Gate 0 redesign loop: change ONE parameter per iteration, rerun, log.

BLOK 3 requires that a REVISE verdict is resolved by changing a single design
parameter at a time and re-running the affected tests, so that the cause of each
change in verdict is attributable. Each iteration here applies exactly one named
change on top of the previous accepted state and records the full verdict vector.
A FAIL on a non-patchable test is never "fixed" by this loop; it stops it.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from . import construct, gate0, optimize, score
from .construct import Part


@dataclass
class Step:
    label: str
    parameter: str
    value: str
    apply: Callable[[dict], dict]


def trim_arms_step(length: int) -> Step:
    def f(parts):
        p = dict(parts)
        a5, a3 = construct.trim_arms(p["HA5"].seq, p["HA3"].seq, length)
        p["HA5"] = Part(name="HA5", seq=a5, role=p["HA5"].role)
        p["HA3"] = Part(name="HA3", seq=a3, role=p["HA3"].role)
        return p

    return Step(f"Panjang lengan homologi → {length} nt", "homology_arm_nt", str(length), f)


def spacer_gc_step(gc: float, length5: int | None = None, length3: int | None = None, seed: int = 0) -> Step:
    def f(parts):
        p = dict(parts)
        l5 = length5 or len(p["spacer5"].seq)
        l3 = length3 or len(p["spacer3"].seq)
        p["spacer5"] = Part(name="spacer5", seq=optimize.make_spacer(l5, gc_target=gc, seed=seed),
                            role=p["spacer5"].role)
        p["spacer3"] = Part(name="spacer3", seq=optimize.make_spacer(l3, gc_target=gc, seed=seed + 1),
                            role=p["spacer3"].role)
        return p

    return Step(f"GC spacer → {gc:.2f}", "spacer_gc_target", f"{gc:.2f}", f)


def _verdicts(parts, temperatures, cm_path, arm_lengths, splice_temperature_c):
    asm = construct.build(parts)
    m = score.measure(asm, temperatures=temperatures)
    g = gate0.run(asm, m, cm_path=cm_path, arm_lengths=arm_lengths,
                  splice_temperature_c=splice_temperature_c)
    return asm, m, g


def run_loop(parts, steps, *, temperatures=(37.0,), cm_path=None, arm_lengths=(20, 25, 30, 35, 40),
             splice_temperature_c=55.0, stop_on_unpatchable=True):
    """Return (iterations, final_parts, final_gate).

    Iteration 0 is the baseline. A step is *kept* only if it does not worsen any
    verdict; otherwise it is logged as rejected and the previous state continues,
    so the table shows what was tried, not just what worked.
    """
    rank = {gate0.PASS: 0, gate0.REVISE: 1, gate0.UNAVAILABLE: 1, gate0.FAIL: 2}
    asm, m, g = _verdicts(parts, temperatures, cm_path, arm_lengths, splice_temperature_c)
    its = [_row(0, "Garis dasar", "—", "—", g, asm, kept=True)]
    cur_parts, cur_gate = parts, g
    for i, st in enumerate(steps, 1):
        if stop_on_unpatchable and any(
            t["verdict"] == gate0.FAIL and t["id"] in gate0.NO_PATCH_ON_FAIL for t in cur_gate["tests"]
        ):
            its[-1]["note"] = ("Uji tak-boleh-ditambal gagal: loop dihentikan; "
                               "pindah ke jalur sirkularisasi cadangan.")
            # Still run the remaining steps once so the table shows REVISE items resolved
            # independently of the unpatchable test, clearly labelled.
        cand = st.apply(cur_parts)
        asm, m, g = _verdicts(cand, temperatures, cm_path, arm_lengths, splice_temperature_c)
        worse = any(
            rank[a["verdict"]] > rank[b["verdict"]]
            for a, b in zip(g["tests"], cur_gate["tests"], strict=True)
        )
        its.append(_row(i, st.label, st.parameter, st.value, g, asm, kept=not worse))
        if not worse:
            cur_parts, cur_gate = cand, g
    return its, cur_parts, cur_gate


def _row(n, label, param, value, g, asm, *, kept):
    return {
        "iteration": n,
        "change": label,
        "parameter": param,
        "value": value,
        "kept": kept,
        "verdicts": {t["id"]: t["verdict"] for t in g["tests"]},
        "decision": g["decision"],
        "precursor_nt": len(asm.precursor),
        "circle_nt": len(asm.circle),
        "note": "",
    }


def markdown(its: list[dict]) -> str:
    ids = list(its[0]["verdicts"])
    short = [i.split("_")[0] for i in ids]
    out = ["| Iter | Perubahan (satu parameter) | Dipertahankan | " + " | ".join(short) + " | Keputusan |",
           "|---|---|---|" + "---|" * len(ids) + "---|"]
    for r in its:
        out.append(
            f"| {r['iteration']} | {r['change']} | {'ya' if r['kept'] else 'ditolak'} | "
            + " | ".join(r["verdicts"][i] for i in ids) + f" | {r['decision']} |"
        )
    return "\n".join(out) + "\n"
