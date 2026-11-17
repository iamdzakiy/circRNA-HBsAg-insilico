"""Shared helpers for scripts 01-12: arguments, fixed seed, parameter/version sidecar.

Every script writes a ``*.params.json`` next to its main output recording the
arguments, seed, tool versions and UTC time, so a result can be traced to the run
that produced it. Nothing here invents data: a missing input raises a message
that says what to provide.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

DEFAULT_SEED = 20260101


def parser(description: str, *, orf: bool = True) -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--parts", default="config/parts.example.yaml",
                   help="parts YAML (default: demo config with placeholder sequences)")
    p.add_argument("--results", default="results", help="results directory")
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    if orf:
        p.add_argument("--orf-fasta", dest="orf_fasta", default=None,
                       help="codon-optimised ORF FASTA (output of 02_codon_optimize.py)")
    return p


def versions() -> dict[str, str]:
    from circhbsag import content

    try:
        from circhbsag import __version__ as v
    except Exception:  # pragma: no cover
        v = "unknown"
    return {"circhbsag": v, "python": sys.version.split()[0], **content.detected_versions()}


def sidecar(output: Path, args: argparse.Namespace, extra: dict | None = None) -> Path:
    path = Path(str(output) + ".params.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "script": Path(sys.argv[0]).name,
        "args": {k: v for k, v in vars(args).items()},
        "seed": getattr(args, "seed", None),
        "versions": versions(),
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        **(extra or {}),
    }
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str))
    return path


def load_parts(args):
    from circhbsag import cli

    orf = getattr(args, "orf_fasta", None) or _default_orf(args)
    return cli._load(args.parts, orf)


def _default_orf(args):
    p = Path(args.results) / "demo" / "orf_optimized.fa"
    alt = Path(args.results) / "orf_optimized.fa"
    for c in (alt, p):
        if c.exists():
            return str(c)
    return None


def need(path: Path, hint: str) -> Path:
    if not Path(path).exists():
        raise SystemExit(f"Missing input: {path}\n  -> {hint}")
    return Path(path)


def out_dir(args, sub: str) -> Path:
    d = Path(args.results) / sub
    d.mkdir(parents=True, exist_ok=True)
    return d
