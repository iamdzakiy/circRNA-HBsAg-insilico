"""Pemuatan bagian konstruk dari config/parts.yaml beserta provenansnya.

Format parts.yaml:

    provenance_policy: strict        # strict | allow_synthetic
    protein:
      hbsag_source: data/reference/hbsag_S_protein.fa   # FASTA protein (1 rekaman)
    parts:
      intron_3half:
        file: data/reference/anabaena_intron_3half.fa
        accession: "X58790.1"
        coordinates: "..."
        features: {IGS: [0, 6]}
      E2: {inline: "..."} ...

Setiap bagian WAJIB punya salah satu dari `file` atau `inline`, plus blok
provenans (`accession`/`source`/`note`). Dengan `provenance_policy: strict`,
bagian yang ditandai `synthetic: true` akan DITOLAK, sehingga hasil yang
dilaporkan tidak mungkin berasal dari sekuens buatan.
"""
from __future__ import annotations

from pathlib import Path

import yaml

from .construct import Part
from .seqio import SequenceError, read_fasta


class ProvenanceError(SequenceError):
    """Bagian tidak punya provenans yang memadai untuk kebijakan yang berlaku."""


def load_parts(
    config_path: str | Path, *, root: str | Path | None = None
) -> tuple[dict[str, Part], dict]:
    """Muat bagian konstruk dan kembalikan (parts, metadata_provenans)."""
    config_path = Path(config_path)
    root = Path(root) if root else config_path.resolve().parent.parent
    cfg = yaml.safe_load(config_path.read_text()) or {}
    policy = cfg.get("provenance_policy", "strict")
    raw = cfg.get("parts", {}) or {}

    parts: dict[str, Part] = {}
    prov: dict[str, dict] = {}
    synthetic_used: list[str] = []

    for name, spec in raw.items():
        spec = spec or {}
        if spec.get("synthetic"):
            if policy == "strict":
                raise ProvenanceError(
                    f"bagian '{name}' ditandai synthetic: true tetapi "
                    "provenance_policy = strict. Ganti dengan sekuens terverifikasi, "
                    "atau setel provenance_policy: allow_synthetic HANYA untuk uji perangkat lunak."
                )
            synthetic_used.append(name)
        if "file" in spec:
            path = (root / spec["file"]).resolve()
            if not path.exists():
                raise SequenceError(
                    f"berkas sekuens untuk '{name}' tidak ada: {path}. "
                    "Lihat data/reference/README.md untuk cara mengunduhnya."
                )
            recs = read_fasta(path)
            seq = recs[0].seq
            src = {"file": str(spec["file"]), "fasta_header": recs[0].seq_id, "sha256": recs[0].sha256}
        elif "inline" in spec:
            seq = str(spec["inline"])
            src = {"inline": True}
        elif spec.get("optional"):
            continue
        else:
            raise SequenceError(f"bagian '{name}' tidak punya 'file' maupun 'inline'")

        features = {k: tuple(v) for k, v in (spec.get("features") or {}).items()}
        parts[name] = Part(name=name, seq=seq, role=spec.get("role", ""), features=features)
        prov[name] = {
            **src,
            "length": len(parts[name].seq),
            "accession": spec.get("accession"),
            "source": spec.get("source"),
            "note": spec.get("note"),
            "synthetic": bool(spec.get("synthetic")),
            "provenance_status": "SINTETIS-UJI" if spec.get("synthetic") else "DISEDIAKAN-PENGGUNA",
        }

    meta = {
        "config": str(config_path),
        "provenance_policy": policy,
        "parts": prov,
        "synthetic_parts": synthetic_used,
        "results_may_be_reported_as_scientific": not synthetic_used,
        "protein": cfg.get("protein", {}),
        "design": cfg.get("design", {}),
    }
    return parts, meta


def load_protein(config_path: str | Path, *, root: str | Path | None = None) -> tuple[str, dict]:
    """Muat sekuens protein sasaran (HBsAg) dari FASTA yang ditunjuk parts.yaml."""
    config_path = Path(config_path)
    root = Path(root) if root else config_path.resolve().parent.parent
    cfg = yaml.safe_load(config_path.read_text()) or {}
    spec = (cfg.get("protein") or {})
    if "inline" in spec:
        return str(spec["inline"]).strip().upper(), {"inline": True, "synthetic": bool(spec.get("synthetic"))}
    if "file" not in spec:
        raise SequenceError(
            "blok 'protein' pada parts.yaml harus menunjuk 'file' (FASTA protein HBsAg) "
            "atau 'inline'."
        )
    path = (root / spec["file"]).resolve()
    if not path.exists():
        raise SequenceError(f"FASTA protein tidak ada: {path}")
    text = path.read_text()
    seqs, header, cur = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            header = header or line[1:].strip()
            cur = []
        elif line.strip():
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    if not seqs:
        raise SequenceError(f"tidak ada rekaman protein di {path}")
    protein = seqs[0].upper().replace(" ", "").rstrip("*")
    return protein, {
        "file": spec["file"],
        "header": header,
        "length_aa": len(protein),
        "accession": spec.get("accession"),
        "genotype": spec.get("genotype"),
        "synthetic": bool(spec.get("synthetic")),
    }
