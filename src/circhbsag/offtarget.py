"""Pemeriksaan homologi off-target dengan BLAST+ (lisensi: domain publik, NCBI).

BLAST+ membutuhkan basis data sekuens dan akses jaringan untuk mengunduhnya,
sehingga modul ini dirancang agar:
  - mendeteksi apakah `blastn`/`makeblastdb` ada di PATH,
  - menjalankan BLAST sungguhan bila biner DAN basis data tersedia,
  - dan bila tidak, mengembalikan status "SKIPPED" yang eksplisit beserta
    perintah yang harus dijalankan pengguna. TIDAK ADA angka yang dikarang.

Selain BLAST, modul ini menyediakan pemeriksaan homologi mandiri berbasis
k-mer terhadap berkas FASTA lokal mana pun yang disediakan pengguna (mis. CDS
manusia), yang berguna untuk penyaringan cepat tanpa basis data BLAST.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

BLAST_LICENSE = "NCBI BLAST+ — domain publik (karya pemerintah AS), gratis"


def blast_available() -> dict:
    return {
        "blastn": shutil.which("blastn"),
        "makeblastdb": shutil.which("makeblastdb"),
        "license": BLAST_LICENSE,
    }


def run_blastn(
    query_fasta: str | Path,
    db: str,
    *,
    evalue: float = 1e-3,
    word_size: int = 11,
    max_target_seqs: int = 50,
    threads: int = 2,
    task: str = "blastn",
) -> dict:
    """Jalankan blastn bila tersedia; kembalikan hasil terparsing atau status SKIPPED."""
    av = blast_available()
    if not av["blastn"]:
        return {
            "status": "SKIPPED",
            "reason": "blastn tidak ada di PATH",
            "how_to_fix": (
                "conda install -c bioconda blast  # atau: apt-get install ncbi-blast+\n"
                "lalu unduh basis data, mis.:\n"
                "  update_blastdb.pl --decompress refseq_rna\n"
                "atau buat basis data lokal dari CDS manusia:\n"
                "  makeblastdb -in GRCh38_cds.fa -dbtype nucl -out human_cds"
            ),
            "license": BLAST_LICENSE,
        }
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
    cmd = [
        av["blastn"], "-task", task, "-query", str(query_fasta), "-db", db,
        "-evalue", str(evalue), "-word_size", str(word_size),
        "-max_target_seqs", str(max_target_seqs), "-num_threads", str(threads),
        "-outfmt", outfmt,
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1800, check=False)
    except (subprocess.TimeoutExpired, OSError) as exc:
        return {"status": "ERROR", "reason": str(exc), "cmd": " ".join(cmd)}
    if proc.returncode != 0:
        return {
            "status": "ERROR",
            "returncode": proc.returncode,
            "stderr": proc.stderr[-2000:],
            "cmd": " ".join(cmd),
            "hint": "Basis data mungkin tidak ada. Periksa $BLASTDB dan nama -db.",
        }
    cols = outfmt.split()[1:]
    hits = []
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < len(cols):
            continue
        rec = dict(zip(cols, parts, strict=False))
        for k in ("pident", "evalue", "bitscore"):
            rec[k] = float(rec[k])
        for k in ("length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send"):
            rec[k] = int(rec[k])
        hits.append(rec)
    return {
        "status": "OK",
        "cmd": " ".join(cmd),
        "n_hits": len(hits),
        "hits": hits,
        "license": BLAST_LICENSE,
    }


def kmer_homology(
    query: str,
    reference_fasta: str | Path,
    *,
    k: int = 21,
    min_shared: int = 1,
) -> dict:
    """Penyaringan homologi mandiri: cari k-mer persis yang dibagi dengan referensi.

    k = 21 dipilih karena k-mer 21 nt persis hampir tidak pernah muncul secara
    kebetulan pada transkriptom manusia (4^21 >> 3e9), sehingga satu k-mer yang
    sama sudah menandakan homologi lokal yang perlu diperiksa. Hasil berupa
    daftar wilayah query yang cocok, bukan nilai statistik BLAST.
    """
    from .seqio import read_fasta, revcomp

    q = query.upper().replace("U", "T")
    qk: dict[str, list[int]] = defaultdict(list)
    for i in range(len(q) - k + 1):
        qk[q[i : i + k]].append(i)

    shared: dict[str, list[dict]] = defaultdict(list)
    for rec in read_fasta(reference_fasta):
        s = rec.seq
        for strand, seq in (("+", s), ("-", revcomp(s))):
            for i in range(len(seq) - k + 1):
                km = seq[i : i + k]
                if km in qk:
                    for qpos in qk[km]:
                        shared[rec.seq_id].append(
                            {"query_pos": qpos, "ref_pos": i, "strand": strand, "kmer": km}
                        )
    out = {
        "status": "OK",
        "k": k,
        "reference": str(reference_fasta),
        "n_reference_sequences_with_hits": len(shared),
        "hits_per_reference": {kk: len(v) for kk, v in sorted(shared.items(), key=lambda x: -len(x[1]))[:25]},
        "flag_off_target": any(len(v) >= min_shared for v in shared.values()),
        "examples": {kk: v[:5] for kk, v in list(shared.items())[:5]},
    }
    return out


def check_construct(
    asm,
    *,
    db: str | None = None,
    reference_fasta: str | Path | None = None,
    out_dir: str | Path = "results/offtarget",
) -> dict:
    """Pemeriksaan off-target untuk lingkaran dan ORF.

    Prioritas: BLAST+ bila `db` diberikan dan biner ada; selain itu penyaringan
    k-mer bila `reference_fasta` diberikan; selain itu status SKIPPED eksplisit.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    queries = {"circle": asm.circle}
    if "ORF" in asm.circle_coords:
        queries["ORF"] = asm.subseq("circle", "ORF")

    results: dict[str, dict] = {}
    for name, seq in queries.items():
        if db and blast_available()["blastn"]:
            with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
                fh.write(f">{name}\n{seq}\n")
                qpath = fh.name
            results[name] = run_blastn(qpath, db)
            results[name]["method"] = "blastn"
        elif reference_fasta:
            results[name] = kmer_homology(seq, reference_fasta)
            results[name]["method"] = "kmer_exact"
        else:
            results[name] = {
                "status": "SKIPPED",
                "method": None,
                "reason": (
                    "tidak ada basis data BLAST (-db) maupun FASTA referensi lokal. "
                    "Pemeriksaan homologi off-target BELUM dilakukan."
                ),
                "how_to_fix": blast_available(),
            }
    (out_dir / "offtarget.json").write_text(json.dumps(results, indent=2))
    return results
