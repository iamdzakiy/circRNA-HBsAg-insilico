"""Filter statistik sambungan balik (BSJ) ala Szabo & Salzman.

Gagasan yang diadopsi [LITERATUR]: sebuah bacaan yang melintasi sambungan balik
hanya meyakinkan bila (i) jangkar (anchor) di kedua sisi sambungan cukup panjang,
(ii) jumlah mismatch rendah, dan (iii) pemetaannya unik — yaitu jangkarnya tidak
sama baiknya menempel di tempat lain pada referensi linear. Modul ini memberi
nomor pada ketiganya lewat posterior Bayes per bacaan, lalu menggabungkannya
menjadi dukungan tingkat-sambungan.

Implementasi mandiri (MIT) dipakai agar pipeline dapat berjalan tanpa jaringan.
Bila biner STAR / CIRI2 / CIRCexplorer tersedia di PATH, `external_available()`
melaporkannya dan `write_star_cmd()` mencetak perintah yang setara untuk
validasi silang; keduanya opsional dan tidak pernah wajib.
"""
from __future__ import annotations

import math
import random
import shutil
from dataclasses import dataclass

from .seqio import revcomp

EXTERNAL_TOOLS = {
    "STAR": "STAR (MIT) — chimeric alignment untuk deteksi BSJ",
    "CIRI2": "CIRI2 (GPL) — pemanggil circRNA berbasis BWA-MEM",
    "CIRCexplorer2": "CIRCexplorer2 (MIT) — anotasi BSJ",
    "bwa": "BWA (GPLv3) — prasyarat CIRI2",
    "art_illumina": "ART (GPL) — simulator bacaan Illumina",
    "wgsim": "wgsim (MIT) — simulator bacaan sederhana",
}


def external_available() -> dict[str, str | None]:
    """Biner eksternal opsional yang terpasang di PATH."""
    return {k: shutil.which(k) for k in EXTERNAL_TOOLS}


# --------------------------------------------------------------------------
# Simulasi bacaan
# --------------------------------------------------------------------------


@dataclass
class Read:
    name: str
    seq: str
    origin: str  # "circular" | "linear"
    true_crosses_bsj: bool
    n_errors: int


def simulate_reads(
    circle: str,
    linear: str,
    *,
    n_circular: int = 2000,
    n_linear: int = 2000,
    read_len: int = 150,
    error_rate: float = 0.001,
    seed: int = 20260101,
) -> list[Read]:
    """Simulator bacaan single-end sederhana dengan galat substitusi i.i.d.

    Bacaan dari lingkaran diambil dengan koordinat modular, sehingga sebagian
    melintasi sambungan balik secara alami dengan proporsi read_len/len(circle).
    Bacaan dari templat linear (prekursor/mRNA linear) berfungsi sebagai latar
    yang tidak boleh menghasilkan dukungan BSJ.
    Bila `art_illumina` atau `wgsim` tersedia, keduanya lebih disukai untuk
    publikasi; fungsi ini dipakai agar pipeline tetap jalan tanpa jaringan.
    """
    rng = random.Random(seed)
    reads: list[Read] = []
    n = len(circle)
    bsj_zone = n  # posisi sambungan ada di batas indeks n-1 -> 0
    for i in range(n_circular):
        start = rng.randrange(n)
        seq = "".join(circle[(start + k) % n] for k in range(read_len))
        crosses = (start + read_len) > bsj_zone
        seq, e = _mutate(seq, error_rate, rng)
        if rng.random() < 0.5:
            seq = revcomp(seq)
        reads.append(Read(f"circ_{i}", seq, "circular", crosses, e))
    for i in range(n_linear):
        if len(linear) <= read_len:
            seq = linear
        else:
            start = rng.randrange(len(linear) - read_len)
            seq = linear[start : start + read_len]
        seq, e = _mutate(seq, error_rate, rng)
        if rng.random() < 0.5:
            seq = revcomp(seq)
        reads.append(Read(f"lin_{i}", seq, "linear", False, e))
    rng.shuffle(reads)
    return reads


def _mutate(seq: str, rate: float, rng: random.Random) -> tuple[str, int]:
    if rate <= 0:
        return seq, 0
    chars, n_err = list(seq), 0
    for i, b in enumerate(chars):
        if rng.random() < rate:
            chars[i] = rng.choice([x for x in "ACGT" if x != b])
            n_err += 1
    return "".join(chars), n_err


# --------------------------------------------------------------------------
# Pemanggilan BSJ
# --------------------------------------------------------------------------


def _best_match(hay: str, needle: str, max_mm: int) -> tuple[int, int] | None:
    """Pencarian kecocokan terbaik needle di hay dengan toleransi mismatch."""
    best, best_mm = None, max_mm + 1
    ln = len(needle)
    for i in range(len(hay) - ln + 1):
        mm = 0
        for a, b in zip(hay[i : i + ln], needle, strict=False):
            if a != b:
                mm += 1
                if mm >= best_mm:
                    break
        if mm < best_mm:
            best, best_mm = i, mm
            if mm == 0:
                break
    return (best, best_mm) if best is not None else None


def call_bsj(
    reads: list[Read],
    circle: str,
    linear: str,
    *,
    anchor: int = 20,
    max_mismatch: int = 2,
) -> dict:
    """Panggil bacaan pendukung BSJ dengan uji jangkar dua sisi dan uji keunikan.

    Untuk setiap bacaan (dan komplemen baliknya), diambil jangkar `anchor` nt
    dari kedua ujung. Bacaan dinyatakan kandidat BSJ bila:
      - jangkar kiri menempel di sisi 5'-sambungan lingkaran (ujung 3' E1),
      - jangkar kanan menempel di sisi 3'-sambungan (awal E2), dan
      - urutan bacaan secara keseluruhan TIDAK menempel ke referensi linear
        dalam satu blok beruntun (uji keunikan/anti-linear).
    """
    n = len(circle)
    # jendela referensi di kedua sisi sambungan (sambungan ada di antara n-1 dan 0)
    left_ref = circle[max(0, n - 120) :]
    right_ref = circle[: min(n, 120)]
    results: list[dict] = []
    ambiguous: list[dict] = []
    for r in reads:
        for orient, s in (("fwd", r.seq), ("rev", revcomp(r.seq))):
            if len(s) < 2 * anchor:
                continue
            la, ra = s[:anchor], s[-anchor:]
            hitl = _best_match(left_ref, la, max_mismatch)
            hitr = _best_match(right_ref, ra, max_mismatch)
            if not hitl or not hitr:
                continue
            # jangkar kiri harus dekat ujung 3' left_ref, kanan dekat awal right_ref
            if hitl[0] + anchor < len(left_ref) - 100 or hitr[0] > 100:
                continue
            lin_hit = _best_match(linear, s, max_mismatch)
            overlap_left = len(left_ref) - hitl[0]
            overlap_right = hitr[0] + anchor
            (ambiguous if lin_hit is not None else results).append(
                {
                    "read": r.name,
                    "orientation": orient,
                    "origin_true": r.origin,
                    "anchor_left_mismatch": hitl[1],
                    "anchor_right_mismatch": hitr[1],
                    "overlap_left_nt": min(overlap_left, len(s) - anchor),
                    "overlap_right_nt": min(overlap_right, len(s) - anchor),
                    "maps_to_linear_contiguously": lin_hit is not None,
                    "linear_mismatch": lin_hit[1] if lin_hit else None,
                }
            )
            break
    return {
        "anchor": anchor,
        "max_mismatch": max_mismatch,
        "n_reads": len(reads),
        "n_candidates": len(results),
        "candidates": results,
        # Uji KEUNIKAN (kriteria ketiga Szabo): bacaan yang seluruhnya dapat
        # dijelaskan oleh referensi LINEAR bukan bukti sirkularisasi, berapa pun
        # panjang jangkarnya, karena satu pengulangan sekuens sudah cukup
        # menghasilkan pola seperti sambungan balik. Bacaan semacam itu
        # DIKELUARKAN dari penghitungan posterior dan dilaporkan terpisah.
        "n_excluded_linear_explained": len(ambiguous),
        "excluded_linear_explained": ambiguous[:50],
        "excluded_by_origin": {
            "circular": sum(1 for c in ambiguous if c["origin_true"] == "circular"),
            "linear": sum(1 for c in ambiguous if c["origin_true"] == "linear"),
        },
    }


# --------------------------------------------------------------------------
# Posterior Bayes per bacaan
# --------------------------------------------------------------------------


def read_posterior(
    cand: dict,
    *,
    prior_circular: float = 0.01,
    anchor: int = 20,
    read_len: int = 150,
    error_rate: float = 0.001,
    genome_complexity: float = 3.1e9,
) -> dict:
    """P(sirkular | bukti) untuk satu bacaan kandidat.

    Model likelihood (sederhana dan eksplisit):

      Bukti  = (overlap jangkar kiri o_L, kanan o_R; mismatch m; menempel linear?)

      Di bawah H1 (bacaan benar dari lingkaran):
        - mismatch mengikuti Binomial(read_len, error_rate)
        - overlap jangkar: bacaan yang melintasi sambungan memberi overlap >= anchor
          dengan probabilitas ~1 (sudah disyaratkan oleh pemanggil)
        - tidak menempel ke referensi linear secara beruntun: probabilitas ~1

      Di bawah H0 (artefak: template switching RT, sirkularisasi ligasi,
      atau pemetaan keliru):
        - peluang dua jangkar sepanjang o_L dan o_R kebetulan bersebelahan dalam
          orientasi sambungan balik ~ 4^-(o_L + o_R) dikalikan kompleksitas
          ruang pencarian `genome_complexity`
        - mismatch cenderung lebih tinggi: dimodelkan Binomial dengan laju 5x

    Keluaran: posterior per bacaan dan log10 Bayes factor. Nilai prior dan
    `genome_complexity` adalah parameter yang WAJIB dilaporkan bersama hasil;
    keduanya [PERKIRAAN] dan dapat diubah.
    """
    o_l = max(anchor, cand["overlap_left_nt"])
    o_r = max(anchor, cand["overlap_right_nt"])
    m = cand["anchor_left_mismatch"] + cand["anchor_right_mismatch"]

    def binom_logpmf(k: int, n: int, p: float) -> float:
        p = min(max(p, 1e-12), 1 - 1e-12)
        return (
            math.lgamma(n + 1)
            - math.lgamma(k + 1)
            - math.lgamma(n - k + 1)
            + k * math.log(p)
            + (n - k) * math.log(1 - p)
        )

    ll_h1 = binom_logpmf(m, read_len, error_rate)
    ll_h0 = binom_logpmf(m, read_len, min(0.5, error_rate * 5))
    # komponen kebetulan: peluang jangkar gabungan muncul secara acak
    ll_h0 += math.log(genome_complexity) - (o_l + o_r) * math.log(4)
    if cand["maps_to_linear_contiguously"]:
        # Bacaan yang dapat dijelaskan templat linear seharusnya sudah
        # dikeluarkan oleh `call_bsj`. Bila fungsi ini dipanggil langsung atas
        # bacaan semacam itu, likelihood H1 ditekan kuat: sebuah bacaan BSJ asli
        # hanya bisa juga memetakan beruntun ke referensi linear jika ada
        # pengulangan sekuens, yang jarang dan harus diperiksa manual.
        ll_h1 += math.log(1e-9)
    log_bf = (ll_h1 - ll_h0) / math.log(10)
    odds_prior = prior_circular / (1 - prior_circular)
    log_post_odds = math.log10(odds_prior) + log_bf
    post = 1.0 / (1.0 + 10 ** (-log_post_odds))
    return {
        "read": cand["read"],
        "overlap_left_nt": o_l,
        "overlap_right_nt": o_r,
        "mismatch_total": m,
        "maps_to_linear": cand["maps_to_linear_contiguously"],
        "log10_bayes_factor": round(log_bf, 3),
        "posterior_circular": round(post, 6),
        "prior_circular": prior_circular,
    }


def junction_support(
    calls: dict,
    *,
    prior_circular: float = 0.01,
    posterior_threshold: float = 0.95,
    read_len: int = 150,
    error_rate: float = 0.001,
) -> dict:
    """Agregasi posterior per bacaan menjadi dukungan tingkat-sambungan.

    Dilaporkan: jumlah bacaan lolos ambang, posterior gabungan (digabung dalam
    ruang log-odds dengan asumsi bacaan independen), dan tabel kebenaran
    terhadap asal bacaan simulasi (memungkinkan pengukuran FPR/TPR filter).
    """
    posts = [
        read_posterior(
            c, prior_circular=prior_circular, anchor=calls["anchor"],
            read_len=read_len, error_rate=error_rate
        )
        for c in calls["candidates"]
    ]
    passing = [p for p in posts if p["posterior_circular"] >= posterior_threshold]
    tp = sum(1 for p, c in zip(posts, calls["candidates"], strict=False) if c["origin_true"] == "circular" and p["posterior_circular"] >= posterior_threshold)
    fp = sum(1 for p, c in zip(posts, calls["candidates"], strict=False) if c["origin_true"] == "linear" and p["posterior_circular"] >= posterior_threshold)
    fn = sum(1 for p, c in zip(posts, calls["candidates"], strict=False) if c["origin_true"] == "circular" and p["posterior_circular"] < posterior_threshold)
    log_odds = 0.0
    for p in passing:
        log_odds += p["log10_bayes_factor"]
    if passing:
        log_odds += math.log10(prior_circular / (1 - prior_circular))
    joint = 1.0 / (1.0 + 10 ** (-log_odds)) if passing else 0.0
    return {
        "n_candidates": len(posts),
        "n_excluded_linear_explained": calls.get("n_excluded_linear_explained", 0),
        "n_passing": len(passing),
        "posterior_threshold": posterior_threshold,
        "prior_circular": prior_circular,
        "joint_posterior_circular": round(joint, 8),
        "log10_joint_bayes_factor": round(log_odds, 3),
        "confusion_on_simulated_origin": {
            "true_positive_circular": tp,
            "false_positive_from_linear": fp,
            "false_negative_circular": fn,
            "precision": round(tp / (tp + fp), 4) if (tp + fp) else None,
        },
        "per_read": posts[:50],
        "note": (
            "Angka presisi/recall di sini mengukur FILTER STATISTIK pada bacaan "
            "SIMULASI, bukan kinerja eksperimen. Pada data nyata, jalankan juga "
            "STAR chimeric atau CIRI2 (lihat external_available()) sebagai validasi silang. "
            "Bacaan yang dapat dijelaskan referensi linear dikeluarkan lebih dahulu; "
            "bila jumlah yang dikeluarkan besar, konstruk kemungkinan memuat "
            "pengulangan sekuens yang menurunkan daya deteksi BSJ dan perlu didesain ulang."
        ),
    }


def write_star_cmd(fastq: str, genome_dir: str, out_prefix: str, *, threads: int = 4) -> str:
    """Perintah STAR chimeric yang setara, untuk dijalankan bila STAR terpasang."""
    return (
        f"STAR --runThreadN {threads} --genomeDir {genome_dir} --readFilesIn {fastq} "
        f"--outFileNamePrefix {out_prefix} --outSAMtype BAM Unsorted "
        "--chimSegmentMin 20 --chimScoreMin 1 --chimScoreJunctionNonGTAG 0 "
        "--chimScoreDropMax 20 --chimScoreSeparation 10 --alignSJoverhangMin 15 "
        "--alignSJDBoverhangMin 10 --alignMatesGapMax 0 --alignIntronMax 0 "
        "--outFilterMultimapNmax 2 --chimOutType Junctions SeparateSAMold"
    )


def run_simulation(
    circle: str,
    linear: str,
    *,
    n_circular: int = 1500,
    n_linear: int = 1500,
    read_len: int = 150,
    error_rate: float = 0.001,
    anchor: int = 20,
    max_mismatch: int = 2,
    prior_circular: float = 0.01,
    seed: int = 20260101,
) -> dict:
    """Pipeline lengkap: simulasi bacaan -> pemanggilan BSJ -> posterior Bayes."""
    reads = simulate_reads(
        circle, linear, n_circular=n_circular, n_linear=n_linear,
        read_len=read_len, error_rate=error_rate, seed=seed,
    )
    calls = call_bsj(reads, circle, linear, anchor=anchor, max_mismatch=max_mismatch)
    support = junction_support(
        calls, prior_circular=prior_circular, read_len=read_len, error_rate=error_rate
    )
    expected_cross = n_circular * min(1.0, read_len / len(circle))
    return {
        "parameters": {
            "n_circular": n_circular,
            "n_linear": n_linear,
            "read_len": read_len,
            "error_rate": error_rate,
            "anchor": anchor,
            "max_mismatch": max_mismatch,
            "prior_circular": prior_circular,
            "seed": seed,
            "circle_length": len(circle),
            "linear_length": len(linear),
        },
        "expected_bsj_crossing_reads": round(expected_cross, 1),
        "calls": {k: v for k, v in calls.items() if k != "candidates"},
        "support": support,
        "external_tools_available": external_available(),
    }
