#!/usr/bin/env python3
"""Simulate paired-end reads from a reference slice with a known variant truth set.

Why: variant calling needs a truth set to be evaluated against. Public benchmark data
(GIAB) is 100+ GB, which is the wrong size for a demo. Simulating reads from a real
reference with a known set of SNVs and indels gives a small, fast, fully reproducible
dataset where precision and recall can actually be measured.

Reads are drawn from a haplotype built by applying the truth variants to the reference,
so every truth variant is present in the data by construction. Sequencing errors are
added at a fixed rate, which is what makes the exercise non-trivial.

Usage:
    python scripts/simulate_reads.py --ref data/ref/ref.fa --outdir data/sim \
        --coverage 30 --snvs 400 --indels 60 --seed 7
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
# sequencing error substitution: A<->G, C<->T transitions dominate, plus transversions
SUBSTITUTE = {ord("A"): "GTC", ord("C"): "TAG", ord("G"): "ACT", ord("T"): "CGA"}


def read_fasta(path: Path) -> tuple[str, np.ndarray]:
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
            else:
                chunks.append(line.strip().upper())
    seq = "".join(chunks)
    return name, np.frombuffer(seq.encode(), dtype=np.uint8)


def make_variants(rng: np.random.Generator, ref: np.ndarray, n_snv: int, n_indel: int):
    """Pick SNVs and indels away from the edges; return list of (pos0, ref_str, alt_str)."""
    usable = np.arange(200, len(ref) - 200)
    snv_pos = rng.choice(usable, size=n_snv, replace=False)
    remaining = np.setdiff1d(usable, snv_pos, assume_unique=False)
    indel_pos = rng.choice(remaining, size=n_indel, replace=False)

    variants = []
    for p in snv_pos:
        r = chr(ref[p])
        variants.append((int(p), r, rng.choice([b for b in "ACGT" if b != r])))
    for p in indel_pos:
        if rng.random() < 0.5:                      # deletion of 1-4 bp
            n = int(rng.integers(1, 5))
            variants.append((int(p), "".join(map(chr, ref[p:p + n + 1])), chr(ref[p])))
        else:                                        # insertion of 1-4 bp
            n = int(rng.integers(1, 5))
            ins = "".join(rng.choice(list("ACGT"), size=n))
            variants.append((int(p), chr(ref[p]), chr(ref[p]) + ins))
    variants.sort()
    return variants


def build_haplotype(ref: np.ndarray, variants) -> np.ndarray:
    out = []
    cursor = 0
    for pos, r, a in variants:
        out.append(ref[cursor:pos])
        out.append(np.frombuffer(a.encode(), dtype=np.uint8))
        cursor = pos + len(r)
    out.append(ref[cursor:])
    return np.concatenate(out)


def wrap(seq: np.ndarray, width: int = 60) -> str:
    s = seq.tobytes().decode()
    return "\n".join(s[i:i + width] for i in range(0, len(s), width))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref", default="data/ref/ref.fa")
    ap.add_argument("--outdir", default="data/sim")
    ap.add_argument("--coverage", type=int, default=30)
    ap.add_argument("--read-length", type=int, default=150)
    ap.add_argument("--fragment", type=int, default=350)
    ap.add_argument("--snvs", type=int, default=400)
    ap.add_argument("--indels", type=int, default=60)
    ap.add_argument("--error-rate", type=float, default=0.003)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    chrom, ref = read_fasta(Path(args.ref))
    rng = np.random.default_rng(args.seed)
    variants = make_variants(rng, ref, args.snvs, args.indels)
    hap = build_haplotype(ref, variants)

    # truth VCF
    with open(outdir / "truth.vcf", "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##contig=<ID={chrom},length={len(ref)}>\n")
        fh.write('##INFO=<ID=SIM,Number=0,Type=Flag,Description="simulated truth variant">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i, (pos, r, a) in enumerate(variants, 1):
            fh.write(f"{chrom}\t{pos + 1}\tSIM{i}\t{r}\t{a}\t100\tPASS\tSIM\n")

    # reads
    n_reads = int(len(ref) * args.coverage / (2 * args.read_length))
    frag_len = rng.normal(args.fragment, 40, size=n_reads).astype(int).clip(args.read_length, 800)
    starts = rng.integers(0, len(hap) - frag_len.max() - 1, size=n_reads)

    def sample_reads(r1_start, r1_len):
        idx = r1_start[:, None] + np.arange(args.read_length)[None, :]
        idx = np.clip(idx, 0, len(hap) - 1)
        return hap[idx].copy()

    r1 = sample_reads(starts, args.read_length)
    r2 = sample_reads(starts + frag_len - args.read_length, args.read_length)
    r2 = r2[:, ::-1].copy()
    r2 = np.where(r2 == ord("A"), ord("T"), np.where(r2 == ord("T"), ord("A"),
                 np.where(r2 == ord("C"), ord("G"), np.where(r2 == ord("G"), ord("C"), r2))))

    # sequencing errors: independent per base
    def add_errors(reads):
        mask = rng.random(reads.shape) < args.error_rate
        rows, cols = np.nonzero(mask)
        for r, c in zip(rows, cols):
            b = chr(reads[r, c])
            if b in SUBSTITUTE:
                reads[r, c] = ord(rng.choice(list(SUBSTITUTE[b])))

    add_errors(r1)
    add_errors(r2)

    with gzip.open(outdir / "reads_1.fastq.gz", "wt") as f1, gzip.open(outdir / "reads_2.fastq.gz", "wt") as f2:
        for i in range(n_reads):
            f1.write(f"@sim{i}/1\n{r1[i].tobytes().decode()}\n+\n{'I' * args.read_length}\n")
            f2.write(f"@sim{i}/2\n{r2[i].tobytes().decode()}\n+\n{'I' * args.read_length}\n")

    with open(outdir / "hap.fa", "w") as fh:
        fh.write(f">{chrom}_hap\n{wrap(hap)}\n")

    n_snv = sum(1 for _, r, a in variants if len(r) == len(a) == 1)
    print(f"chrom {chrom}: ref {len(ref):,} bp, haplotype {len(hap):,} bp")
    print(f"truth: {n_snv} SNVs, {len(variants) - n_snv} indels -> {outdir/'truth.vcf'}")
    print(f"reads: {n_reads:,} pairs x {args.read_length} bp (~{args.coverage}x) -> {outdir}/reads_*.fastq.gz")


if __name__ == "__main__":
    main()
