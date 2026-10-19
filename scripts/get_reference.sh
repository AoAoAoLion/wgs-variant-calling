#!/usr/bin/env bash
# Fetch a small reference slice and index it for bwa/samtools.
# Everything comes from public sources; no data ships with this repository.
set -euo pipefail

REGION="${REGION:-chr20:1000000-3000000}"
OUTDIR="${OUTDIR:-data/ref}"
mkdir -p "$OUTDIR"

if [ ! -s "$OUTDIR/ref.fa" ]; then
  echo "[ref] downloading hg38 chr20 from UCSC (~63 MB compressed)"
  curl -sSL "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz" \
    | gunzip -c > "$OUTDIR/chr20.fa"
  echo "[ref] slicing $REGION"
  samtools faidx "$OUTDIR/chr20.fa"
  # rename the slice and upper-case soft-masked bases so alleles compare cleanly
  samtools faidx "$OUTDIR/chr20.fa" "$REGION" \
    | sed '1s/.*/>chr20_slice/' \
    | awk '/^>/{print; next} {print toupper($0)}' > "$OUTDIR/ref.fa"
  rm -f "$OUTDIR/chr20.fa" "$OUTDIR/chr20.fa.fai"
fi

if [ ! -s "$OUTDIR/ref.fa.bwt" ]; then
  echo "[ref] indexing with bwa"
  bwa index "$OUTDIR/ref.fa"
fi
samtools faidx "$OUTDIR/ref.fa"

echo "[ref] ready:"
samtools faidx "$OUTDIR/ref.fa" | awk '{print "  " $1 " length=" $2}'
