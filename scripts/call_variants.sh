#!/usr/bin/env bash
# Alignment and variant calling for the simulated WGS demo.
# bwa-mem -> sorted BAM -> bcftools mpileup/call -> normalise -> hard filter.
set -euo pipefail

REF="${REF:-data/ref/ref.fa}"
SIM="${SIM:-data/sim}"
OUT="${OUT:-results}"
THREADS="${THREADS:-2}"
MIN_QUAL="${MIN_QUAL:-30}"
MIN_DP="${MIN_DP:-10}"
SAMPLE="${SAMPLE:-SIM1}"
mkdir -p "$OUT"

if [ ! -s "$REF.bwt" ]; then
  echo "[call] reference not indexed, run scripts/get_reference.sh first" >&2
  exit 1
fi

echo "[call] aligning with bwa-mem ($THREADS threads)"
bwa mem -t "$THREADS" -R "@RG\tID:sim1\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1" \
  "$REF" "$SIM/reads_1.fastq.gz" "$SIM/reads_2.fastq.gz" 2> "$OUT/bwa.log" \
  | samtools sort -@ "$THREADS" -o "$OUT/aln.bam" -
samtools index "$OUT/aln.bam"
samtools flagstat "$OUT/aln.bam" > "$OUT/qc_flagstat.txt"

echo "[call] coverage"
samtools depth -a "$OUT/aln.bam" | awk '{s+=$3; n++; if($3>=10) c++} END {printf "mean_coverage\t%.2f\nbases_ge10x\t%d\nbases_total\t%d\n", s/n, c, n}' \
  > "$OUT/qc_coverage.txt"

echo "[call] mpileup + call"
bcftools mpileup -f "$REF" -a AD,DP -q 20 -Q 20 -Ou "$OUT/aln.bam" \
  | bcftools call -mv -Oz -o "$OUT/raw.vcf.gz" 2> "$OUT/call.log"
bcftools index -f "$OUT/raw.vcf.gz"

echo "[call] normalise (left-align, split multiallelics)"
bcftools norm -f "$REF" -m -any "$OUT/raw.vcf.gz" -Oz -o "$OUT/norm.vcf.gz" 2>> "$OUT/call.log"
bcftools index -f "$OUT/norm.vcf.gz"

echo "[call] hard filter: QUAL>=${MIN_QUAL}, DP>=${MIN_DP}"
bcftools filter -i "QUAL>=${MIN_QUAL} && INFO/DP>=${MIN_DP}" "$OUT/norm.vcf.gz" -Oz -o "$OUT/calls.vcf.gz" 2>> "$OUT/call.log"
bcftools index -f "$OUT/calls.vcf.gz"

bcftools stats "$OUT/calls.vcf.gz" > "$OUT/qc_bcftools_stats.txt"
echo "[call] raw:    $(bcftools view -H "$OUT/raw.vcf.gz" | wc -l) records"
echo "[call] passed: $(bcftools view -H "$OUT/calls.vcf.gz" | wc -l) records -> $OUT/calls.vcf.gz"
