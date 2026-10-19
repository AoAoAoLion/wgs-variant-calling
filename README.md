# WGS variant calling, reads to evaluated call set

A small, complete germline variant-calling pipeline: reference prep, read simulation with a
known truth set, alignment, variant calling, normalisation, hard filtering, and a measured
accuracy report.

Runs on a 2 Mb slice of human chr20 and finishes in under a minute on a laptop, so the whole
path is inspectable instead of being buried in a 30-hour run.

```bash
make ref     # download chr20 slice from UCSC, index for bwa (~63 MB download, once)
make sim     # simulate 30x paired-end reads + truth VCF
make call    # bwa-mem -> bcftools mpileup/call -> norm -> filter
make eval    # compare against truth, write results/report.html
```

## Result on the bundled simulation

| Class | Truth | Called | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---|---|---|---|---|---|
| all | 460 | 458 | 447 | 11 | 13 | 0.976 | 0.972 | 0.974 |
| SNV | 400 | 399 | 399 | 0 | 1 | 1.000 | 0.998 | 0.999 |
| INDEL | 60 | 59 | 48 | 11 | 12 | 0.814 | 0.800 | 0.807 |

SNVs are essentially solved at 30x; indels are the honest weak spot, which is what indel
calling looks like with a default bcftools setup. Reporting that instead of tuning the demo
until the number looks good is the point.

## Why simulate reads instead of shipping a public dataset

GIAB truth data is the right benchmark and the wrong size for a demo: 100+ GB of FASTQ, hours
of compute, and a reference bundle larger than most people will download to look at a repo.

Simulating from a real reference gives:

- a small, fast dataset that anyone can regenerate byte-for-byte (`--seed 7`);
- a truth set that is guaranteed complete, so precision and recall are real measurements
  rather than assertions;
- a documented error model (0.3% substitution errors) so the difficulty is explicit.

The trade-off is real: simulated reads lack the artefacts of real sequencing (GC bias,
duplicate structure, mapping ambiguity in repeats). The pipeline code paths are the same, and
pointing it at real FASTQ is a matter of swapping `data/sim/reads_*.fastq.gz`.

## Pipeline steps

1. **Reference** — `scripts/get_reference.sh` downloads hg38 chr20, slices
   `chr20:1000000-3000000`, renames the contig to `chr20_slice`, upper-cases the sequence and
   builds the bwa index.
2. **Simulation** — `scripts/simulate_reads.py` draws 400 SNVs and 60 indels, builds a
   haplotype by applying them to the reference, then samples 200,000 read pairs (150 bp,
   ~350 bp fragments) from that haplotype with 0.3% substitution errors.
3. **Alignment** — `bwa mem` with a read group, sorted and indexed with `samtools`.
4. **Calling** — `bcftools mpileup` (min mapQ 20, min baseQ 20) piped into `bcftools call -mv`.
5. **Normalisation** — `bcftools norm -f ref -m -any`: left-align indels and split
   multiallelic records so alleles can be compared one to one.
6. **Filtering** — `QUAL >= 30` and `DP >= 10`. Raw calls are kept
   (`results/raw.vcf.gz`), so the filter is auditable rather than destructive.
7. **Evaluation** — `scripts/evaluate.py` trims shared prefixes/suffixes before matching, so
   `AT > A` and `A > ATT` at the same site are recognised as the same indel. Precision,
   recall and F1 are computed per variant class, and the HTML report is written from them.

## Outputs

| File | Content |
|---|---|
| `results/aln.bam`, `.bai` | sorted alignments |
| `results/raw.vcf.gz` → `norm.vcf.gz` → `calls.vcf.gz` | call sets at each stage |
| `results/qc_flagstat.txt`, `qc_coverage.txt` | alignment and depth summary |
| `results/qc_bcftools_stats.txt` | bcftools stats |
| `results/eval.json` | metrics, machine-readable |
| `results/report.html` | self-contained report with the accuracy table |

## Requirements

`bwa`, `samtools`, `bcftools` (apt: `sudo apt install bwa samtools bcftools`) and Python with
`numpy matplotlib jinja2 pandas` (`pip install -r requirements.txt`).

## What this pipeline deliberately does not do

- **No duplicate marking.** Simulated reads have no PCR duplicates; adding the step would be
  cargo cult. Real data needs it.
- **No BQSR.** Base quality recalibration needs a known-sites resource and a real error
  profile. It matters for real data and is dishonest to fake here.
- **No joint genotyping.** Single simulated sample, so there is nothing to joint-call.
- **No structural variants, no CNVs.** Different tool chain, different truth set.
- **No truth filtering of calls.** Every call is compared; no "known sites only" shortcut.

Each omission is a documented decision, not an oversight, and each is the obvious next step
for a production germline pipeline.

## Layout

```
scripts/get_reference.sh   reference download, slice, index
scripts/simulate_reads.py  truth set + paired-end read simulation
scripts/call_variants.sh   alignment, calling, normalisation, filtering
scripts/evaluate.py        truth comparison + HTML report
data/                      reference and simulated reads (gitignored, generated)
results/                   outputs (gitignored, generated)
tests/test_pipeline.py     end-to-end test on a reduced dataset
```

## License

MIT.
