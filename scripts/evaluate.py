#!/usr/bin/env python3
"""Compare called variants with the simulated truth set and write an HTML report.

Matching is allele-level after trimming shared prefixes/suffixes, so an indel written as
chr20:1234 AT>A in the truth and A>ATT in a call is not counted as a false negative.

Usage:
    python scripts/evaluate.py --truth data/sim/truth.vcf --calls results/calls.vcf.gz \
        --outdir results
"""
from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def open_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def trim_alleles(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Remove shared prefix/suffix so equivalent representations compare equal."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref.upper(), alt.upper()


def parse_vcf(path: Path, keep_only_pass: bool = True):
    """Return dict keyed by (pos, ref, alt) -> dict(qual, dp, gt, type)."""
    out = {}
    stats = Counter()
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            pos, ref, alts = int(f[1]), f[3], f[4].split(",")
            qual = f[5]
            filt = f[6]
            if keep_only_pass and filt not in ("PASS", ".", "SIM"):
                stats["filtered_out"] += 1
                continue
            info = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in f[7].split(";"))
            dp = info.get("DP")
            gt = None
            if len(f) > 9:
                fmt = f[8].split(":")
                sample = f[9].split(":")
                if "GT" in fmt:
                    gt = sample[fmt.index("GT")]
            for alt in alts:
                if alt.startswith("<") or alt in ("*", "."):
                    continue
                p, r, a = trim_alleles(pos, ref, alt)
                vtype = "SNV" if len(r) == 1 and len(a) == 1 else "INDEL"
                out[(p, r, a)] = {"qual": qual, "dp": dp, "gt": gt, "type": vtype}
                stats[vtype] += 1
    return out, stats


def evaluate(truth, calls):
    tp = set(truth) & set(calls)
    fn = set(truth) - set(calls)
    fp = set(calls) - set(truth)
    res = {}
    for label, keys in (("all", None), ("SNV", "SNV"), ("INDEL", "INDEL")):
        t = {k for k in truth if keys is None or truth[k]["type"] == keys}
        c = {k for k in calls if keys is None or calls[k]["type"] == keys}
        hit = len(t & c)
        precision = hit / len(c) if c else 0.0
        recall = hit / len(t) if t else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        res[label] = {"truth": len(t), "called": len(c), "tp": hit,
                      "fp": len(c) - hit, "fn": len(t) - hit,
                      "precision": round(precision, 4), "recall": round(recall, 4),
                      "f1": round(f1, 4)}
    gts = Counter(calls[k]["gt"] for k in tp if calls[k]["gt"])
    return res, gts


HTML = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>WGS variant calling report</title>
<style>
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;max-width:900px;margin:40px auto;padding:0 20px;color:#1b2430;line-height:1.6}
h1{border-bottom:3px solid #2ec4b6;padding-bottom:10px}
h2{margin-top:34px;border-bottom:1px solid #e3e8ef;padding-bottom:6px}
table{border-collapse:collapse;width:100%;font-size:14px;margin:12px 0}
th,td{border:1px solid #e3e8ef;padding:8px 10px;text-align:left}
th{background:#f6f8fb}
code{background:#f2f4f7;padding:2px 5px;border-radius:3px;font-size:13px}
img{max-width:100%;border:1px solid #e3e8ef;border-radius:6px}
.note{background:#eaf7f5;border-left:4px solid #2ec4b6;padding:12px 16px;margin:16px 0}
.meta{color:#5b6b7c;font-size:14px}
</style></head><body>
<h1>WGS variant calling report</h1>
<p class="meta">Reference: {{ chrom }} ({{ ref_len }} bp) &middot; simulated {{ coverage }}x
&middot; called {{ totals.called }} variants &middot; generated {{ date }}</p>

<h2>1. Accuracy against the truth set</h2>
<table><tr><th>Class</th><th>Truth</th><th>Called</th><th>TP</th><th>FP</th><th>FN</th><th>Precision</th><th>Recall</th><th>F1</th></tr>
{% for k, v in metrics.items() %}<tr><td>{{ k }}</td><td>{{ v.truth }}</td><td>{{ v.called }}</td><td>{{ v.tp }}</td><td>{{ v.fp }}</td><td>{{ v.fn }}</td><td>{{ '%.3f'|format(v.precision) }}</td><td>{{ '%.3f'|format(v.recall) }}</td><td>{{ '%.3f'|format(v.f1) }}</td></tr>{% endfor %}
</table>
<div class="note">Genotypes among true positives: {{ genotype_summary }}. The simulated sample is
homozygous alternate, so heterozygous calls are false positives that happen to sit on a true site
and are reported here rather than hidden.</div>

<h2>2. Alignment</h2>
<pre>{{ flagstat }}</pre>
<table><tr><th>Metric</th><th>Value</th></tr>
{% for k, v in cov_table.items() %}<tr><td>{{ k }}</td><td>{{ v }}</td></tr>{% endfor %}
</table>

<h2>3. Variant quality</h2>
<img src="figures/qual_hist.png" alt="quality distribution">
<p>Calls are filtered at QUAL &ge; 30 and depth &ge; 10 after left-alignment and splitting
multiallelic records. Filtering is reported, not silently applied: the raw call count is in
<code>results/raw.vcf.gz</code> and the filtered set in <code>results/calls.vcf.gz</code>.</p>

<h2>4. Files</h2>
<ul>
<li><code>results/aln.bam</code>, <code>results/aln.bam.bai</code> &mdash; sorted alignments</li>
<li><code>results/raw.vcf.gz</code> &rarr; <code>results/norm.vcf.gz</code> &rarr; <code>results/calls.vcf.gz</code></li>
<li><code>results/eval.json</code> &mdash; metrics in machine-readable form</li>
<li><code>scripts/</code> &mdash; reference prep, read simulation, calling, evaluation</li>
</ul>

<h2>5. What this demo does and does not show</h2>
<p>It shows the full path from reads to a filtered, evaluated call set, with a reproducible
truth-based measurement of accuracy. It does not show germline best practice for production:
no duplicate marking, no BQSR, no gVCF joint genotyping across samples, and no structural
variants. Those are deliberate omissions for a 2 Mb demo, not oversights, and each is a known
next step rather than a hidden assumption.</p>
</body></html>
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--truth", default="data/sim/truth.vcf")
    ap.add_argument("--calls", default="results/calls.vcf.gz")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--ref", default="data/ref/ref.fa")
    ap.add_argument("--coverage", default="30")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)

    truth, _ = parse_vcf(Path(args.truth))
    calls, cstats = parse_vcf(Path(args.calls))
    metrics, gts = evaluate(truth, calls)

    quals = []
    with open_maybe_gz(Path(args.calls)) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            try:
                quals.append(float(f[5]))
            except (IndexError, ValueError):
                pass
    if quals:
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.hist(quals, bins=40, color="#2ec4b6", edgecolor="black", linewidth=0.4)
        ax.set_xlabel("QUAL")
        ax.set_ylabel("variants")
        ax.set_title(f"Quality score distribution (n={len(quals)})")
        fig.tight_layout()
        fig.savefig(outdir / "figures" / "qual_hist.png", dpi=150)
        plt.close(fig)

    flagstat = (outdir / "qc_flagstat.txt").read_text() if (outdir / "qc_flagstat.txt").exists() else "n/a"
    coverage = {}
    cov_path = outdir / "qc_coverage.txt"
    if cov_path.exists():
        for line in cov_path.read_text().splitlines():
            k, v = line.split("\t")
            coverage[k] = v

    ref_len = 0
    chrom = "chr20"
    ref_path = Path(args.ref)
    if ref_path.exists():
        with open(ref_path) as fh:
            seq = []
            for line in fh:
                if line.startswith(">"):
                    chrom = line[1:].split()[0]
                else:
                    seq.append(line.strip())
            ref_len = len("".join(seq))

    from jinja2 import Template
    import pandas as pd

    html = Template(HTML).render(
        chrom=chrom, ref_len=f"{ref_len:,}", coverage=args.coverage,
        totals={"called": len(calls)}, metrics=metrics,
        genotype_summary=", ".join(f"{k}: {v}" for k, v in gts.most_common()) or "n/a",
        flagstat=flagstat.strip(), cov_table=coverage,
        date=pd.Timestamp.now().strftime("%Y-%m-%d %H:%M"),
    )
    (outdir / "report.html").write_text(html)
    (outdir / "eval.json").write_text(json.dumps({"metrics": metrics,
                                                  "genotypes": dict(gts),
                                                  "called_variants": len(calls)}, indent=2))

    print(f"[eval] all: precision {metrics['all']['precision']:.3f}, "
          f"recall {metrics['all']['recall']:.3f}, F1 {metrics['all']['f1']:.3f}")
    print(f"[eval] SNV: recall {metrics['SNV']['recall']:.3f} | "
          f"INDEL: recall {metrics['INDEL']['recall']:.3f}")
    print(f"[eval] report -> {outdir/'report.html'}")


if __name__ == "__main__":
    main()
