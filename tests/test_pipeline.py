"""End-to-end test for the WGS variant-calling demo.

Runs on a reduced simulation (8x, 60 variants) so the suite stays fast; the production
example in the README uses 30x and 460 variants.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
PY = sys.executable


@pytest.fixture(scope="module")
def called_run(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("wgs")
    ref_dir = REPO / "data" / "ref"
    if not (ref_dir / "ref.fa").exists():
        pytest.skip("reference not built; run 'make ref' first")

    sim = tmp / "sim"
    subprocess.run([PY, str(REPO / "scripts" / "simulate_reads.py"),
                    "--ref", str(ref_dir / "ref.fa"), "--outdir", str(sim),
                    "--coverage", "8", "--snvs", "60", "--indels", "10", "--seed", "11"],
                   check=True, capture_output=True)

    out = tmp / "results"
    env = {"REF": str(ref_dir / "ref.fa"), "SIM": str(sim), "OUT": str(out),
           "THREADS": "2", "MIN_DP": "5", "PATH": "/usr/local/bin:/usr/bin:/bin"}
    proc = subprocess.run(["bash", str(REPO / "scripts" / "call_variants.sh")],
                          env=env, check=True, capture_output=True, text=True)

    subprocess.run([PY, str(REPO / "scripts" / "evaluate.py"),
                    "--truth", str(sim / "truth.vcf"),
                    "--calls", str(out / "calls.vcf.gz"),
                    "--outdir", str(out), "--ref", str(ref_dir / "ref.fa"),
                    "--coverage", "8"],
                   check=True, capture_output=True)
    return tmp, out, proc.stdout


def test_variants_were_called(called_run):
    _, out, _ = called_run
    metrics = json.loads((out / "eval.json").read_text())["metrics"]
    assert metrics["all"]["called"] > 30


def test_accuracy_is_reasonable(called_run):
    _, out, _ = called_run
    metrics = json.loads((out / "eval.json").read_text())["metrics"]
    assert metrics["all"]["precision"] > 0.8, metrics
    assert metrics["all"]["recall"] > 0.8, metrics
    assert metrics["SNV"]["recall"] > 0.85, metrics


def test_pipeline_stages_exist(called_run):
    _, out, _ = called_run
    for name in ["aln.bam", "raw.vcf.gz", "norm.vcf.gz", "calls.vcf.gz",
                 "qc_flagstat.txt", "qc_coverage.txt"]:
        assert (out / name).exists(), f"missing {name}"


def test_report_written(called_run):
    _, out, _ = called_run
    html = (out / "report.html").read_text()
    assert "WGS variant calling report" in html
    assert "Precision" in html
    assert (out / "figures" / "qual_hist.png").exists()


def test_truth_set_is_well_formed():
    truth = (REPO / "data" / "sim" / "truth.vcf")
    if not truth.exists():
        pytest.skip("run 'make sim' first")
    lines = [l for l in truth.read_text().splitlines() if not l.startswith("#")]
    assert len(lines) > 100
    assert all(len(l.split("\t")) >= 5 for l in lines)
