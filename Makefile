.PHONY: ref sim call eval all test clean

PY ?= ../venv/bin/python
THREADS ?= 2

all: ref sim call eval

ref:
	bash scripts/get_reference.sh

sim:
	$(PY) scripts/simulate_reads.py --ref data/ref/ref.fa --outdir data/sim \
		--coverage 30 --snvs 400 --indels 60 --seed 7

call:
	THREADS=$(THREADS) bash scripts/call_variants.sh

eval:
	$(PY) scripts/evaluate.py --truth data/sim/truth.vcf --calls results/calls.vcf.gz \
		--outdir results --ref data/ref/ref.fa --coverage 30

test:
	$(PY) -m pytest tests -q

clean:
	rm -rf data results __pycache__ scripts/__pycache__ tests/__pycache__
