SHELL := /bin/bash
ENV_NAME := Morpheus

.PHONY: help install env env-mamba env-conda env-staged check test clean

help:
	@echo "Morpheus Makefile"
	@echo "  make install      Create the environment and verify the toolchain"
	@echo "  make env          Create/update the environment only (auto backend)"
	@echo "  make env-mamba    Create/update with mamba"
	@echo "  make env-conda    Create/update with conda"
	@echo "  make env-staged   Staged install, for memory-limited login nodes"
	@echo "  make check        Verify tools and run the self test; create nothing"
	@echo "  make test         Run the smoke test"
	@echo "  make clean        Remove test scratch and __pycache__"

install:
	bash install.sh

env:
	bash install.sh --env-only --backend=auto

env-mamba:
	bash install.sh --env-only --backend=mamba

env-conda:
	bash install.sh --env-only --backend=conda

env-staged:
	bash install.sh --env-only --backend=staged

check:
	bash install.sh --check-only

test:
	bash tests/smoke_test.sh

clean:
	rm -rf tests/tmp_smoke
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	find . -name '*.pyc' -delete
