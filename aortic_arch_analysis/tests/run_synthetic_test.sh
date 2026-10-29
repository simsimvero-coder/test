#!/usr/bin/env bash
# Runs the full pipeline on SYNTHETIC data (code testing only; not study results).
# Usage, from aortic_arch_analysis/:   bash tests/run_synthetic_test.sh
set -euo pipefail
Rscript tests/make_synthetic_data.R tests/synthetic/synthetic_arch.csv
export ARCH_RAW_DATA=tests/synthetic/synthetic_arch.csv
export ARCH_RESULTS=tests/synthetic_results
export ARCH_DERIVED=tests/synthetic_derived
export ARCH_M=${ARCH_M:-5}
export ARCH_NSIM=${ARCH_NSIM:-200}
export ARCH_GLMM_FALLBACK=${ARCH_GLMM_FALLBACK:-none}
Rscript run_all.R
