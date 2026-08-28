#!/usr/bin/env bash
# Shared preamble for the Morpheus SLURM scripts: load the site's conda module,
# activate the environment, move to the project directory, load the config.
#
# Apocrita provides conda through `ml miniforge`. Set MORPHEUS_MODULE to your
# site's module name, or to the empty string if conda is already on PATH.

MODULE="${MORPHEUS_MODULE-miniforge}"
if [[ -n "${MODULE}" ]] && command -v module >/dev/null 2>&1; then
    module load "${MODULE}" || true
fi

if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${MORPHEUS_ENV:-Morpheus}" || true
fi

MORPHEUS_HOME="${MORPHEUS_HOME:?set MORPHEUS_HOME to the Morpheus checkout}"
PROJECT_DIR="${MORPHEUS_PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "${PROJECT_DIR}"
mkdir -p logs

# shellcheck disable=SC1091
source "${MORPHEUS_HOME}/config/config.sh"
export PYTHONPATH="${MORPHEUS_HOME}${PYTHONPATH:+:${PYTHONPATH}}"
export THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-4}}"

morpheus_run() { bash "${MORPHEUS_HOME}/scripts/run_pipeline.sh" "$@"; }
