#!/usr/bin/env bash
# Shared preamble for the Morpheus Slurm scripts: load the site's conda module,
# activate the environment, move to the project directory, load the config.
#
# Apocrita provides conda through `ml miniforge`. Set MORPHEUS_MODULE to your
# site's module name, or to the empty string if conda is already on PATH.
#
# These scripts are normally submitted for you by `morpheus run --mode slurm`,
# which passes every resolved setting through --export. Submitting them by hand
# works too; then they read paths.txt from the submission directory as usual.

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

# Carry --search and --skip_plot into every job, so one submitted chain runs the
# same analysis end to end. A step the flags exclude is skipped by run_pipeline
# itself, which keeps the decision in exactly one place.
declare -a MORPHEUS_FLAGS=(--search "${MORPHEUS_SEARCH:-both}")
[[ "${MORPHEUS_SKIP_PLOT:-0}" == "1" ]] && MORPHEUS_FLAGS+=(--skip_plot)

morpheus_run() {
    bash "${MORPHEUS_HOME}/scripts/run_pipeline.sh" "${MORPHEUS_FLAGS[@]}" "$@"
}
