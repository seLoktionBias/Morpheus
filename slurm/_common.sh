#!/usr/bin/env bash
# Shared preamble for the Morpheus Slurm scripts.
#
# The one rule here: never guess, and never swallow a failure. An environment
# that half-activated produces a job that dies several steps later with an error
# about a missing file, and the real cause -- `conda activate` doing nothing in a
# non-interactive shell -- is nowhere in the log.
#
# Environment modes, chosen with --env_mode:
#   inherit  (default)  use the PATH of the shell that submitted the job.
#                       sbatch --export=ALL carries it in. If you ran
#                       `mamba activate Morpheus` before submitting, the job
#                       already has everything.
#   conda|mamba         run the work through `<manager> run -n <name>`, which
#                       works in a non-interactive shell; `conda activate` does
#                       not, and is the usual reason a batch job fails where the
#                       same command worked interactively.
#   venv                source a Python venv's activate script.
#
# --slurm_module is opt-in and never assumed: sites name their modules
# differently, and loading the wrong one silently is worse than not loading one.

MORPHEUS_HOME="${MORPHEUS_HOME:?set MORPHEUS_HOME to the Morpheus checkout}"
PROJECT_DIR="${MORPHEUS_PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "${PROJECT_DIR}"
mkdir -p logs

_fail() { echo "[Morpheus] ERROR: $*" >&2; exit 1; }

# ---- site module (only if the submitter asked for one) --------------------
if [[ -n "${MORPHEUS_MODULE:-}" ]]; then
    if type module >/dev/null 2>&1; then
        module load "${MORPHEUS_MODULE}" \
            || _fail "module load ${MORPHEUS_MODULE} failed"
    elif command -v ml >/dev/null 2>&1; then
        ml "${MORPHEUS_MODULE}" || _fail "ml ${MORPHEUS_MODULE} failed"
    else
        _fail "--slurm_module ${MORPHEUS_MODULE} was requested but neither 'module' nor 'ml' is available"
    fi
fi

# ---- environment ----------------------------------------------------------
MORPHEUS_ENV_MODE="${MORPHEUS_ENV_MODE:-inherit}"
MORPHEUS_ENV_NAME="${MORPHEUS_ENV_NAME:-${MORPHEUS_ENV:-Morpheus}}"
MORPHEUS_VENV_PATH="${MORPHEUS_VENV_PATH:-}"

declare -a ENV_PREFIX=()
case "${MORPHEUS_ENV_MODE}" in
    inherit)
        ;;
    conda|mamba)
        command -v "${MORPHEUS_ENV_MODE}" >/dev/null 2>&1 \
            || _fail "--env_mode ${MORPHEUS_ENV_MODE} was requested but '${MORPHEUS_ENV_MODE}' is not on PATH in this job. Add --slurm_module <your conda module>, or use --env_mode inherit after activating the environment before you submit."
        # `<manager> run -n` needs no interactive shell, unlike `activate`.
        ENV_PREFIX=("${MORPHEUS_ENV_MODE}" run -n "${MORPHEUS_ENV_NAME}")
        # Without this, `conda run` buffers everything until the process exits:
        # a twelve-hour job would show an empty log the whole time it ran, and
        # nothing at all if it were killed. Older conda lacks the flag, so probe.
        if "${MORPHEUS_ENV_MODE}" run --help 2>&1 | grep -q -- --no-capture-output; then
            ENV_PREFIX+=(--no-capture-output)
        else
            echo "[Morpheus] note: ${MORPHEUS_ENV_MODE} has no --no-capture-output;" \
                 "job output will only appear when each step finishes" >&2
        fi
        ;;
    venv)
        [[ -n "${MORPHEUS_VENV_PATH}" ]] || _fail "--env_mode venv needs --venv_path"
        [[ -f "${MORPHEUS_VENV_PATH}/bin/activate" ]] \
            || _fail "no activation script at ${MORPHEUS_VENV_PATH}/bin/activate"
        # shellcheck disable=SC1091
        source "${MORPHEUS_VENV_PATH}/bin/activate"
        ;;
    *)
        _fail "unknown --env_mode '${MORPHEUS_ENV_MODE}'" ;;
esac

# Fail here, with a message naming the cause, rather than several steps later
# with one that does not. Under `<manager> run` the tools are inside the
# environment, so probe there rather than on this shell's PATH.
_need=(python3 mafft blastx)
[[ "${MORPHEUS_SKIP_PLOT:-0}" == "1" ]] || _need+=(Rscript)
_missing=()
for _exe in "${_need[@]}"; do
    "${ENV_PREFIX[@]+"${ENV_PREFIX[@]}"}" command -v "${_exe}" >/dev/null 2>&1 \
        || _missing+=("${_exe}")
done
if [[ ${#_missing[@]} -gt 0 ]]; then
    echo "[Morpheus] ERROR: not on PATH in this job: ${_missing[*]}" >&2
    echo "  env_mode  ${MORPHEUS_ENV_MODE}" >&2
    echo "  env_name  ${MORPHEUS_ENV_NAME}" >&2
    echo "  module    ${MORPHEUS_MODULE:-<none requested>}" >&2
    echo >&2
    echo "  Either activate the environment before submitting and keep the" >&2
    echo "  default --env_mode inherit, or submit with:" >&2
    echo "    --env_mode mamba --env_name ${MORPHEUS_ENV_NAME} --slurm_module <your conda module>" >&2
    exit 1
fi
unset _need _missing _exe

# shellcheck disable=SC1091
source "${MORPHEUS_HOME}/config/config.sh"
export PYTHONPATH="${MORPHEUS_HOME}${PYTHONPATH:+:${PYTHONPATH}}"
export THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-4}}"
# config.sh's own conda search would fight the mode chosen above.
export MORPHEUS_ENV_RESOLVED=1

# Carry --search and --skip_plot into every job, so one submitted chain runs the
# same analysis end to end. A step the flags exclude is skipped by run_pipeline
# itself, which keeps that decision in exactly one place.
declare -a MORPHEUS_FLAGS=(--search "${MORPHEUS_SEARCH:-both}")
[[ "${MORPHEUS_SKIP_PLOT:-0}" == "1" ]] && MORPHEUS_FLAGS+=(--skip_plot)

# Run anything inside the chosen environment.
morpheus_exec() { "${ENV_PREFIX[@]+"${ENV_PREFIX[@]}"}" "$@"; }

morpheus_run() {
    morpheus_exec bash "${MORPHEUS_HOME}/scripts/run_pipeline.sh" \
        "${MORPHEUS_FLAGS[@]}" "$@"
}
