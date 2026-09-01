#!/usr/bin/env bash
# slurm/_common.sh must fail loudly, at the top of the log, naming the cause.
#
# The failure this guards against: `conda activate` does nothing useful in a
# non-interactive batch shell, and the old preamble ran it with `|| true`. The
# job then died several steps later complaining about a missing file, and the
# real cause appeared nowhere. A dependent chain turns that into
# `DependencyNeverSatisfied` on every downstream job, with no clue why.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d)" || { echo "mktemp failed" >&2; exit 1; }
# Guard the cleanup and never let it decide the exit status: a trap whose last
# command fails silently overwrites the result of the tests themselves.
cleanup() { [[ -n "${TMP:-}" && -d "${TMP}" ]] && rm -rf "${TMP}"; return 0; }
trap cleanup EXIT
pass=0; fail=0

# Run _common.sh in a clean environment, so the host's PATH cannot mask a bug.
sandbox() {
    env -i PATH=/usr/bin:/bin HOME="${HOME}" \
        MORPHEUS_HOME="${HERE}" MORPHEUS_PROJECT_DIR="${TMP}" "$@" \
        bash -c 'set -euo pipefail; source "${MORPHEUS_HOME}/slurm/_common.sh"' 2>&1
}

# Accepts several alternative phrases: the contract under test is "it fails
# loudly, naming something the user can act on", not "it emits this exact
# string". A job whose environment is missing in two ways can legitimately
# report either. On failure it dumps the sandbox's own view of PATH, because a
# rare disagreement here has meant the sandbox was not as isolated as intended.
want_fail() {  # label, phrase [, phrase ...], -- , extra env...
    local label="$1"; shift
    local patterns=()
    while [[ $# -gt 0 && "$1" != "--" ]]; do patterns+=("$1"); shift; done
    [[ "${1-}" == "--" ]] && shift
    local out rc
    out="$(sandbox "$@")"; rc=$?
    if [[ ${rc} -eq 0 ]]; then
        echo "  FAIL  ${label}: exited 0, should have failed"; fail=$((fail+1))
        return
    fi
    local p
    for p in "${patterns[@]}"; do
        if grep -q -- "${p}" <<<"${out}"; then
            echo "  PASS  ${label}"; pass=$((pass+1)); return
        fi
    done
    echo "  FAIL  ${label}: message matched none of: ${patterns[*]}"
    echo "        got : $(head -1 <<<"${out}")"
    echo "        PATH: $(env -i PATH=/usr/bin:/bin bash -c "echo \$PATH")"
    echo "        conda in sandbox: $(env -i PATH=/usr/bin:/bin bash -c \
                     'command -v conda || echo ABSENT')"
    fail=$((fail+1))
}

echo "== _common.sh failure modes =="
want_fail "missing tools are named"          "not on PATH in this job" --
want_fail "and the remedy is spelled out"    "env_mode inherit" --
want_fail "unknown env_mode is refused"      "unknown --env_mode" -- MORPHEUS_ENV_MODE=bogus
# Either message is a correct outcome: conda missing, or conda present but the
# environment it names not usable. Both name the cause and both exit non-zero.
want_fail "conda absent from the job" \
    "'conda' is not on PATH in this job" "not on PATH in this job" -- MORPHEUS_ENV_MODE=conda
want_fail "mamba absent from the job" \
    "'mamba' is not on PATH in this job" "not on PATH in this job" -- MORPHEUS_ENV_MODE=mamba
want_fail "venv without a path"              "needs --venv_path" -- MORPHEUS_ENV_MODE=venv
want_fail "venv path that does not exist"    "no activation script" -- MORPHEUS_ENV_MODE=venv MORPHEUS_VENV_PATH=/nonexistent
want_fail "module requested, none available" "neither 'module' nor 'ml'" -- MORPHEUS_MODULE=miniforge

# MORPHEUS_HOME is the one thing that cannot be defaulted.
out="$(env -i PATH=/usr/bin:/bin HOME="${HOME}" bash -c \
      'set -euo pipefail; source "'"${HERE}"'/slurm/_common.sh"' 2>&1)"; rc=$?
if [[ ${rc} -ne 0 ]] && grep -q "MORPHEUS_HOME" <<<"${out}"; then
    echo "  PASS  unset MORPHEUS_HOME is refused"; pass=$((pass+1))
else
    echo "  FAIL  unset MORPHEUS_HOME is refused"; fail=$((fail+1))
fi

echo "== ${pass} passed, ${fail} failed =="
[[ ${fail} -eq 0 ]] && exit 0
exit 1
