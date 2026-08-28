#!/usr/bin/env bash
###############################################################################
# Morpheus installer.
#
# Creates the `Morpheus` environment and checks that every external tool the
# pipeline shells out to is present. It never calls `module load` for you: on a
# managed cluster (Apocrita and similar) load your site's conda module first,
# e.g. `ml miniforge`, then run this.
#
#   bash install.sh                     # auto: mamba, then conda, then staged
#   bash install.sh --backend=mamba
#   bash install.sh --backend=conda
#   bash install.sh --backend=staged    # lower solver pressure, for tight nodes
#   bash install.sh --backend=current   # no env; check the active interpreter
#   bash install.sh --env-only          # create the environment, skip the checks
#   bash install.sh --check-only        # check tools, create nothing
#
# Morpheus itself needs no third-party Python packages, so a failed solve only
# costs you the external tools, never the pipeline logic.
###############################################################################
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="${MORPHEUS_ENV:-Morpheus}"
BACKEND="${MORPHEUS_INSTALL_BACKEND:-auto}"
ENV_ONLY=0
CHECK_ONLY=0

usage() { sed -n '3,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

for arg in "$@"; do
    case "$arg" in
        --backend=*)  BACKEND="${arg#--backend=}" ;;
        --env=*)      ENV_NAME="${arg#--env=}" ;;
        --env-only)   ENV_ONLY=1 ;;
        --check-only) CHECK_ONLY=1 ;;
        -h|--help)    usage; exit 0 ;;
        *) echo "[Morpheus] unknown option: $arg" >&2; usage >&2; exit 2 ;;
    esac
done

say()  { printf '[Morpheus] %s\n' "$*"; }
warn() { printf '[Morpheus] WARNING: %s\n' "$*" >&2; }
have() { command -v "$1" >/dev/null 2>&1; }

run() { say "\$ $*"; "$@"; }

env_exists_with() {
    local mgr="$1"
    "$mgr" env list 2>/dev/null | awk '{print $1}' | grep -qx "${ENV_NAME}"
}

# --------------------------------------------------------------------------
# environment creation
# --------------------------------------------------------------------------

try_full() {
    local mgr="$1"
    have "$mgr" || return 1
    if env_exists_with "$mgr"; then
        say "updating existing '${ENV_NAME}' with ${mgr}"
        run "$mgr" env update -n "${ENV_NAME}" -f "${HERE}/environment.yml" --prune && return 0
    else
        say "creating '${ENV_NAME}' with ${mgr}"
        run "$mgr" env create -n "${ENV_NAME}" -f "${HERE}/environment.yml" && return 0
    fi
    return 1
}

# One solve over the whole file can exhaust memory on a login node. Installing
# in groups keeps each solve small; the result is the same environment.
try_staged() {
    local mgr="$1"
    have "$mgr" || return 1
    say "staged install with ${mgr} (small solves, for memory-limited nodes)"
    env_exists_with "$mgr" || run "$mgr" create -y -n "${ENV_NAME}" \
        -c conda-forge -c bioconda python=3.11 pip git || return 1
    run "$mgr" install -y -n "${ENV_NAME}" -c conda-forge -c bioconda blast mafft || return 1
    run "$mgr" install -y -n "${ENV_NAME}" -c conda-forge -c bioconda hyphy bedtools || return 1
    run "$mgr" install -y -n "${ENV_NAME}" -c conda-forge r-base r-ape r-ggplot2 r-paletteer || return 1
    return 0
}

install_env() {
    case "${BACKEND}" in
        auto)
            try_full mamba && return 0
            warn "mamba solve failed or unavailable; trying conda"
            try_full conda && return 0
            warn "conda solve failed; trying a staged install"
            try_staged mamba || try_staged conda || return 1
            ;;
        mamba)   try_full mamba || try_staged mamba || return 1 ;;
        conda)   try_full conda || try_staged conda || return 1 ;;
        staged)  try_staged mamba || try_staged conda || return 1 ;;
        current) say "using the active interpreter; creating no environment" ;;
        *) echo "[Morpheus] invalid backend '${BACKEND}'" >&2; return 2 ;;
    esac
    return 0
}

# --------------------------------------------------------------------------
# verification
# --------------------------------------------------------------------------

check_tools() {
    local missing=()
    say "checking external tools"

    # Look inside the env we just built, without requiring `conda activate` to
    # work in a non-interactive shell.
    local envbin=""
    for mgr in mamba conda; do
        have "$mgr" || continue
        local base
        base="$("$mgr" info --base 2>/dev/null)" || continue
        [[ -d "${base}/envs/${ENV_NAME}/bin" ]] && { envbin="${base}/envs/${ENV_NAME}/bin"; break; }
    done
    [[ -n "${envbin}" ]] && export PATH="${envbin}:${PATH}"

    for exe in python3 blastx blastp makeblastdb mafft hyphy Rscript; do
        if have "$exe"; then
            printf '  %-14s %s\n' "$exe" "$(command -v "$exe")"
        else
            printf '  %-14s MISSING\n' "$exe"; missing+=("$exe")
        fi
    done

    if have Rscript; then
        local rmissing
        rmissing="$(Rscript -e 'p <- c("ape","ggplot2","paletteer"); cat(paste(p[!sapply(p, requireNamespace, quietly=TRUE)], collapse=" "))' 2>/dev/null)"
        if [[ -n "${rmissing}" ]]; then
            warn "R packages missing: ${rmissing}"
            warn "install with: Rscript -e 'install.packages(c(\"${rmissing// /\",\"}\"), repos=\"https://cloud.r-project.org\")'"
        else
            say "R packages present: ape, ggplot2, paletteer"
        fi
    fi

    # The pipeline logic is standard-library only; confirm the interpreter is
    # new enough rather than checking for packages that are not needed.
    if have python3; then
        python3 - <<'EOF' || missing+=("python>=3.9")
import sys
if sys.version_info < (3, 9):
    print(f"  python3        TOO OLD ({sys.version.split()[0]}); need >= 3.9")
    raise SystemExit(1)
print(f"  python3        {sys.version.split()[0]} (no third-party packages needed)")
EOF
    fi

    if [[ ${#missing[@]} -gt 0 ]]; then
        warn "missing: ${missing[*]}"
        return 1
    fi
    say "all external tools found"
    return 0
}

self_test() {
    say "running the built-in self test"
    PYTHONPATH="${HERE}${PYTHONPATH:+:${PYTHONPATH}}" python3 - <<'EOF'
import itertools, random, sys
from morpheus.matching import max_weight_matching
from morpheus.common import translate, revcomp, orf_report
from morpheus import newick

random.seed(0)
for _ in range(60):
    r, c = random.randint(1, 5), random.randint(1, 5)
    P = [[round(random.random(), 3) for _ in range(c)] for _ in range(r)]
    got = sum(P[i][j] for i, j in max_weight_matching(P, min_profit=-1))
    best = max(sum(P[a][b] for a, b in zip(rows, cols))
               for rows in itertools.permutations(range(r), min(r, c))
               for cols in itertools.permutations(range(c), min(r, c)))
    assert abs(got - best) < 1e-9, "matching is not optimal"

assert translate("ATGGTGCATCTGACTCCTGAGGAGAAGTAA") == "MVHLTPEEK*"
assert revcomp("ATGCN") == "NGCAT"
assert orf_report("ATG" + "AAA" * 750 + "TGA" + "AAAA")["orf_fraction"] > 0.95

t = newick.parse("((A:1,B:2)n1:3,(C:4,D:5)n2:6)root;")
assert sorted(t.leaf_names()) == ["A", "B", "C", "D"]
assert newick.write(newick.prune(t, {"A", "C"})) == "(A:4,C:10)root;"
print("  optimal matching, codon table, ORF report and Newick pruning all pass")
EOF
}

# --------------------------------------------------------------------------
say "Morpheus $(cat "${HERE}/VERSION" 2>/dev/null || echo '?')  home=${HERE}"

if [[ ${CHECK_ONLY} -eq 1 ]]; then
    check_tools; rc=$?; self_test || rc=1; exit $rc
fi

install_env || { warn "environment creation did not complete"; }

if [[ ${ENV_ONLY} -eq 1 ]]; then
    say "environment step finished (--env-only)"
    exit 0
fi

rc=0
check_tools || rc=1
self_test || rc=1

chmod +x "${HERE}/bin/morpheus" 2>/dev/null || true

cat <<EOF

[Morpheus] Installation finished (backend: ${BACKEND})

Activate and run:
  conda activate ${ENV_NAME}
  export PATH="${HERE}/bin:\$PATH"
  morpheus env          # confirm what was resolved
  morpheus run --dry-run

On a cluster, load your site's conda module first, e.g.:
  ml miniforge && conda activate ${ENV_NAME}
EOF
exit $rc
