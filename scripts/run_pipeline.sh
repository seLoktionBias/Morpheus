#!/usr/bin/env bash
###############################################################################
# run_pipeline.sh - run the whole analysis, or any part of it.
#
#   morpheus run --gene MX1
#   morpheus run --gene_list genes.txt --search both
#   morpheus run --path_file paths.txt --mode slurm
#   morpheus run --from assign                  # resume from a step
#   morpheus run --only copy-number             # one step
#
# Steps, in order:
#   cds-table  human-reference  families  bat-search  screen  assign
#   compare  copy-number  deliverables  align  plots  summary
#
# Every step writes into results/ and reads only what earlier steps produced,
# so re-running one step never invalidates the others.
###############################################################################
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MORPHEUS_HOME="${MORPHEUS_HOME:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
export MORPHEUS_HOME

ALL_STEPS=(cds-table human-reference families bat-search screen assign
           compare copy-number deliverables align plots summary)

FROM=""; SPECIES=(); OVERWRITE=""; DRY_RUN=0
GENE_ARG=""; GENE_LIST_ARG=""; MODE="local"; SEARCH="both"
TREE_ARG=""; BAT_ARG=""; REF_ARG=""; OUT_ARG=""; PATH_FILE_ARG=""
SKIP_PLOT=0; SLURM_MODULE_SET=0
SLURM_PARTITION=""; SLURM_ACCOUNT=""; SLURM_MODULE_ARG=""; SLURM_EXTRA=""
ARRAY_THROTTLE=""; SLURM_TIME=""
ENV_MODE="inherit"; ENV_NAME=""; VENV_PATH=""
PER_GENE="yes"; SINGLE_GENE=""

usage() {
    cat <<'EOF'
morpheus run - region-restricted transcript recovery, orthologous transcript
               assignment and gene copy number across TOGA2 genomes.

  morpheus run --gene_list gene_list.txt        # the usual invocation
  morpheus run --gene_list gene_list.txt --dry-run     # check first

Inputs come from paths.txt in the current directory unless you say otherwise.
Full examples at the bottom of this message.

Genes (give one of these):
  --gene NAME            a single gene symbol
  --gene_list FILE       a text file, one gene symbol per line

Inputs (give --path_file, or the directories directly):
  --path_file FILE       paths.txt holding the directories (see README)
  --bat_dir DIR          TOGA2 annotation directory, one sub-directory per genome
  --ref_dir DIR          human genome + Ensembl annotation directory
  --output DIR           working/output directory (default: current directory)
  --tree FILE            species tree in Newick, full path including filename

How much to run:
  --mode local|slurm     run here, or submit the whole chain to Slurm
                           (default local)
  --search region|similarity|both
                         region      - the gene's own syntenic locus, ranked by
                                       synteny and exon structure
                         similarity  - the gene anywhere, ranked by sequence
                                       identity
                         both        - both, side by side (default)
  --skip_plot            produce no figures, only the tables
  --per_gene yes|no      yes (default): each gene is run on its own -- one
                           Slurm job per gene, or one gene at a time locally,
                           so a long list cannot swamp the machine and one
                           gene failing cannot spoil the rest. Results land in
                           <output>/genes/<GENE>/ and are merged afterwards
                           into <output>/combined/.
                         no: one pass over the whole list, sharing the search
                           across genes. Faster for many genes, because each
                           query genome is read once instead of once per gene.

Environment (how a job finds python3, blast, mafft and R):
  --env_mode MODE        inherit | conda | mamba | venv   (default inherit)
                           inherit     use the PATH of the submitting shell.
                                       Activate the environment first, then
                                       submit; --export=ALL carries it in.
                           conda|mamba run the work through `<mgr> run -n NAME`,
                                       which works in a non-interactive shell.
                                       `conda activate` does not, and that is
                                       the usual reason a batch job fails where
                                       the same command worked interactively.
                           venv        source --venv_path/bin/activate
  --env_name NAME        environment for --env_mode conda|mamba (default Morpheus)
  --venv_path DIR        virtualenv for --env_mode venv

Slurm (only meaningful with --mode slurm):
  --time SPEC            wall clock for every job: hours (12), HH:MM:SS or
                           D-HH:MM:SS. Omit to keep each job's own default,
                           which differs by job (see slurm/*.sbatch)
  --slurm_partition NAME queue to submit to, e.g. compute
  --slurm_account NAME   account to charge
  --slurm_module NAME    site module providing conda (default miniforge;
                           pass "" if conda is already on PATH)
  --array_throttle N     at most N search tasks running at once
  --slurm_extra "..."    passed verbatim to every sbatch, e.g. "--qos=short".
                           Per-job --mem and --time defaults live in the sbatch
                           files; anything given here overrides all four jobs.

Partial runs:
  --from STEP            start at STEP and run everything after it
  --only STEP            run just STEP (repeatable)
  --species NAME...      restrict the search to these Genus_species
  --overwrite            redo per-species searches that are already complete
  --dry-run              print the plan, then stop
  -h, --help             this message

Steps: cds-table human-reference families bat-search screen assign compare
       copy-number deliverables align plots summary

Settings (export before running): THREADS, MIN_SPECIES_FRACTION,
       MIN_ASSIGNMENT_SCORE, MIN_ASSIGNMENT_PIDENT, PLOT_MIN_SPECIES
       (PLOT_MIN_SPECIES defaults to MIN_SPECIES_FRACTION of the species searched)

EXAMPLES

  A paths.txt saves repeating the directories. Three lines is the whole file:

      human_genome_dir=/data/reference/human
      bat_annotation_dir=/data/toga2/bat1k
      primary_working_dir=/scratch/me/isg_project

  Run from the directory holding it:

      morpheus run --gene_list gene_list.txt

  One gene, no list file needed:

      morpheus run --gene OAS1

  No paths.txt at all - everything on the command line:

      morpheus run --gene_list gene_list.txt \
          --ref_dir /data/reference/human \
          --bat_dir /data/toga2/bat1k \
          --tree    /data/reference/bat1k_tree.nwk \
          --output  /scratch/me/isg_project

  Shared paths.txt, but write this run somewhere else:

      morpheus run --gene_list gene_list.txt \
          --path_file /shared/paths.txt --output ./run2

  On a cluster. Submits four dependent jobs and returns immediately; pass the
  whole gene list, not one gene at a time - the search is shared across genes
  and the array already parallelises over genomes:

      MORPHEUS=/data/home/me/software/Morpheus/bin/morpheus
      "$MORPHEUS" run \
          --gene_list "${wdr}/gene_list.txt" \
          --path_file "${wdr}/paths.txt" \
          --mode slurm \
          --slurm_partition compute \
          --slurm_module miniforge \
          --array_throttle 20

  Half the work: the syntenic locus only, tables but no figures:

      morpheus run --gene_list gene_list.txt --search region --skip_plot

  Resume after a failure, or redo one step:

      morpheus run --gene_list gene_list.txt --from assign
      morpheus run --only copy-number

  Try two species before committing to all of them:

      morpheus run --gene OAS1 --species Myotis_myotis Phyllostomus_discolor

WHEN IT WILL NOT START

  morpheus env      prints every resolved path with ok / MISSING beside it.
                    That is almost always the answer.
EOF
}

declare -a ONLY_STEPS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --gene)        GENE_ARG="${2:?--gene needs a gene symbol}"; shift 2 ;;
        --gene_list|--gene-list)
                       GENE_LIST_ARG="${2:?--gene_list needs a file}"; shift 2 ;;
        --mode)        MODE="${2:?--mode needs local or slurm}"; shift 2 ;;
        --search)      SEARCH="${2:?--search needs region, similarity or both}"; shift 2 ;;
        --tree)        TREE_ARG="${2:?--tree needs a file}"; shift 2 ;;
        --bat_dir|--bat-dir)
                       BAT_ARG="${2:?--bat_dir needs a directory}"; shift 2 ;;
        --ref_dir|--ref-dir)
                       REF_ARG="${2:?--ref_dir needs a directory}"; shift 2 ;;
        --output)      OUT_ARG="${2:?--output needs a directory}"; shift 2 ;;
        --path_file|--path-file)
                       PATH_FILE_ARG="${2:?--path_file needs a file}"; shift 2 ;;
        --skip_plot|--skip-plot) SKIP_PLOT=1; shift ;;
        --per_gene|--per-gene)
                       PER_GENE="${2:?--per_gene needs yes or no}"; shift 2 ;;
        # internal: the per-gene driver re-enters this script for one gene
        --_single_gene) SINGLE_GENE="${2:?}"; PER_GENE="no"; shift 2 ;;
        --env_mode|--env-mode)
                       ENV_MODE="${2:?--env_mode needs inherit, conda, mamba or venv}"; shift 2 ;;
        --env_name|--env-name)
                       ENV_NAME="${2:?--env_name needs a name}"; shift 2 ;;
        --venv_path|--venv-path)
                       VENV_PATH="${2:?--venv_path needs a directory}"; shift 2 ;;
        --time)        SLURM_TIME="${2:?--time needs hours, HH:MM:SS or D-HH:MM:SS}"; shift 2 ;;
        --slurm_partition|--slurm-partition)
                       SLURM_PARTITION="${2:?--slurm_partition needs a name}"; shift 2 ;;
        --slurm_account|--slurm-account)
                       SLURM_ACCOUNT="${2:?--slurm_account needs a name}"; shift 2 ;;
        --slurm_module|--slurm-module)
                       SLURM_MODULE_ARG="${2-}"; SLURM_MODULE_SET=1; shift 2 ;;
        --array_throttle|--array-throttle)
                       ARRAY_THROTTLE="${2:?--array_throttle needs a number}"; shift 2 ;;
        --slurm_extra|--slurm-extra)
                       SLURM_EXTRA="${2:?--slurm_extra needs a string}"; shift 2 ;;
        --from)        FROM="${2:?--from needs a step}"; shift 2 ;;
        --only)        ONLY_STEPS+=("${2:?--only needs a step}"); shift 2 ;;
        --species)     shift; while [[ $# -gt 0 && "$1" != --* ]]; do SPECIES+=("$1"); shift; done ;;
        --overwrite)   OVERWRITE="--overwrite"; shift ;;
        --dry-run)     DRY_RUN=1; shift ;;
        -h|--help)     usage; exit 0 ;;
        *) echo "ERROR: unknown option '$1'" >&2; echo "Try: morpheus run --help" >&2; exit 1 ;;
    esac
done

die() { echo "ERROR: $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# validate the flags before anything touches the filesystem
# ---------------------------------------------------------------------------
case "${MODE}" in local|slurm) ;; *) die "--mode must be 'local' or 'slurm', not '${MODE}'" ;; esac
case "${SEARCH}" in region|similarity|both) ;;
    *) die "--search must be 'region', 'similarity' or 'both', not '${SEARCH}'" ;;
esac

[[ -n "${GENE_ARG}" && -n "${GENE_LIST_ARG}" ]] && \
    die "give --gene or --gene_list, not both"

if [[ "${MODE}" != "slurm" ]]; then
    warn_unused() { [[ -n "$2" ]] && echo "WARNING: $1 is ignored without --mode slurm" >&2; return 0; }
    warn_unused --slurm_partition "${SLURM_PARTITION}"
    warn_unused --slurm_account   "${SLURM_ACCOUNT}"
    warn_unused --array_throttle  "${ARRAY_THROTTLE}"
    warn_unused --slurm_extra     "${SLURM_EXTRA}"
    warn_unused --time            "${SLURM_TIME}"
    [[ ${SLURM_MODULE_SET} -eq 1 ]] && \
        echo "WARNING: --slurm_module is ignored without --mode slurm" >&2
fi
[[ -n "${ARRAY_THROTTLE}" && ! "${ARRAY_THROTTLE}" =~ ^[0-9]+$ ]] && \
    die "--array_throttle must be a number, not '${ARRAY_THROTTLE}'"

case "${PER_GENE}" in yes|no) ;;
    *) die "--per_gene must be 'yes' or 'no', not '${PER_GENE}'" ;;
esac
case "${ENV_MODE}" in
    inherit|conda|mamba|venv) ;;
    *) die "--env_mode must be inherit, conda, mamba or venv, not '${ENV_MODE}'" ;;
esac
[[ "${ENV_MODE}" == "venv" && -z "${VENV_PATH}" ]] && die "--env_mode venv needs --venv_path"
[[ -n "${VENV_PATH}" && ! -f "${VENV_PATH}/bin/activate" ]] && \
    die "no activation script at ${VENV_PATH}/bin/activate"
[[ -n "${ENV_NAME}" && "${ENV_MODE}" == "inherit" ]] && \
    echo "WARNING: --env_name is ignored with --env_mode inherit" >&2
: "${ENV_NAME:=${MORPHEUS_ENV:-Morpheus}}"

# Accept plain hours (Pensieve's spelling), HH:MM:SS, or D-HH:MM:SS.
normalise_time() {
    local spec="$1"
    if [[ "${spec}" =~ ^[0-9]+$ ]]; then
        local d=$(( spec / 24 )) h=$(( spec % 24 ))
        (( d > 0 )) && printf '%d-%02d:00:00' "${d}" "${h}" || printf '%02d:00:00' "${h}"
        return 0
    fi
    [[ "${spec}" =~ ^([0-9]+-)?[0-9]{1,2}:[0-9]{2}:[0-9]{2}$ ]] || return 1
    printf '%s' "${spec}"
}
if [[ -n "${SLURM_TIME}" ]]; then
    SLURM_TIME="$(normalise_time "${SLURM_TIME}")" \
        || die "--time must be hours (12), HH:MM:SS or D-HH:MM:SS"
fi

for s in ${ONLY_STEPS[@]+"${ONLY_STEPS[@]}"} ${FROM:+"${FROM}"}; do
    ok=0; for k in "${ALL_STEPS[@]}"; do [[ "$s" == "$k" ]] && ok=1; done
    [[ $ok -eq 1 ]] || die "unknown step '$s'. Steps: ${ALL_STEPS[*]}"
done

# ---------------------------------------------------------------------------
# flags win over paths.txt, which wins over the built-in defaults. config.sh
# uses `: "${VAR:=default}"` throughout, so exporting here is all it takes.
# ---------------------------------------------------------------------------
if [[ -n "${PATH_FILE_ARG}" ]]; then
    [[ -s "${PATH_FILE_ARG}" ]] || die "--path_file not found or empty: ${PATH_FILE_ARG}"
    export PATHS_FILE="$(cd "$(dirname "${PATH_FILE_ARG}")" && pwd)/$(basename "${PATH_FILE_ARG}")"
    # A paths.txt elsewhere means the project lives there too, so gene_list.txt
    # and the tree resolve next to it rather than next to wherever you happen
    # to be standing. --output still decides where results are written.
    : "${MORPHEUS_PROJECT_DIR:=$(dirname "${PATHS_FILE}")}"
    export MORPHEUS_PROJECT_DIR
fi
[[ -n "${BAT_ARG}"  ]] && export BAT_ANNOTATION_DIR="${BAT_ARG}"
[[ -n "${REF_ARG}"  ]] && export HUMAN_GENOME_DIR="${REF_ARG}"
[[ -n "${OUT_ARG}"  ]] && export WORKING_DIR="${OUT_ARG}"
[[ -n "${TREE_ARG}" ]] && export SPECIES_TREE="${TREE_ARG}"
[[ -n "${GENE_LIST_ARG}" ]] && export GENE_LIST="${GENE_LIST_ARG}"

CACHE_DIR_FROM_CALLER="${CACHE_DIR:-}"
source "${MORPHEUS_HOME}/config/config.sh"

# --gene: a one-line gene list, written where the results live so a later
# --only step finds the same list without the flag being repeated.
if [[ -n "${GENE_ARG}" ]]; then
    mkdir -p "${RESULTS_DIR}"
    GENE_LIST="${RESULTS_DIR}/gene_list.single.txt"
    printf '%s\n' "${GENE_ARG}" > "${GENE_LIST}"
fi

# ---------------------------------------------------------------------------
# --search selects both the scopes to assign under and the rankings to build.
#
# The two are separate axes, and the flag deliberately couples them:
#   region      the gene's own syntenic locus, judged by position and exon
#               structure - the conservative orthology answer
#   similarity  the gene wherever it sits, judged by sequence identity - the
#               answer that finds the functional copy when the locus diverged
# SEQUENCE_SCOPE is the scope the FASTAs are drawn from; STATUS_SCOPES are the
# scopes that get a transcript-status table and figure.
# ---------------------------------------------------------------------------
case "${SEARCH}" in
    region)
        SCOPES="region_restricted";  POLICIES="synteny_aware"
        SEQUENCE_SCOPE="region_restricted"; STATUS_SCOPES="region_restricted"
        DEFAULT_POLICY="synteny_aware" ;;
    similarity)
        SCOPES="unrestricted";       POLICIES="sequence_similarity"
        SEQUENCE_SCOPE="unrestricted";      STATUS_SCOPES="unrestricted"
        DEFAULT_POLICY="sequence_similarity" ;;
    both)
        SCOPES="region_restricted unrestricted"
        POLICIES="sequence_similarity synteny_aware"
        SEQUENCE_SCOPE="unrestricted"
        STATUS_SCOPES="region_restricted unrestricted" ;;
esac
# The scope comparison needs both scopes; there is nothing to compare otherwise.
DO_COMPARE=0; [[ "${SEARCH}" == "both" ]] && DO_COMPARE=1

should_run() {
    local step="$1"
    [[ "$step" == "plots"   && ${SKIP_PLOT} -eq 1 ]] && return 1
    [[ "$step" == "compare" && ${DO_COMPARE} -eq 0 ]] && return 1
    if [[ ${#ONLY_STEPS[@]-0} -gt 0 ]]; then
        for s in "${ONLY_STEPS[@]}"; do [[ "$s" == "$step" ]] && return 0; done
        return 1
    fi
    if [[ -n "${FROM}" ]]; then
        local seen=0
        for s in "${ALL_STEPS[@]}"; do
            [[ "$s" == "${FROM}" ]] && seen=1
            [[ "$s" == "$step" ]] && { [[ $seen -eq 1 ]] && return 0 || return 1; }
        done
        return 1
    fi
    return 0
}

banner() {
    printf '\n%s\n== %-58s ==\n%s\n' \
        "==============================================================" "$1" \
        "=============================================================="
}

# ---------------------------------------------------------------------------
# --per_gene yes: each gene gets its own complete, self-contained results tree.
#
# Isolation is the point. Nothing is shared but the flattened GTF cache, which
# is read-only after the reference job writes it, so concurrent gene jobs cannot
# corrupt a common table and nothing needs locking. One gene failing leaves
# every other gene's results intact and re-runnable on its own.
#
# The cost is real and worth stating: each gene job reads every query genome's
# annotation, so N genes read them N times. --per_gene no runs one pass over the
# whole list instead, reading each genome once.
# ---------------------------------------------------------------------------
GENES_ROOT="${WORKING_DIR}/genes"
COMBINED_DIR="${WORKING_DIR}/combined"
# One flattened-GTF cache for every gene: it is the single most expensive thing
# in the pipeline and it does not depend on which gene is being analysed. It
# sits beside genes/ and combined/, not inside any one gene's results.
if [[ "${PER_GENE}" == "yes" ]]; then
    export CACHE_DIR="${CACHE_DIR_FROM_CALLER:-${WORKING_DIR}/cache}"
fi

read_genes() { grep -vE '^[[:space:]]*(#|$)' "${GENE_LIST}" | tr -d '\r' | awk '{print $1}'; }

if [[ "${PER_GENE}" == "yes" ]]; then
    [[ -s "${GENE_LIST}" ]] || die "missing gene list: ${GENE_LIST} (use --gene or --gene_list)"
    # A while-read loop, not mapfile: mapfile is bash 4+, and macOS still
    # ships bash 3.2 as /bin/bash. Local runs have to work there.
    GENES=()
    while IFS= read -r _g; do [[ -n "${_g}" ]] && GENES+=("${_g}"); done < <(read_genes)
    [[ ${#GENES[@]} -gt 0 ]] || die "no gene symbols in ${GENE_LIST}"

    # Flags every per-gene invocation inherits. Paths are already resolved, so a
    # gene job never re-reads paths.txt from a directory it is not standing in.
    declare -a PASS=(--search "${SEARCH}"
                     --bat_dir "${BAT_ANNOTATION_DIR}" --ref_dir "${HUMAN_GENOME_DIR}"
                     --tree "${SPECIES_TREE}")
    [[ ${SKIP_PLOT} -eq 1 ]] && PASS+=(--skip_plot)
    [[ -n "${OVERWRITE}" ]] && PASS+=(--overwrite)
    [[ -n "${FROM}" ]] && PASS+=(--from "${FROM}")
    for s in ${ONLY_STEPS[@]+"${ONLY_STEPS[@]}"}; do PASS+=(--only "$s"); done
    [[ ${#SPECIES[@]-0} -gt 0 ]] && PASS+=(--species "${SPECIES[@]}")

    if [[ "${MODE}" == "slurm" ]]; then
        [[ ${DRY_RUN} -eq 1 ]] || command -v sbatch >/dev/null 2>&1 \
            || die "--mode slurm, but 'sbatch' is not on PATH"
        mkdir -p "${WORKING_DIR}/logs" "${GENES_ROOT}"
        printf '%s\n' "${GENES[@]}" > "${GENES_ROOT}/.gene_order.txt"

        SB=()
        [[ -n "${SLURM_PARTITION}" ]] && SB+=(--partition="${SLURM_PARTITION}")
        [[ -n "${SLURM_ACCOUNT}"   ]] && SB+=(--account="${SLURM_ACCOUNT}")
        [[ -n "${SLURM_TIME}"      ]] && SB+=(--time="${SLURM_TIME}")
        # shellcheck disable=SC2206
        [[ -n "${SLURM_EXTRA}"     ]] && SB+=(${SLURM_EXTRA})
        ARRAY_SPEC="1-${#GENES[@]}"
        [[ -n "${ARRAY_THROTTLE}" ]] && ARRAY_SPEC="${ARRAY_SPEC}%${ARRAY_THROTTLE}"

        EXPORTS="ALL,MORPHEUS_HOME=${MORPHEUS_HOME},MORPHEUS_PROJECT_DIR=${WORKING_DIR}"
        EXPORTS+=",GENE_LIST=${GENE_LIST},SPECIES_TREE=${SPECIES_TREE}"
        EXPORTS+=",BAT_ANNOTATION_DIR=${BAT_ANNOTATION_DIR},HUMAN_GENOME_DIR=${HUMAN_GENOME_DIR}"
        EXPORTS+=",WORKING_DIR=${WORKING_DIR},CACHE_DIR=${CACHE_DIR}"
        EXPORTS+=",MORPHEUS_SEARCH=${SEARCH},MORPHEUS_SKIP_PLOT=${SKIP_PLOT}"
        EXPORTS+=",MORPHEUS_ENV_MODE=${ENV_MODE},MORPHEUS_ENV_NAME=${ENV_NAME}"
        [[ -n "${VENV_PATH}" ]] && EXPORTS+=",MORPHEUS_VENV_PATH=${VENV_PATH}"
        [[ ${SLURM_MODULE_SET} -eq 1 ]] && EXPORTS+=",MORPHEUS_MODULE=${SLURM_MODULE_ARG}"

        if [[ ${DRY_RUN} -eq 1 ]]; then
            cat <<EOF
would submit ${#GENES[@]} genes, one job each${SB[0]+, with ${SB[*]}}:
  00_reference.sbatch          flatten the GTF once (shared cache)
  01_gene_array.sbatch         --array=${ARRAY_SPEC}   one task per gene
  04_collect.sbatch            merge, figures, summary
genes: ${GENES[*]}
with --export=${EXPORTS}
EOF
            exit 0
        fi

        S="${MORPHEUS_HOME}/slurm"
        submit() { sbatch --parsable --export="${EXPORTS}" ${SB[@]+"${SB[@]}"} "$@"; }
        ref=$(submit "${S}/00_reference.sbatch")
        arr=$(submit --dependency=afterok:"${ref}" --array="${ARRAY_SPEC}" \
                     "${S}/01_gene_array.sbatch")
        # afterany, not afterok: one gene failing must not stop the other
        # twenty being merged and plotted. merge-genes names any gene that
        # produced nothing.
        coll=$(submit --dependency=afterany:"${arr}" "${S}/04_collect.sbatch")
        cat <<EOF

Submitted ${#GENES[@]} genes, one job each.

  ${ref}   reference   flatten the GTF once, into ${CACHE_DIR}
  ${arr}   genes       array ${ARRAY_SPEC}, one task per gene
  ${coll}   collect     merge, figures, summary (runs even if a gene fails)

  squeue -j ${ref},${arr},${coll}
  logs        ${WORKING_DIR}/logs/
  per gene    ${GENES_ROOT}/<GENE>/results/
  combined    ${COMBINED_DIR}/
EOF
        exit 0
    fi

    # ---- local: one gene at a time -----------------------------------------
    # Sequential on purpose. A long gene list run in parallel on a workstation
    # is how you end up with a machine you cannot type on.
    mkdir -p "${GENES_ROOT}" "${WORKING_DIR}/logs"
    echo "Running ${#GENES[@]} genes one at a time. Each gets its own results tree."
    failed=()
    for i in "${!GENES[@]}"; do
        g="${GENES[$i]}"
        printf '\n===== [%d/%d] %s =====\n' "$((i + 1))" "${#GENES[@]}" "${g}"
        if bash "${BASH_SOURCE[0]}" --_single_gene "${g}" --gene "${g}" \
                --output "${GENES_ROOT}/${g}" "${PASS[@]}"; then
            :
        else
            echo "WARNING: ${g} failed; continuing with the rest" >&2
            failed+=("${g}")
        fi
    done

    echo
    if [[ ${#failed[@]} -gt 0 ]]; then
        echo "${#failed[@]} of ${#GENES[@]} genes failed: ${failed[*]}" >&2
        echo "Re-run one with: morpheus run --gene <GENE> --output ${GENES_ROOT}/<GENE>" >&2
    fi

    MORPHEUS_SKIP_PLOT="${SKIP_PLOT}" \
        bash "${MORPHEUS_HOME}/scripts/merge_and_plot.sh" "${GENES_ROOT}" "${COMBINED_DIR}"

    echo
    echo "Per gene   ${GENES_ROOT}/<GENE>/results/"
    echo "Combined   ${COMBINED_DIR}/"
    [[ ${#failed[@]} -eq 0 ]] || exit 1
    exit 0
fi

# ---------------------------------------------------------------------------
# --mode slurm: submit the whole chain and return. Every resolved setting is
# handed to the jobs, so the batch run sees exactly what was resolved here
# rather than re-reading paths.txt from a different working directory.
# ---------------------------------------------------------------------------
if [[ "${MODE}" == "slurm" ]]; then
    [[ -d "${BAT_ANNOTATION_DIR}" ]] || die "missing bat annotation directory: ${BAT_ANNOTATION_DIR}"
    [[ -s "${GENE_LIST}" ]] || die "missing gene list: ${GENE_LIST} (use --gene or --gene_list)"
    [[ -s "${SPECIES_TREE}" ]] || die "missing species tree: ${SPECIES_TREE} (use --tree)"
    n_genomes=$(find "${BAT_ANNOTATION_DIR}" -maxdepth 1 -mindepth 1 -type d ! -name '.*' | wc -l | tr -d ' ')
    [[ "${n_genomes}" -gt 0 ]] || die "no genome sub-directories in ${BAT_ANNOTATION_DIR}"

    mkdir -p "${WORKING_DIR}/logs"

    # Site options, applied to all four jobs. Anything here overrides the
    # #SBATCH defaults in the scripts, which is what sbatch's own precedence
    # rules already do for command-line flags.
    SB=()
    [[ -n "${SLURM_PARTITION}" ]] && SB+=(--partition="${SLURM_PARTITION}")
    [[ -n "${SLURM_ACCOUNT}"   ]] && SB+=(--account="${SLURM_ACCOUNT}")
    [[ -n "${SLURM_TIME}"      ]] && SB+=(--time="${SLURM_TIME}")
    # shellcheck disable=SC2206
    [[ -n "${SLURM_EXTRA}"     ]] && SB+=(${SLURM_EXTRA})
    ARRAY_SPEC="1-${n_genomes}"
    [[ -n "${ARRAY_THROTTLE}" ]] && ARRAY_SPEC="${ARRAY_SPEC}%${ARRAY_THROTTLE}"

    EXPORTS="ALL,MORPHEUS_HOME=${MORPHEUS_HOME},MORPHEUS_PROJECT_DIR=${WORKING_DIR}"
    [[ ${SLURM_MODULE_SET} -eq 1 ]] && EXPORTS+=",MORPHEUS_MODULE=${SLURM_MODULE_ARG}"
    EXPORTS+=",MORPHEUS_ENV_MODE=${ENV_MODE},MORPHEUS_ENV_NAME=${ENV_NAME}"
    [[ -n "${VENV_PATH}" ]] && EXPORTS+=",MORPHEUS_VENV_PATH=${VENV_PATH}"
    EXPORTS+=",GENE_LIST=${GENE_LIST},SPECIES_TREE=${SPECIES_TREE}"
    EXPORTS+=",BAT_ANNOTATION_DIR=${BAT_ANNOTATION_DIR},HUMAN_GENOME_DIR=${HUMAN_GENOME_DIR}"
    EXPORTS+=",WORKING_DIR=${WORKING_DIR},MORPHEUS_SEARCH=${SEARCH}"
    EXPORTS+=",MORPHEUS_SKIP_PLOT=${SKIP_PLOT}"

    if [[ ${DRY_RUN} -eq 1 ]]; then
        cat <<EOF
would submit${SB[0]+, with ${SB[*]}}:
  00_reference.sbatch
  02_genome_array.sbatch       --array=${ARRAY_SPEC}
  03_after_search.sbatch
  04_collect.sbatch
with --export=${EXPORTS}
EOF
        exit 0
    fi

    command -v sbatch >/dev/null 2>&1 || die "--mode slurm, but 'sbatch' is not on PATH"
    S="${MORPHEUS_HOME}/slurm"
    submit() { sbatch --parsable --export="${EXPORTS}" ${SB[@]+"${SB[@]}"} "$@"; }
    ref=$(submit "${S}/00_reference.sbatch")
    arr=$(submit --dependency=afterok:"${ref}" --array="${ARRAY_SPEC}" \
                 "${S}/02_genome_array.sbatch")
    post=$(submit --dependency=afterok:"${arr}" "${S}/03_after_search.sbatch")
    coll=$(submit --dependency=afterok:"${post}" "${S}/04_collect.sbatch")
    cat <<EOF

Submitted the whole chain; nothing else to do.

  ${ref}   reference     steps 1-3
  ${arr}   search        array ${ARRAY_SPEC}, one task per genome
  ${post}   post-search   screen, assign, copy number, sequences, align
  ${coll}   collect       $( [[ ${SKIP_PLOT} -eq 1 ]] && echo "summary only (--skip_plot)" || echo "figures and summary" )

  squeue -j ${ref},${arr},${post},${coll}
  logs      ${WORKING_DIR}/logs/
  results   ${RESULTS_DIR}/
EOF
    exit 0
fi

# ---------------------------------------------------------------------------
# checks
# ---------------------------------------------------------------------------
activate_env || exit 1
# Create every step directory up front. Otherwise a step that has not run yet
# simply leaves no directory, and the numbering reads as though a stage were
# missing or the run had broken.
STEP_DIRS=("${REFERENCE_DIR}" "${SEARCH_DIR}" "${ASSIGN_DIR}" "${COPY_DIR}" "${GENES_DIR}")
[[ ${SKIP_PLOT} -eq 1 ]] || STEP_DIRS+=("${PLOTS_DIR}")
mkdir -p "${RESULTS_DIR}" "${LOG_DIR}" "${CACHE_DIR}" "${STEP_DIRS[@]}"

# Mark the step directories that are still empty, so the numbering never reads
# as though a stage were missing or the run had broken. Called again at the end,
# because a directory filled during this run must lose its marker.
mark_empty_steps() {
    local d
    for d in "${STEP_DIRS[@]}"; do
        [[ -d "$d" ]] || continue
        # "Empty" must ignore the marker itself, or the marker counts as content
        # on the next invocation and deletes itself.
        if [[ -z "$(find "$d" -mindepth 1 ! -name '.not_run_yet' -print -quit 2>/dev/null)" ]]; then
            printf '%s\n' \
              "This step has not produced results yet." \
              "Run: morpheus run --only <step>   (see morpheus run --help)" \
              > "$d/.not_run_yet"
        else
            rm -f "$d/.not_run_yet"
        fi
    done
}
mark_empty_steps

missing=0
for f in "${GENE_LIST}" "${SPECIES_TREE}" "${HUMAN_GTF}" "${HUMAN_CDS_FASTA}" "${HUMAN_LONGEST_PC_TX}"; do
    [[ -s "$f" ]] || { echo "ERROR: missing input: ${f:-<unset>}" >&2; missing=1; }
done
[[ -d "${BAT_ANNOTATION_DIR}" ]] || { echo "ERROR: missing ${BAT_ANNOTATION_DIR}" >&2; missing=1; }
NEEDED=(python3 blastx blastp makeblastdb mafft)
[[ ${SKIP_PLOT} -eq 1 ]] || NEEDED+=(Rscript)
for exe in "${NEEDED[@]}"; do
    command -v "$exe" >/dev/null || { echo "ERROR: '$exe' not on PATH" >&2; missing=1; }
done
if [[ $missing -ne 0 ]]; then
    echo "
Inputs come from --path_file/paths.txt or from --bat_dir/--ref_dir/--output.
Run 'morpheus env' to see what was resolved." >&2
    exit 1
fi

n_species=$(find "${BAT_ANNOTATION_DIR}" -maxdepth 1 -mindepth 1 -type d ! -name '.*' | wc -l | tr -d ' ')
[[ ${#SPECIES[@]-0} -gt 0 ]] && n_species=${#SPECIES[@]}
PLOT_MIN_SPECIES="$(plot_min_species_for "${n_species}")"
cat <<EOF
Morpheus $(cat "${MORPHEUS_HOME}/VERSION" 2>/dev/null || echo "")
  genes          $(grep -cve '^\s*$' "${GENE_LIST}") from ${GENE_LIST}
  bat genomes    ${n_species} in ${BAT_ANNOTATION_DIR}
  human GTF      $(basename "${HUMAN_GTF}")
  search         ${SEARCH}  (scopes: ${SCOPES} | rankings: ${POLICIES})
  min species    ${MIN_SPECIES_FRACTION} of ${n_species} = ${PLOT_MIN_SPECIES}
  figures        $( [[ ${SKIP_PLOT} -eq 1 ]] && echo "skipped (--skip_plot)" || echo "on" )
  results        ${RESULTS_DIR}
  threads        ${THREADS}
EOF
[[ $DRY_RUN -eq 1 ]] && { echo; for s in "${ALL_STEPS[@]}"; do should_run "$s" && echo "would run: $s"; done; exit 0; }

run() { "${MORPHEUS[@]}" "$@"; }

# ---------------------------------------------------------------------------
# 1. flatten the Ensembl GTF into a CDS exon table (cached)
# ---------------------------------------------------------------------------
if should_run cds-table; then
    banner "cds-table"
    run cds-table --gtf "${HUMAN_GTF}" --out "${CDS_EXON_TABLE}" \
        2>&1 | tee "${LOG_DIR}/01_cds_table.log"
fi

# ---------------------------------------------------------------------------
# 2. human reference for the genes of interest
# ---------------------------------------------------------------------------
if should_run human-reference; then
    banner "human-reference"
    run human-reference --genes "${GENE_LIST}" --cds-exons "${CDS_EXON_TABLE}" \
        --longest-pc-tx "${HUMAN_LONGEST_PC_TX}" --cds "${HUMAN_CDS_FASTA}" \
        --outdir "${REFERENCE_DIR}" 2>&1 | tee "${LOG_DIR}/02_human_reference.log"
fi

# ---------------------------------------------------------------------------
# 3. paralog family of each gene, from the human proteome
# ---------------------------------------------------------------------------
if should_run families; then
    banner "families"
    run families --reference "${REFERENCE_DIR}" --outdir "${REFERENCE_DIR}" \
        --min-fraction "${FAMILY_MIN_FRACTION}" --threads "${THREADS}" \
        2>&1 | tee "${LOG_DIR}/03_families.log"
fi

# ---------------------------------------------------------------------------
# 4. region-restricted transcript search in every bat genome
# ---------------------------------------------------------------------------
if should_run bat-search; then
    banner "bat-search"
    species_args=()
    [[ ${#SPECIES[@]-0} -gt 0 ]] && species_args=(--species "${SPECIES[@]}")
    run bat-search --annotations "${BAT_ANNOTATION_DIR}" \
        --reference "${REFERENCE_DIR}" \
        --families "${REFERENCE_DIR}/gene_families.tsv" \
        --tree "${SPECIES_TREE}" \
        --outdir "${SEARCH_DIR}" --locus-slop "${LOCUS_SLOP}" \
        ${OVERWRITE} ${species_args[@]+"${species_args[@]}"} \
        2>&1 | tee "${LOG_DIR}/04_bat_search.log"
fi

CANDIDATES="${SEARCH_DIR}/all_species_candidates.tsv"
CANDIDATE_FASTA="${SEARCH_DIR}/all_species_candidate_cds.fa"
LOCI="${SEARCH_DIR}/all_species_loci.tsv"

# ---------------------------------------------------------------------------
# 5. paralog screen against the whole human proteome
# ---------------------------------------------------------------------------
if should_run screen; then
    banner "screen"
    run screen --candidates "${CANDIDATES}" --candidate-fasta "${CANDIDATE_FASTA}" \
        --reference "${REFERENCE_DIR}" --outdir "${SEARCH_DIR}" --threads "${THREADS}" \
        2>&1 | tee "${LOG_DIR}/05_screen.log"
fi

SCREEN="${SEARCH_DIR}/paralog_screen.tsv"

# ---------------------------------------------------------------------------
# 6. pairwise BLAST and one-to-one transcript assignment
#      region_restricted - the gene's own syntenic locus only
#      unrestricted      - the gene anywhere, paralogs still excluded
# ---------------------------------------------------------------------------
if should_run assign; then
    banner "assign"
    # shellcheck disable=SC2086
    run assign --candidates "${CANDIDATES}" --candidate-fasta "${CANDIDATE_FASTA}" \
        --reference "${REFERENCE_DIR}" --screen "${SCREEN}" --outdir "${ASSIGN_DIR}" \
        --min-score "${MIN_ASSIGNMENT_SCORE}" --min-pident "${MIN_ASSIGNMENT_PIDENT}" \
        --scopes ${SCOPES} --policies ${POLICIES} \
        --threads "${THREADS}" 2>&1 | tee "${LOG_DIR}/06_assign.log"
fi

# ---------------------------------------------------------------------------
# 6b. where do the two scopes disagree?  (--search both only)
# ---------------------------------------------------------------------------
if should_run compare; then
    banner "compare"
    : > "${LOG_DIR}/06b_compare.log"
    for policy in ${POLICIES}; do
        run compare \
            --region-restricted "${ASSIGN_DIR}/transcript_assignments_region_restricted__${policy}.tsv" \
            --unrestricted "${ASSIGN_DIR}/transcript_assignments_unrestricted__${policy}.tsv" \
            --candidate-fasta "${CANDIDATE_FASTA}" --outdir "${ASSIGN_DIR}" \
            --policy "${policy}" 2>&1 | tee -a "${LOG_DIR}/06b_compare.log"
    done
fi

# ---------------------------------------------------------------------------
# 7. region-aware gene copy number
# ---------------------------------------------------------------------------
if should_run copy-number; then
    banner "copy-number"
    run copy-number --candidates "${CANDIDATES}" --candidate-fasta "${CANDIDATE_FASTA}" \
        --screen "${SCREEN}" --loci "${LOCI}" --outdir "${COPY_DIR}" \
        --families "${REFERENCE_DIR}/gene_families.tsv" \
        --cluster-slop "${COPY_CLUSTER_SLOP}" 2>&1 | tee "${LOG_DIR}/07_copy_number.log"
fi

# ---------------------------------------------------------------------------
# 8. sequence directories, and a transcript-status table for each scope
# ---------------------------------------------------------------------------
if should_run deliverables; then
    banner "deliverables"
    : > "${LOG_DIR}/08_deliverables.log"
    # One sequence tree per ranking policy, so the views can be compared gene by
    # gene instead of one being chosen on the reader's behalf.
    for policy in ${POLICIES}; do
        run deliverables \
            --assignments "${ASSIGN_DIR}/transcript_assignments_${SEQUENCE_SCOPE}__${policy}.tsv" \
            --candidate-fasta "${CANDIDATE_FASTA}" --reference "${REFERENCE_DIR}" \
            --tree "${SPECIES_TREE}" --outdir "${GENES_DIR}" --policy "${policy}" \
            --min-species-fraction "${MIN_SPECIES_FRACTION}" \
            --total-species "${n_species}" --scope "${SEQUENCE_SCOPE}" \
            2>&1 | tee -a "${LOG_DIR}/08_deliverables.log"
        for scope in ${STATUS_SCOPES}; do
            # the sequence scope already wrote its own status table above
            [[ "${scope}" == "${SEQUENCE_SCOPE}" ]] && continue
            run deliverables \
                --assignments "${ASSIGN_DIR}/transcript_assignments_${scope}__${policy}.tsv" \
                --candidate-fasta "${CANDIDATE_FASTA}" --reference "${REFERENCE_DIR}" \
                --tree "${SPECIES_TREE}" --outdir "${GENES_DIR}" --policy "${policy}" \
                --min-species-fraction "${MIN_SPECIES_FRACTION}" \
                --total-species "${n_species}" --scope "${scope}" \
                --status-only 2>&1 | tee -a "${LOG_DIR}/08_deliverables.log"
        done
    done
fi

# ---------------------------------------------------------------------------
# 9. codon-aware alignments
# ---------------------------------------------------------------------------
if should_run align; then
    banner "align"
    : > "${LOG_DIR}/09_align.log"
    for policy in ${POLICIES}; do
        run align --manifest "${GENES_DIR}/${policy}/manifest.tsv" \
            --outdir "${GENES_DIR}/${policy}" --threads "${THREADS}" \
            2>&1 | tee -a "${LOG_DIR}/09_align.log"
    done
fi

# ---------------------------------------------------------------------------
# 10. plots
# ---------------------------------------------------------------------------
PLOTS_FAILED=0
# Figures are secondary to the tables. A threshold set too high for this
# particular subset, or a missing R package, must not turn a run that produced
# every table into a failed run -- especially under --per_gene, where the driver
# would then report the gene itself as failed.
draw() {
    if ! Rscript "$@"; then
        echo "WARNING: figure failed: Rscript $*" >&2
        PLOTS_FAILED=$((PLOTS_FAILED + 1))
    fi
}
if should_run plots; then
    banner "plots"
    mkdir -p "${PLOTS_DIR}"
    : > "${LOG_DIR}/10_plot_status.log"; : > "${LOG_DIR}/10b_plot_scope.log"
    # Transcript status depends on the ranking, so each policy gets its own
    # figures. Copy number does not depend on it and stays at the top level.
    for policy in ${POLICIES}; do
        for scope in ${STATUS_SCOPES}; do
            draw "${R_DIR}/plot_transcript_status.R" \
                --status "${GENES_DIR}/${policy}/transcript_status_${scope}.tsv" \
                --tree "${SPECIES_TREE}" --outdir "${PLOTS_DIR}/${policy}" \
                --min-species "${PLOT_MIN_SPECIES}" \
                --stem "transcript_status_${scope}" \
                2>&1 | tee -a "${LOG_DIR}/10_plot_status.log"
        done
        if [[ ${DO_COMPARE} -eq 1 ]]; then
            draw "${R_DIR}/plot_scope_comparison.R" \
                --comparison "${ASSIGN_DIR}/scope_comparison__${policy}.tsv" \
                --tree "${SPECIES_TREE}" --outdir "${PLOTS_DIR}/${policy}" \
                --min-species "${PLOT_MIN_SPECIES}" \
                2>&1 | tee -a "${LOG_DIR}/10b_plot_scope.log"
        fi
    done
    draw "${R_DIR}/plot_copy_number.R" \
        --matrix "${COPY_DIR}/copy_number_matrix.tsv" \
        --tree "${SPECIES_TREE}" --outdir "${PLOTS_DIR}" \
        2>&1 | tee "${LOG_DIR}/11_plot_copy_number.log"
    [[ ${PLOTS_FAILED} -eq 0 ]] || echo "WARNING: ${PLOTS_FAILED} figure(s) failed; every table was still written" >&2
fi

# ---------------------------------------------------------------------------
# 11. one-page summary joining every stage
# ---------------------------------------------------------------------------
if should_run summary; then
    banner "summary"
    : > "${LOG_DIR}/12_summary.log"
    for policy in ${POLICIES}; do
        run summary --results "${RESULTS_DIR}" --policy "${policy}" \
            2>&1 | tee -a "${LOG_DIR}/12_summary.log"
    done
fi

mark_empty_steps

banner "done"
cat <<EOF
Results:
  human reference       ${REFERENCE_DIR}
  bat search            ${SEARCH_DIR}
  transcript assignment ${ASSIGN_DIR}/transcript_assignments_<scope>__<policy>.tsv
  copy number           ${COPY_DIR}/copy_number_matrix.tsv
EOF
for policy in ${POLICIES}; do
    echo "  sequences             ${GENES_DIR}/${policy}/{all_gene_files,gene_files_selected}/"
done
[[ ${DO_COMPARE} -eq 1 ]] && echo "  scope comparison      ${ASSIGN_DIR}/scope_comparison__<policy>.tsv"
if [[ ${SKIP_PLOT} -eq 1 ]]; then
    echo "  plots                 skipped (--skip_plot)"
else
    echo "  plots                 ${PLOTS_DIR}/<policy>/  (copy number at the top level)"
fi
echo
echo "Start here:"
echo "  ${RESULTS_DIR}/SUMMARY__${DEFAULT_POLICY}.md"
