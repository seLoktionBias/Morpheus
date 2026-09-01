#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Morpheus configuration
#
# Everything the pipeline needs to know about this machine lives here. Input
# locations are read from paths.txt so they stay in one place; override any
# value by exporting it before sourcing this file.
# ---------------------------------------------------------------------------

# MORPHEUS_HOME is where the package lives; PROJECT_DIR is where the data and
# results live. They are usually different, so a checkout can be shared and each
# analysis keeps its own working directory.
MORPHEUS_HOME="${MORPHEUS_HOME:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
PROJECT_DIR="${MORPHEUS_PROJECT_DIR:-$(pwd)}"

# ---- input locations (from paths.txt) -------------------------------------
# Anything already exported wins, so a command-line flag always beats the file.
# Keys are matched case-insensitively with synonyms, because "what exactly do I
# call this variable" should not be something anyone has to guess or grep for.
PATHS_FILE="${PATHS_FILE:-${PROJECT_DIR}/paths.txt}"
if [[ -s "${PATHS_FILE}" ]]; then
    _lineno=0
    while IFS= read -r _raw || [[ -n "${_raw}" ]]; do
        _lineno=$((_lineno + 1))
        _line="${_raw%$'\r'}"
        _line="${_line%%#*}"                       # strip trailing comments
        _line="${_line#"${_line%%[![:space:]]*}"}" # ltrim
        _line="${_line%"${_line##*[![:space:]]}"}" # rtrim
        [[ -z "${_line}" ]] && continue
        if [[ "${_line}" != *=* ]]; then
            echo "ERROR: ${PATHS_FILE} line ${_lineno}: expected key=value, got: ${_raw}" >&2
            return 1 2>/dev/null || exit 1
        fi
        _key="${_line%%=*}"; _val="${_line#*=}"
        _key="${_key//[[:space:]]/}"
        _key="$(printf '%s' "${_key}" | tr '[:upper:]' '[:lower:]')"
        _val="${_val#"${_val%%[![:space:]]*}"}"
        _val="${_val%"${_val##*[![:space:]]}"}"
        _val="${_val%\"}"; _val="${_val#\"}"       # tolerate quoting
        _val="${_val%\'}"; _val="${_val#\'}"
        _val="${_val/#\~/${HOME}}"                  # expand a leading ~
        case "${_key}" in
            human_genome_dir|human_dir|ref_dir|reference_dir)
                : "${HUMAN_GENOME_DIR:=${_val}}" ;;
            bat_annotation_dir|bat_dir|toga_dir|toga2_dir|annotation_dir)
                : "${BAT_ANNOTATION_DIR:=${_val}}" ;;
            primary_working_dir|working_dir|output|output_dir|outdir)
                : "${WORKING_DIR:=${_val}}" ;;
            tree|species_tree|tree_file)
                : "${SPECIES_TREE:=${_val}}" ;;
            gene_list|genes|gene_list_file)
                : "${GENE_LIST:=${_val}}" ;;
            *)
                echo "ERROR: ${PATHS_FILE} line ${_lineno}: unrecognised key '${_key}'." >&2
                echo "       Recognised keys: human_genome_dir, bat_annotation_dir," >&2
                echo "       primary_working_dir, tree (optional), gene_list (optional)." >&2
                return 1 2>/dev/null || exit 1 ;;
        esac
    done < "${PATHS_FILE}"
    unset _raw _line _key _val _lineno
fi

# Fall back to conventional sub-directories of the project, so a tidy layout
# needs no paths.txt at all.
: "${HUMAN_GENOME_DIR:=${PROJECT_DIR}/human_genome}"
: "${BAT_ANNOTATION_DIR:=${PROJECT_DIR}/bat1k_genomes}"
: "${WORKING_DIR:=${PROJECT_DIR}}"

# ---- project inputs -------------------------------------------------------
: "${GENE_LIST:=${PROJECT_DIR}/gene_list.txt}"
: "${SPECIES_TREE:=${PROJECT_DIR}/bat1k_tree.nwk}"

# Human annotation files. Names differ between Ensembl releases, so each is
# resolved by pattern with an explicit override available.
_pick() { for f in "$@"; do [[ -s "$f" ]] && { printf '%s' "$f"; return 0; }; done; return 1; }
: "${HUMAN_GTF:=$(_pick "${HUMAN_GENOME_DIR}"/Homo_sapiens.*.gtf.gz "${HUMAN_GENOME_DIR}"/*.gtf.gz)}"
: "${HUMAN_CDS_FASTA:=$(_pick "${HUMAN_GENOME_DIR}"/Homo_sapiens.*.cds.all.fa "${HUMAN_GENOME_DIR}"/*cds*.fa)}"
: "${HUMAN_LONGEST_PC_TX:=${HUMAN_GENOME_DIR}/longest_pc_tx.tsv}"

# ---- outputs --------------------------------------------------------------
: "${RESULTS_DIR:=${WORKING_DIR}/results}"
: "${LOG_DIR:=${RESULTS_DIR}/logs}"
: "${CACHE_DIR:=${RESULTS_DIR}/cache}"

REFERENCE_DIR="${RESULTS_DIR}/01_human_reference"
SEARCH_DIR="${RESULTS_DIR}/02_bat_search"
ASSIGN_DIR="${RESULTS_DIR}/03_transcript_assignment"
COPY_DIR="${RESULTS_DIR}/04_copy_number"
GENES_DIR="${RESULTS_DIR}/05_genes"
PLOTS_DIR="${RESULTS_DIR}/06_plots"

CDS_EXON_TABLE="${CACHE_DIR}/human_cds_exons.tsv.gz"

# ---- analysis parameters --------------------------------------------------
: "${THREADS:=$( (sysctl -n hw.logicalcpu 2>/dev/null || nproc 2>/dev/null || echo 4) )}"
: "${MIN_SPECIES_FRACTION:=0.5}"   # a transcript needs this share of query species
                                   # to enter gene_files_selected and be analysed
# Ranking policies. sequence_similarity asks which query sequence most resembles
# the transcript; synteny_aware asks which is the ortholog, judging by syntenic
# position and exon structure, so an intronless retrocopy cannot win on identity
# alone. Both are built; DEFAULT_POLICY drives the status plots and the summary.
: "${POLICIES:=sequence_similarity synteny_aware}"
: "${DEFAULT_POLICY:=synteny_aware}"
: "${PLOT_MIN_SPECIES:=50}"        # a transcript column needs this many species to be plotted
: "${MIN_ASSIGNMENT_SCORE:=0.30}"  # below this a human transcript stays unassigned
: "${MIN_ASSIGNMENT_PIDENT:=40}"   # protein identity floor for any assignment
: "${FAMILY_MIN_FRACTION:=0.25}"   # paralog family membership threshold
: "${COPY_CLUSTER_SLOP:=10000}"    # bp gap tolerated when merging projections into one copy
: "${LOCUS_SLOP:=0}"               # bp gap tolerated when clustering seed projections

# ---- environment ----------------------------------------------------------
: "${CONDA_ENV:=${MORPHEUS_ENV:-Morpheus}}"
: "${CONDA_ROOT:=}"

activate_env() {
    # Already inside the right environment, or a site module has already put the
    # tools on PATH (Apocrita: `ml miniforge`)? Nothing to do.
    if [[ "${CONDA_DEFAULT_ENV:-}" == "${CONDA_ENV}" ]]; then return 0; fi
    if command -v mafft >/dev/null 2>&1 && command -v blastx >/dev/null 2>&1 \
       && command -v Rscript >/dev/null 2>&1; then
        return 0
    fi

    local root="${CONDA_ROOT}"
    if [[ -z "${root}" ]]; then
        for candidate in "${CONDA_PREFIX_1:-}" "${HOME}/miniforge3" "${HOME}/mambaforge" \
                         "${HOME}/miniconda3" "${HOME}/anaconda3" \
                         "${HOME}/Downloads/miniconda3" "/opt/miniconda3" \
                         "/opt/homebrew/Caskroom/miniforge/base"; do
            [[ -z "${candidate}" ]] && continue
            [[ -x "${candidate}/bin/conda" ]] && { root="${candidate}"; break; }
        done
    fi
    if [[ -z "${root}" ]]; then
        echo "ERROR: no conda installation found. Set CONDA_ROOT." >&2
        return 1
    fi
    # Put the environment on PATH directly; this works the same in an
    # interactive shell and in a non-interactive script.
    local envdir="${root}/envs/${CONDA_ENV}"
    [[ -d "${envdir}" ]] || { echo "ERROR: conda env '${CONDA_ENV}' not found in ${root}/envs" >&2; return 1; }
    export PATH="${envdir}/bin:${PATH}"
    export CONDA_PREFIX="${envdir}"
    export CONDA_DEFAULT_ENV="${CONDA_ENV}"
    return 0
}

export PYTHONPATH="${MORPHEUS_HOME}${PYTHONPATH:+:${PYTHONPATH}}"
MORPHEUS=(python3 -m morpheus)
R_DIR="${MORPHEUS_HOME}/R"
