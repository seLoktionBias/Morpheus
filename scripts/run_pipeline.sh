#!/usr/bin/env bash
###############################################################################
# run_pipeline.sh - run the whole analysis, or any part of it.
#
#   morpheus run                            # everything
#   morpheus run --from assign              # resume from a step
#   morpheus run --only copy-number         # one step
#   morpheus run --species Myotis_myotis Desmodus_rotundus
#
# Steps, in order:
#   cds-table  human-reference  families  bat-search  screen  assign
#   compare  copy-number  deliverables  align  hyphy  plots  summary
#
# Every step writes into results/ and reads only what earlier steps produced,
# so re-running one step never invalidates the others.
###############################################################################
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MORPHEUS_HOME="${MORPHEUS_HOME:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
source "${MORPHEUS_HOME}/config/config.sh"

ALL_STEPS=(cds-table human-reference families bat-search screen assign
           compare copy-number deliverables align hyphy plots summary)

FROM=""; ONLY=""; SPECIES=(); OVERWRITE=""; DRY_RUN=0

usage() {
    sed -n '3,17p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    cat <<'EOF'

Options:
  --from STEP        start at STEP and run everything after it
  --only STEP        run just STEP (repeatable)
  --species NAME...  restrict the bat search to these Genus_species
  --overwrite        redo per-species searches that are already complete
  --dry-run          print what would run, then stop
  -h, --help         this message

Key settings (override by exporting before running):
  THREADS, MIN_SPECIES_FRACTION, MIN_ASSIGNMENT_SCORE, HYPHY_METHODS, HYPHY_TIMEOUT
EOF
}

declare -a ONLY_STEPS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --from) FROM="$2"; shift 2 ;;
        --only) ONLY_STEPS+=("$2"); shift 2 ;;
        --species) shift; while [[ $# -gt 0 && "$1" != --* ]]; do SPECIES+=("$1"); shift; done ;;
        --overwrite) OVERWRITE="--overwrite"; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option '$1'" >&2; usage >&2; exit 1 ;;
    esac
done

should_run() {
    local step="$1"
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
# checks
# ---------------------------------------------------------------------------
activate_env || exit 1
# Create every step directory up front. Otherwise a step that has not run yet
# simply leaves no directory, and the numbering reads as though a stage were
# missing or the run had broken.
mkdir -p "${RESULTS_DIR}" "${LOG_DIR}" "${CACHE_DIR}" \
         "${REFERENCE_DIR}" "${SEARCH_DIR}" "${ASSIGN_DIR}" "${COPY_DIR}" \
         "${GENES_DIR}" "${SELECTION_DIR}" "${PLOTS_DIR}"
for d in "${REFERENCE_DIR}" "${SEARCH_DIR}" "${ASSIGN_DIR}" "${COPY_DIR}" \
         "${GENES_DIR}" "${SELECTION_DIR}" "${PLOTS_DIR}"; do
    # "Empty" must ignore the marker itself, or the marker counts as content on
    # the next invocation and deletes itself.
    if [[ -z "$(find "$d" -mindepth 1 ! -name '.not_run_yet' -print -quit 2>/dev/null)" ]]; then
        printf '%s\n' \
          "This step has not produced results yet." \
          "Run: morpheus run --only <step>   (see morpheus run --help)" \
          > "$d/.not_run_yet"
    else
        rm -f "$d/.not_run_yet"
    fi
done

missing=0
for f in "${GENE_LIST}" "${SPECIES_TREE}" "${HUMAN_GTF}" "${HUMAN_CDS_FASTA}" "${HUMAN_LONGEST_PC_TX}"; do
    [[ -s "$f" ]] || { echo "ERROR: missing input: ${f:-<unset>}" >&2; missing=1; }
done
[[ -d "${BAT_ANNOTATION_DIR}" ]] || { echo "ERROR: missing ${BAT_ANNOTATION_DIR}" >&2; missing=1; }
for exe in python3 blastx blastp makeblastdb mafft hyphy Rscript; do
    command -v "$exe" >/dev/null || { echo "ERROR: '$exe' not on PATH" >&2; missing=1; }
done
[[ $missing -eq 0 ]] || exit 1

n_species=$(find "${BAT_ANNOTATION_DIR}" -maxdepth 1 -mindepth 1 -type d ! -name '.*' | wc -l | tr -d ' ')
cat <<EOF
Morpheus $(cat "${MORPHEUS_HOME}/VERSION" 2>/dev/null || echo "")
  genes          $(grep -cve '^\s*$' "${GENE_LIST}") from ${GENE_LIST}
  bat genomes    ${n_species} in ${BAT_ANNOTATION_DIR}
  human GTF      $(basename "${HUMAN_GTF}")
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
# 6. pairwise BLAST and one-to-one transcript assignment, under both scopes
#      region_restricted - the gene's own syntenic locus only
#      unrestricted      - the gene anywhere, paralogs still excluded
# ---------------------------------------------------------------------------
if should_run assign; then
    banner "assign"
    run assign --candidates "${CANDIDATES}" --candidate-fasta "${CANDIDATE_FASTA}" \
        --reference "${REFERENCE_DIR}" --screen "${SCREEN}" --outdir "${ASSIGN_DIR}" \
        --min-score "${MIN_ASSIGNMENT_SCORE}" --min-pident "${MIN_ASSIGNMENT_PIDENT}" \
        --scopes region_restricted unrestricted \
        --policies ${POLICIES} \
        --threads "${THREADS}" 2>&1 | tee "${LOG_DIR}/06_assign.log"
fi

# The scope comparison and the status plots read the synteny-aware ranking:
# orthology is positional, so that is the conservative default for "what does
# this locus produce". The sequence-similarity ranking is kept alongside it.
REGION_ASSIGN="${ASSIGN_DIR}/transcript_assignments_region_restricted__${DEFAULT_POLICY}.tsv"
FREE_ASSIGN="${ASSIGN_DIR}/transcript_assignments_unrestricted__${DEFAULT_POLICY}.tsv"

# ---------------------------------------------------------------------------
# 6b. where do the two scopes disagree?
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
# 8. sequence directories (unrestricted, for alignment/selection) and a
#    transcript-status table for each scope
# ---------------------------------------------------------------------------
if should_run deliverables; then
    banner "deliverables"
    : > "${LOG_DIR}/08_deliverables.log"
    # One sequence tree per ranking policy, so the two views can be compared
    # gene by gene instead of one being chosen on the reader's behalf.
    for policy in ${POLICIES}; do
        # sequences from the unrestricted scope ...
        run deliverables \
            --assignments "${ASSIGN_DIR}/transcript_assignments_unrestricted__${policy}.tsv" \
            --candidate-fasta "${CANDIDATE_FASTA}" --reference "${REFERENCE_DIR}" \
            --tree "${SPECIES_TREE}" --outdir "${GENES_DIR}" --policy "${policy}" \
            --min-species-fraction "${MIN_SPECIES_FRACTION}" \
            --total-species "${n_species}" --scope unrestricted \
            2>&1 | tee -a "${LOG_DIR}/08_deliverables.log"
        # ... and a status table for each scope, inside the same policy tree
        run deliverables \
            --assignments "${ASSIGN_DIR}/transcript_assignments_region_restricted__${policy}.tsv" \
            --candidate-fasta "${CANDIDATE_FASTA}" --reference "${REFERENCE_DIR}" \
            --tree "${SPECIES_TREE}" --outdir "${GENES_DIR}" --policy "${policy}" \
            --min-species-fraction "${MIN_SPECIES_FRACTION}" \
            --total-species "${n_species}" --scope region_restricted \
            --status-only 2>&1 | tee -a "${LOG_DIR}/08_deliverables.log"
    done
fi

MANIFEST="${GENES_DIR}/${DEFAULT_POLICY}/manifest.tsv"

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
# 10. HyPhy selection analyses (no branch labelling required)
# ---------------------------------------------------------------------------
if should_run hyphy; then
    banner "hyphy"
    hyphy_extra=()
    [[ -n "${HYPHY_TIMEOUT}" ]] && hyphy_extra+=(--timeout "${HYPHY_TIMEOUT}")
    [[ "${HYPHY_LONGEST_ONLY}" == "1" ]] && hyphy_extra+=(--longest-only)
    : > "${LOG_DIR}/10_hyphy.log"
    for policy in ${POLICIES}; do
        # shellcheck disable=SC2086
        run hyphy --manifest "${GENES_DIR}/${policy}/manifest.tsv" \
            --outdir "${SELECTION_DIR}/${policy}" \
            --methods ${HYPHY_METHODS} --jobs "${HYPHY_JOBS}" \
            --hyphy-threads "${HYPHY_THREADS}" \
            ${hyphy_extra[@]+"${hyphy_extra[@]}"} \
            2>&1 | tee -a "${LOG_DIR}/10_hyphy.log"
    done
fi

# ---------------------------------------------------------------------------
# 11. plots
# ---------------------------------------------------------------------------
if should_run plots; then
    banner "plots"
    mkdir -p "${PLOTS_DIR}"
    : > "${LOG_DIR}/11_plot_status.log"; : > "${LOG_DIR}/11b_plot_scope.log"
    # Transcript status and the scope comparison depend on the ranking, so each
    # policy gets its own figures. Copy number does not depend on it at all and
    # stays at the top level.
    for policy in ${POLICIES}; do
        for scope in region_restricted unrestricted; do
            Rscript "${R_DIR}/plot_transcript_status.R" \
                --status "${GENES_DIR}/${policy}/transcript_status_${scope}.tsv" \
                --tree "${SPECIES_TREE}" --outdir "${PLOTS_DIR}/${policy}" \
                --min-species "${PLOT_MIN_SPECIES}" \
                --stem "transcript_status_${scope}" \
                2>&1 | tee -a "${LOG_DIR}/11_plot_status.log"
        done
        Rscript "${R_DIR}/plot_scope_comparison.R" \
            --comparison "${ASSIGN_DIR}/scope_comparison__${policy}.tsv" \
            --tree "${SPECIES_TREE}" --outdir "${PLOTS_DIR}/${policy}" \
            --min-species "${PLOT_MIN_SPECIES}" \
            2>&1 | tee -a "${LOG_DIR}/11b_plot_scope.log"
    done
    Rscript "${R_DIR}/plot_copy_number.R" \
        --matrix "${COPY_DIR}/copy_number_matrix.tsv" \
        --tree "${SPECIES_TREE}" --outdir "${PLOTS_DIR}" \
        2>&1 | tee "${LOG_DIR}/12_plot_copy_number.log"
fi

# ---------------------------------------------------------------------------
# 12. one-page summary joining every stage
# ---------------------------------------------------------------------------
if should_run summary; then
    banner "summary"
    : > "${LOG_DIR}/13_summary.log"
    for policy in ${POLICIES}; do
        run summary --results "${RESULTS_DIR}" --policy "${policy}" \
            2>&1 | tee -a "${LOG_DIR}/13_summary.log"
    done
fi

banner "done"
cat <<EOF
Results:
  human reference       ${REFERENCE_DIR}
  bat search            ${SEARCH_DIR}
  transcript assignment ${ASSIGN_DIR}/transcript_assignments_<scope>__<policy>.tsv
  scope comparison      ${ASSIGN_DIR}/scope_comparison.tsv
  copy number           ${COPY_DIR}/copy_number_matrix.tsv
  sequence-similarity   ${GENES_DIR}/sequence_similarity/{all_gene_files,gene_files_selected}/
  synteny-aware         ${GENES_DIR}/synteny_aware/{all_gene_files,gene_files_selected}/
  selection             ${SELECTION_DIR}/<policy>/hyphy_gene_level.tsv
  plots                 ${PLOTS_DIR}/<policy>/  (copy number at the top level)

Start here:
  ${RESULTS_DIR}/SUMMARY__${DEFAULT_POLICY}.md
EOF
