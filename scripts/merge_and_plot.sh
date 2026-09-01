#!/usr/bin/env bash
###############################################################################
# Recombine per-gene runs and draw the whole-list figures.
#
#   merge_and_plot.sh <genes_root> <combined_dir>
#
# Called by both paths that produce per-gene results -- the local sequential
# driver and the Slurm collect job -- so the merge happens one way, not two.
###############################################################################
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MORPHEUS_HOME="${MORPHEUS_HOME:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
source "${MORPHEUS_HOME}/config/config.sh"

GENES_ROOT="${1:?usage: merge_and_plot.sh <genes_root> <combined_dir>}"
COMBINED="${2:?usage: merge_and_plot.sh <genes_root> <combined_dir>}"
SKIP_PLOT="${MORPHEUS_SKIP_PLOT:-0}"

# Which rankings and scopes exist depends on --search, so read them from what
# the gene runs actually wrote rather than assuming the default pair.
# while-read, not mapfile: mapfile is bash 4+ and macOS ships bash 3.2.
FOUND_POLICIES=(); FOUND_SCOPES=()
while IFS= read -r x; do [[ -n "$x" ]] && FOUND_POLICIES+=("$x"); done < <(
    find "${GENES_ROOT}" -maxdepth 4 -type d -path '*/results/05_genes/*' \
        -exec basename {} \; 2>/dev/null | sort -u)
# Restricted to 05_genes: 06_plots holds transcript_status_<scope>_matrix.tsv,
# a plotting artefact with species as rows and transcripts as columns. Matching
# it here would invent a scope called "<scope>_matrix" and, worse, feed a table
# with a completely different layout to anything that counted its columns.
while IFS= read -r x; do [[ -n "$x" ]] && FOUND_SCOPES+=("$x"); done < <(
    find "${GENES_ROOT}" -path '*/results/05_genes/*' \
         -name 'transcript_status_*.tsv' 2>/dev/null \
        | sed 's|.*/transcript_status_||; s|\.tsv$||' | sort -u)
[[ ${#FOUND_POLICIES[@]} -gt 0 ]] || { echo "no per-gene results under ${GENES_ROOT}" >&2; exit 1; }

echo "merging: rankings=${FOUND_POLICIES[*]}  scopes=${FOUND_SCOPES[*]}"
python3 -m morpheus merge-genes --genes-root "${GENES_ROOT}" --outdir "${COMBINED}" \
    --policies "${FOUND_POLICIES[@]}" --scopes "${FOUND_SCOPES[@]}"

# Count the species that actually appear, rather than assuming the whole
# annotation directory was searched -- a --species run legitimately has fewer.
# Read it off the merged table by column name, and drop the human reference row:
# the threshold is about how many query species carry the transcript.
count_species() {
    local f
    for f in "${COMBINED}"/transcript_status_*.tsv; do
        [[ -s "$f" ]] || continue
        awk -F'\t' 'NR==1 { for (i = 1; i <= NF; i++) if ($i == "species") c = i; next }
                    c && $c != "Homo_sapiens" { print $c }' "$f"
    done | sort -u | wc -l | tr -d ' '
}
N_SPECIES="$(count_species)"
if [[ "${N_SPECIES}" -lt 1 ]]; then
    N_SPECIES=$(find "${BAT_ANNOTATION_DIR}" -maxdepth 1 -mindepth 1 -type d ! -name '.*' | wc -l | tr -d ' ')
fi
PLOT_MIN_SPECIES="$(plot_min_species_for "${N_SPECIES}")"
echo "${N_SPECIES} query species; plotting transcripts seen in >= ${PLOT_MIN_SPECIES}"

# The merged tables are the point; a figure that will not draw must not throw
# them away. Report and carry on.
failed=0
draw() { Rscript "$@" || { echo "WARNING: figure failed: Rscript $*" >&2; failed=$((failed + 1)); }; }

if [[ "${SKIP_PLOT}" == "1" ]]; then
    echo "figures skipped (--skip_plot)"
    exit 0
fi

for policy in "${FOUND_POLICIES[@]}"; do
    for scope in "${FOUND_SCOPES[@]}"; do
        st="${COMBINED}/transcript_status_${scope}__${policy}.tsv"
        [[ -s "${st}" ]] || continue
        draw "${R_DIR}/plot_transcript_status.R" \
            --status "${st}" --tree "${SPECIES_TREE}" \
            --outdir "${COMBINED}/plots/${policy}" \
            --min-species "${PLOT_MIN_SPECIES}" \
            --stem "transcript_status_${scope}"
    done
done

cn="${COMBINED}/copy_number_matrix.tsv"
if [[ -s "${cn}" ]]; then
    draw "${R_DIR}/plot_copy_number.R" --matrix "${cn}" \
        --tree "${SPECIES_TREE}" --outdir "${COMBINED}/plots"
fi

[[ ${failed} -eq 0 ]] || echo "WARNING: ${failed} combined figure(s) failed; the merged tables are complete" >&2
echo "combined results in ${COMBINED}"
