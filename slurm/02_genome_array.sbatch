#!/usr/bin/env bash
###############################################################################
# Region-restricted search, one array task per query genome.
#
# This is the only embarrassingly parallel step: each genome is read
# independently, so an array turns ~2 minutes x N genomes into ~2 minutes.
#
# Normally submitted for you, with the array size already worked out:
#   morpheus run --mode slurm --gene_list genes.txt --path_file paths.txt
#
# By hand, from the project directory:
#   n=$(find "$(awk -F= '/bat_annotation_dir/{print $2}' paths.txt)" \
#         -maxdepth 1 -mindepth 1 -type d ! -name '.*' | wc -l)
#   sbatch --array=1-$n $MORPHEUS_HOME/slurm/02_genome_array.sbatch
#
# Used only by --per_gene no. 00_reference must have run first. Array tasks
# pass --no-merge and 03_after_search merges them; merging inside each task would
# race and leave partial combined tables.
###############################################################################
#SBATCH --job-name=morpheus_search
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail
source "${MORPHEUS_HOME:?set MORPHEUS_HOME}/slurm/_common.sh"

# One genome per array index, in a stable sorted order so a rerun of task N
# always picks the same genome.
GENOMES=()
while IFS= read -r _d; do GENOMES+=("${_d}"); done < <(
    find "${BAT_ANNOTATION_DIR}" -maxdepth 1 -mindepth 1 -type d ! -name '.*' | sort)
idx="${SLURM_ARRAY_TASK_ID:?this script must be submitted as an array}"
if (( idx < 1 || idx > ${#GENOMES[@]} )); then
    echo "array index ${idx} outside 1..${#GENOMES[@]}" >&2
    exit 1
fi
GENOME="${GENOMES[$((idx - 1))]}"
SPECIES="$(basename "${GENOME}" | awk -F'__' '{print $1}')"

echo "[$(date)] task ${idx}/${#GENOMES[@]}  ${SPECIES}"
morpheus_exec python3 -m morpheus bat-search \
    --annotations "${BAT_ANNOTATION_DIR}" \
    --reference "${REFERENCE_DIR}" \
    --families "${REFERENCE_DIR}/gene_families.tsv" \
    --tree "${SPECIES_TREE}" \
    --outdir "${SEARCH_DIR}" \
    --locus-slop "${LOCUS_SLOP}" \
    --species "${SPECIES}" \
    --no-merge \
    --threads "${SLURM_CPUS_PER_TASK:-2}"
echo "[$(date)] done ${SPECIES}"
