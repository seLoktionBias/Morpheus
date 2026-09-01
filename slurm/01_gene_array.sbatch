#!/usr/bin/env bash
###############################################################################
# One array task per gene: the whole per-gene pipeline, in its own results tree.
#
# Each task writes only to <output>/genes/<GENE>/, so tasks cannot corrupt one
# another's tables and nothing needs locking. The one shared thing is the
# flattened GTF cache, written once by 00_reference and read-only here.
#
# Normally submitted for you, with the array size worked out:
#   morpheus run --gene_list genes.txt --path_file paths.txt --mode slurm
#
# Memory and time are per gene, not per list: a gene job reads every query
# genome's annotation, which is the slow part.
###############################################################################
#SBATCH --job-name=morpheus_gene
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail
source "${MORPHEUS_HOME:?set MORPHEUS_HOME}/slurm/_common.sh"

GENES_ROOT="${WORKING_DIR}/genes"
ORDER="${GENES_ROOT}/.gene_order.txt"
# Written by `morpheus run --mode slurm` at submission time, so task N always
# means the same gene no matter how the gene list is edited afterwards.
[[ -s "${ORDER}" ]] || { echo "missing ${ORDER}; submit with 'morpheus run --mode slurm'" >&2; exit 1; }

idx="${SLURM_ARRAY_TASK_ID:?this script must be submitted as an array}"
GENE="$(sed -n "${idx}p" "${ORDER}")"
[[ -n "${GENE}" ]] || { echo "array index ${idx} is past the end of ${ORDER}" >&2; exit 1; }

echo "[$(date)] task ${idx}: ${GENE}"
morpheus_run --per_gene no --gene "${GENE}" --output "${GENES_ROOT}/${GENE}" \
    --bat_dir "${BAT_ANNOTATION_DIR}" --ref_dir "${HUMAN_GENOME_DIR}" \
    --tree "${SPECIES_TREE}"
echo "[$(date)] done ${GENE}"
