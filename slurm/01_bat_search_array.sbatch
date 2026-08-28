#!/usr/bin/env bash
###############################################################################
# Region-restricted search, one array task per query genome.
#
# This is the only embarrassingly parallel step: each genome is read
# independently, so an array turns ~2 minutes x N genomes into ~2 minutes.
#
#   cd /path/to/your/project          # where paths.txt and gene_list.txt live
#   sbatch --array=1-$(ls -d "$(awk -F= '/bat_annotation_dir/{print $2}' paths.txt)"/*/ | wc -l) \
#          /path/to/Morpheus/slurm/01_bat_search_array.sbatch
#
# Steps 1-3 (cds-table, human-reference, families) must have run first; see
# 00_reference.sbatch. Merge the array output with 02_after_search.sbatch.
###############################################################################
#SBATCH --job-name=morpheus_search
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

# Apocrita: the conda module is `miniforge`. Override MORPHEUS_MODULE for a
# site that names it differently, or leave it empty if conda is already on PATH.
MODULE="${MORPHEUS_MODULE-miniforge}"
if [[ -n "${MODULE}" ]] && command -v module >/dev/null 2>&1; then
    module load "${MODULE}" || true
fi
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${MORPHEUS_ENV:-Morpheus}"

MORPHEUS_HOME="${MORPHEUS_HOME:?set MORPHEUS_HOME to the Morpheus checkout}"
PROJECT_DIR="${MORPHEUS_PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "${PROJECT_DIR}"
mkdir -p logs

source "${MORPHEUS_HOME}/config/config.sh"
export PYTHONPATH="${MORPHEUS_HOME}${PYTHONPATH:+:${PYTHONPATH}}"

# One genome per array index, in a stable sorted order so a rerun of task N
# always picks the same genome.
mapfile -t GENOMES < <(find "${BAT_ANNOTATION_DIR}" -maxdepth 1 -mindepth 1 -type d \
                       ! -name '.*' | sort)
idx="${SLURM_ARRAY_TASK_ID:?this script must be submitted as an array}"
if (( idx < 1 || idx > ${#GENOMES[@]} )); then
    echo "array index ${idx} outside 1..${#GENOMES[@]}" >&2
    exit 1
fi
GENOME="${GENOMES[$((idx - 1))]}"
SPECIES="$(basename "${GENOME}" | awk -F'__' '{print $1}')"

echo "[$(date)] task ${idx}/${#GENOMES[@]}  ${SPECIES}"
python3 -m morpheus bat-search \
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
