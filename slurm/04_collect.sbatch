#!/usr/bin/env bash
###############################################################################
# Last job in either chain: merge what ran, draw the whole-list figures, and
# write the summary.
#
# With --per_gene yes it concatenates every gene's tables into <output>/combined/
# and plots those. With --per_gene no there is nothing to merge and it plots the
# single shared results tree.
#
# Submitted with afterany, not afterok: one gene failing must not stop the other
# twenty from being merged and plotted. merge-genes names any gene that produced
# nothing, so a partial result is never mistaken for a complete one.
###############################################################################
#SBATCH --job-name=morpheus_collect
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
source "${MORPHEUS_HOME:?set MORPHEUS_HOME}/slurm/_common.sh"

GENES_ROOT="${WORKING_DIR}/genes"
if [[ -d "${GENES_ROOT}" ]]; then
    morpheus_exec bash "${MORPHEUS_HOME}/scripts/merge_and_plot.sh" \
        "${GENES_ROOT}" "${WORKING_DIR}/combined"
else
    morpheus_run --only plots
    morpheus_run --only summary
fi
