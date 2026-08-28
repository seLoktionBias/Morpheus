#!/usr/bin/env bash
###############################################################################
# Collect the HyPhy JSON the array produced into tidy tables, draw the figures,
# and write SUMMARY.md.
#
#   sbatch --dependency=afterany:<hyphy_array_jobid> \
#          /path/to/Morpheus/slurm/04_collect.sbatch
#
# `afterany`, not `afterok`: a single analysis timing out should not stop the
# rest from being collected. Failures are recorded in hyphy_gene_level.tsv.
###############################################################################
#SBATCH --job-name=morpheus_collect
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
source "${MORPHEUS_HOME:?set MORPHEUS_HOME}/slurm/_common.sh"

morpheus_run --only hyphy     # reuses every JSON the array already wrote
morpheus_run --only plots
morpheus_run --only summary
