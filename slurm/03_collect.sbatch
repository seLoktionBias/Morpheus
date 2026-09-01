#!/usr/bin/env bash
###############################################################################
# Draw the figures and write the summary. Last job in the chain.
#
#   sbatch --dependency=afterok:<post_search_jobid> \
#          /path/to/Morpheus/slurm/03_collect.sbatch
#
# With --skip_plot the plots step is skipped and only the summary is written;
# the job still runs, so the chain always ends the same way.
###############################################################################
#SBATCH --job-name=morpheus_collect
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
source "${MORPHEUS_HOME:?set MORPHEUS_HOME}/slurm/_common.sh"

morpheus_run --only plots
morpheus_run --only summary
