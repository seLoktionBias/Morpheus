#!/usr/bin/env bash
###############################################################################
# Steps 1-3: flatten the Ensembl GTF, build the human reference, define the
# paralog families. Run once; the GTF table is cached and reused thereafter.
#
#   sbatch /path/to/Morpheus/slurm/00_reference.sbatch
#
# Memory is set for the GTF parse, which is the heaviest part of the pipeline
# on a login node and the usual reason a naive run gets killed there.
###############################################################################
#SBATCH --job-name=morpheus_reference
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
source "${MORPHEUS_HOME:?set MORPHEUS_HOME}/slurm/_common.sh"

morpheus_run --only cds-table
morpheus_run --only human-reference
morpheus_run --only families
