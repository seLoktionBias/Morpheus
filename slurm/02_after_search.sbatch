#!/usr/bin/env bash
###############################################################################
# Everything between the search array and the selection analyses: merge the
# per-species output, screen against the proteome, assign transcripts under both
# scopes, compare them, count copies, build the deliverables, align.
#
#   sbatch --dependency=afterok:<search_array_jobid> \
#          /path/to/Morpheus/slurm/02_after_search.sbatch
#
# The screen is the long pole here (one BLAST of every candidate against ~20k
# human proteins); it is cached against a content hash of the query, so a rerun
# that changes only downstream logic does not repeat it.
###############################################################################
#SBATCH --job-name=morpheus_post_search
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
source "${MORPHEUS_HOME:?set MORPHEUS_HOME}/slurm/_common.sh"

python3 -m morpheus merge-search --outdir "${SEARCH_DIR}"
for step in screen assign compare copy-number deliverables align; do
    morpheus_run --only "${step}"
done
