#!/usr/bin/env bash
###############################################################################
# Selection analyses, one array task per alignment.
#
# aBSREL over a hundred taxa can take hours, so the array is the difference
# between a day and a week. Each task runs every requested method on one
# alignment and writes its JSON beside it; a final `--only hyphy` pass collects
# the results, reusing every JSON already on disk.
#
# The task list covers every ranking policy, so one array does both. Count the
# tasks with the helper rather than by hand:
#
#   n=$(bash /path/to/Morpheus/slurm/03_hyphy_array.sbatch --count)
#   sbatch --array=1-${n}%20 /path/to/Morpheus/slurm/03_hyphy_array.sbatch
#   sbatch --dependency=afterany:<array_jobid> /path/to/Morpheus/slurm/04_collect.sbatch
#
# %20 throttles to 20 concurrent tasks; raise or lower to suit the queue.
###############################################################################
#SBATCH --job-name=morpheus_hyphy
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail
source "${MORPHEUS_HOME:?set MORPHEUS_HOME}/slurm/_common.sh"

# Enumerate (policy, gene, transcript) triples in a stable order. Emitting one
# line per task keeps the counting mode and the running mode reading the same
# list, so `--count` can never disagree with what the array indexes into.
task_list() {
    python3 - "${GENES_DIR}" "${POLICIES}" <<'PY'
import csv, sys
from pathlib import Path
genes_dir, policies = Path(sys.argv[1]), sys.argv[2].split()
for policy in policies:
    manifest = genes_dir / policy / "manifest.tsv"
    if not manifest.is_file():
        continue
    rows = [r for r in csv.DictReader(manifest.open(), delimiter="\t")
            if r.get("passes_species_threshold") == "1"]
    for r in sorted(rows, key=lambda r: (r["gene"], r["human_transcript"])):
        fasta = Path(r.get("selected_fasta") or r["fasta"])
        aln = fasta.with_name(fasta.name.replace(".cds.fa", ".codon.aln.fa"))
        tree = r.get("selected_newick") or r["newick"]
        print("\t".join([str(aln), tree, r["gene"], r["human_transcript"], policy]))
PY
}

if [[ "${1-}" == "--count" ]]; then task_list | wc -l; exit 0; fi
if [[ "${1-}" == "--list" ]];  then task_list | cat -n;  exit 0; fi

idx="${SLURM_ARRAY_TASK_ID:?submit this as an array, or pass --count}"
LINE="$(task_list | sed -n "${idx}p")"
[[ -n "${LINE}" ]] || { echo "array index ${idx} outside the task list" >&2; exit 1; }
IFS=$'\t' read -r ALIGNMENT TREE GENE TX POLICY <<<"${LINE}"

echo "[$(date)] task ${idx}: ${GENE} ${TX} (${POLICY})"
[[ -s "${ALIGNMENT}" ]] || { echo "missing alignment ${ALIGNMENT}; run 02_after_search first" >&2; exit 1; }
WORK="$(dirname "${ALIGNMENT}")/hyphy"
mkdir -p "${WORK}"

python3 - "${ALIGNMENT}" "${TREE}" "${WORK}" "${SLURM_CPUS_PER_TASK:-4}" "${HYPHY_TIMEOUT:-}" <<PY
import sys
from pathlib import Path
from morpheus import hyphy
aln, tree, work, threads, timeout = sys.argv[1:6]
prepared = hyphy._tree_for_alignment(Path(aln), Path(tree), Path(work) / "tree_for_hyphy.nwk")
if prepared is None:
    raise SystemExit("alignment and tree taxa do not match")
for method in "${HYPHY_METHODS}".split():
    hyphy.run_one(Path(aln), prepared, Path(work), method,
                  threads=int(threads), timeout=int(timeout) if timeout else None)
PY
echo "[$(date)] done ${GENE} ${TX} (${POLICY})"
