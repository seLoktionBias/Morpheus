#!/usr/bin/env bash
###############################################################################
# Selection analyses, one array task per alignment.
#
# aBSREL over a hundred taxa can take hours, so the array is the difference
# between a day and a week. Each task runs every requested method on one
# alignment and writes its JSON beside it; a final `--only hyphy` pass collects
# the results, reusing every JSON already on disk.
#
#   n=$(awk -F'\t' 'NR>1 && $16=="1"' results/05_genes/manifest.tsv | wc -l)
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

MANIFEST="${GENES_DIR}/manifest.tsv"
[[ -s "${MANIFEST}" ]] || { echo "missing ${MANIFEST}; run 02_after_search first" >&2; exit 1; }

idx="${SLURM_ARRAY_TASK_ID:?submit this as an array}"
read -r ALIGNMENT TREE GENE TX < <(
  python3 - "${MANIFEST}" "${idx}" <<'PY'
import csv, sys
from pathlib import Path
rows = [r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")
        if r.get("passes_species_threshold") == "1"]
rows.sort(key=lambda r: (r["gene"], r["human_transcript"]))
i = int(sys.argv[2]) - 1
if i < 0 or i >= len(rows):
    raise SystemExit(f"array index {i+1} outside 1..{len(rows)}")
r = rows[i]
fasta = Path(r.get("selected_fasta") or r["fasta"])
aln = fasta.with_name(fasta.name.replace(".cds.fa", ".codon.aln.fa"))
tree = r.get("selected_newick") or r["newick"]
print(aln, tree, r["gene"], r["human_transcript"])
PY
)

echo "[$(date)] task ${idx}: ${GENE} ${TX}"
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
echo "[$(date)] done ${GENE} ${TX}"
