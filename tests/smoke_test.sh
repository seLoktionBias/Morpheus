#!/usr/bin/env bash
###############################################################################
# Morpheus smoke test. Exercises the parts that do not need genome data: the
# optimal matcher against brute force, the codon table, the ORF report, Newick
# pruning, the CLI surface, and every R figure script against synthetic input.
#
#   bash tests/smoke_test.sh
###############################################################################
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Unique per run: a fixed path means two concurrent invocations rm -rf each
# other's scratch mid-test, and the loser reports a failure that has nothing to
# do with the code.
TMP="${HERE}/tests/tmp_smoke.$$"
rm -rf "${TMP}"; mkdir -p "${TMP}"
trap 'rm -rf "${TMP}"' EXIT
# The `bash -c` checks below run in child shells; without exporting these they
# expand to empty there and a check can pass for entirely the wrong reason.
export HERE TMP
export PYTHONPATH="${HERE}${PYTHONPATH:+:${PYTHONPATH}}"

pass=0; fail=0
check() { if "$@" >>"${TMP}/log" 2>&1; then echo "  PASS  $*"; pass=$((pass+1));
          else echo "  FAIL  $*"; fail=$((fail+1)); fi; }

echo "== Morpheus smoke test =="

echo "-- python core --"
check python3 - <<'PY'
import itertools, random
from morpheus.matching import max_weight_matching
random.seed(1)
for _ in range(200):
    r, c = random.randint(1, 6), random.randint(1, 6)
    P = [[round(random.random(), 3) for _ in range(c)] for _ in range(r)]
    pairs = max_weight_matching(P, min_profit=-1)
    got = sum(P[i][j] for i, j in pairs)
    k = min(r, c)
    best = max(sum(P[a][b] for a, b in zip(rows, cols))
               for rows in itertools.permutations(range(r), k)
               for cols in itertools.permutations(range(c), k))
    assert abs(got - best) < 1e-9
    assert len({i for i, _ in pairs}) == len(pairs)
    assert len({j for _, j in pairs}) == len(pairs)
PY

check python3 - <<'PY'
from morpheus.common import translate, revcomp, orf_report, clean_nt
assert translate("ATGGTGCATCTGACTCCTGAGGAGAAGTAA") == "MVHLTPEEK*"
assert revcomp("ATGCN") == "NGCAT"
assert clean_nt("aug u") == "ATGT"
r = orf_report("ATG" + "AAA" * 750 + "TGA" + "AAAA")
assert r["orf_fraction"] > 0.95 and r["has_start"]
r = orf_report("ATG" + "AAA" * 10 + "TGA" + "AAA" * 88 + "TAA")
assert r["orf_fraction"] < 0.2
PY

check python3 - <<'PY'
from morpheus import newick
t = newick.parse("((A:1,B:2)n1:3,(C:4,(D:5,E:6)n2:7)n3:8)root;")
assert sorted(t.leaf_names()) == list("ABCDE")
assert newick.write(newick.prune(t, {"A", "C"})) == "(A:4,C:12)root;"
assert newick.prune(t, {"A"}) is not None
PY

check python3 - <<'PY'
from morpheus.deliverables import cds_status
from morpheus.copy_number import _cds_integrity
assert cds_status("ATG" + "AAA" * 10 + "TGA") == "complete"
assert cds_status("ATG" + "AAA" * 750 + "TGA" + "AAAA") == "complete"
assert cds_status("ATG" + "AAA" * 10 + "TGA" + "AAA" * 88 + "TAA") == "pseudogenized"
assert _cds_integrity("ATG" + "AAA" * 10 + "TGA") == "complete_orf"
PY

check python3 - <<'PY'
from morpheus.pairwise import structure_similarity, _lcs_length
assert _lcs_length(list("abcd"), list("axcd")) == 3
sm = structure_similarity(["exon1", "exon2"], ["exon1", "exon2", "exon3"])
assert 0.6 < sm["structure_lcs_coverage"] < 0.7
PY

echo "-- CLI --"
check python3 -m morpheus --help
for step in cds-table human-reference families bat-search merge-search screen \
            assign compare copy-number deliverables align merge-genes summary; do
    check python3 -m morpheus "${step}" --help
done
check bash "${HERE}/bin/morpheus" version
check bash "${HERE}/bin/morpheus" help
check bash "${HERE}/scripts/run_pipeline.sh" --help

# Bad flag values must be refused, not silently coerced into a default that
# quietly runs the wrong analysis.
refuses() {  # label, then the args to pass to run_pipeline
    local label="$1"; shift
    if bash "${HERE}/scripts/run_pipeline.sh" "$@" >>"${TMP}/log" 2>&1; then
        echo "  FAIL  refuses ${label}"; fail=$((fail+1))
    else
        echo "  PASS  refuses ${label}"; pass=$((pass+1))
    fi
}
refuses "--mode nonsense"         --mode nonsense --dry-run
refuses "--search nonsense"       --search nonsense --dry-run
refuses "--only bogus-step"       --only bogus-step --dry-run
refuses "--from bogus-step"       --from bogus-step --dry-run
refuses "an unknown flag"         --not-a-flag --dry-run
refuses "--gene with --gene_list" --gene A --gene_list /dev/null --dry-run
refuses "a missing --path_file"   --path_file /nonexistent/paths.txt --dry-run
refuses "a non-numeric --array_throttle" --mode slurm --array_throttle lots --dry-run
refuses "--per_gene maybe"        --per_gene maybe --dry-run
refuses "--env_mode conda-forge"  --env_mode conda-forge --dry-run
refuses "--env_mode venv with no --venv_path" --env_mode venv --dry-run
refuses "--time 2h"               --mode slurm --time 2h --dry-run

echo "-- selection analysis is gone --"
check bash -c '! grep -rniE "hyphy|absrel|busted" \
    --include="*.py" --include="*.sh" --include="*.sbatch" --include="*.yml" \
    --include="*.R" "${HERE}/morpheus" "${HERE}/scripts" "${HERE}/slurm" \
    "${HERE}/config" "${HERE}/R" "${HERE}/bin" "${HERE}/environment.yml"'
check bash -c '[[ ! -e "${HERE}/morpheus/hyphy.py" ]]'
check bash -c '! python3 -m morpheus hyphy --help 2>/dev/null'

echo "-- paths.txt --"
check bash -c '
    d="${TMP}/paths_ok"; mkdir -p "$d"
    printf "%s\n" "# a comment" "HUMAN_DIR = /tmp/hg   # trailing" \
        "bat_dir=\"/tmp/bats\"" "outdir=/tmp/out" > "$d/paths.txt"
    cd "$d"
    source "${HERE}/config/config.sh"
    [[ "${HUMAN_GENOME_DIR}"   == /tmp/hg   ]] || { echo "ref  wrong: ${HUMAN_GENOME_DIR}";   exit 1; }
    [[ "${BAT_ANNOTATION_DIR}" == /tmp/bats ]] || { echo "bat  wrong: ${BAT_ANNOTATION_DIR}"; exit 1; }
    [[ "${WORKING_DIR}"        == /tmp/out  ]] || { echo "work wrong: ${WORKING_DIR}";        exit 1; }'
# An unknown key is a typo, and a typo that is silently ignored becomes a run
# against the wrong genome directory.
check bash -c '
    d="${TMP}/paths_bad"; mkdir -p "$d"; echo "hooman_dir=/tmp/x" > "$d/paths.txt"
    cd "$d"
    ! ( set -euo pipefail; source "${HERE}/config/config.sh" ) 2>/dev/null'
check bash -c '
    d="${TMP}/paths_example"; mkdir -p "$d"
    cp "${HERE}/examples/paths.txt" "$d/paths.txt"; cd "$d"
    set -euo pipefail; source "${HERE}/config/config.sh"'

echo "-- env launcher --"
# `conda activate Morpheus` must be enough to get `morpheus` on PATH, so the
# installer has to write a working shim into the environment's bin/.
check bash -c '
    root="${TMP}/fakeconda"; mkdir -p "${root}/envs/SmokeEnv/bin"
    MORPHEUS_ENV=SmokeEnv MORPHEUS_CONDA_ROOT="${root}" \
        bash "${HERE}/install.sh" --env-only --backend=current --no-launcher >/dev/null 2>&1
    [[ ! -e "${root}/envs/SmokeEnv/bin/morpheus" ]] || { echo "--no-launcher still wrote one"; exit 1; }
    MORPHEUS_ENV=SmokeEnv MORPHEUS_CONDA_ROOT="${root}" \
        bash "${HERE}/install.sh" --env-only --backend=current >/dev/null 2>&1
    shim="${root}/envs/SmokeEnv/bin/morpheus"
    [[ -x "${shim}" ]]            || { echo "no launcher written"; exit 1; }
    bash -n "${shim}"             || { echo "launcher is not valid bash"; exit 1; }
    [[ "$("${shim}" version)" == "$(cat "${HERE}/VERSION")" ]] \
        || { echo "launcher does not reach the checkout"; exit 1; }
    # A moved checkout must say so, not fail with a bare ENOENT.
    sed -i.bak "s|MORPHEUS_HOME=\".*\"|MORPHEUS_HOME=\"/gone\"|" "${shim}"
    out="$("${shim}" version 2>&1)"; rc=$?
    [[ ${rc} -eq 127 ]]           || { echo "stale launcher exit was ${rc}, wanted 127"; exit 1; }
    grep -q "/gone" <<<"${out}"   || { echo "stale launcher did not name the missing checkout"; exit 1; }'

echo "-- help matches the program --"
check python3 "${HERE}/tests/check_help.py" "${HERE}"

echo "-- shell portability --"
check python3 "${HERE}/tests/check_portability.py" "${HERE}"
# The idiom the banned one was standing in for must actually be safe.
check bash -c 'set -euo pipefail; A=(); [[ ${#A[@]} -eq 0 ]] && printf "%s" "${A[@]+"${A[@]}"}"'

echo "-- results numbering is contiguous --"
check python3 "${HERE}/tests/check_numbering.py" "${HERE}"

check python3 - <<'PY'
# Per-gene runs are merged by column name. A gene whose table lacks a column
# must get NA there, not the next column's value sliding into its place.
import tempfile
from pathlib import Path
from morpheus.merge_genes import _concat
from morpheus.common import read_tsv
d = Path(tempfile.mkdtemp())
(d / "a.tsv").write_text("gene\tx\ty\nMX1\t1\t2\n")
(d / "b.tsv").write_text("y\tgene\tz\n9\tOAS1\t7\n")
out = _concat([d / "a.tsv", d / "b.tsv"], d / "m.tsv")
rows = read_tsv(out)
assert [r["gene"] for r in rows] == ["MX1", "OAS1"], rows
assert rows[0]["y"] == "2" and rows[1]["y"] == "9", "column order not honoured"
assert rows[0]["z"] == "NA" and rows[1]["x"] == "NA", "missing column should be NA"
assert _concat([d / "nope.tsv"], d / "empty.tsv") is None
PY

echo "-- slurm environment modes --"
check bash "${HERE}/tests/check_env_modes.sh"

echo "-- slurm scripts --"
check python3 "${HERE}/tests/check_slurm.py" "${HERE}"

echo "-- shell syntax --"
check bash -n "${HERE}/scripts/run_pipeline.sh"
check bash -n "${HERE}/config/config.sh"
check bash -n "${HERE}/install.sh"
check bash -n "${HERE}/bin/morpheus"
for f in "${HERE}"/slurm/*.sbatch "${HERE}"/slurm/_common.sh; do check bash -n "$f"; done

echo "-- R figures on synthetic input --"
if command -v Rscript >/dev/null 2>&1; then
    python3 - "${TMP}" <<'PY'
import random, sys
from pathlib import Path
tmp = Path(sys.argv[1]); random.seed(2)
tips = [f"Sp{i}" for i in range(1, 21)]
Path(tmp / "tree.nwk").write_text("(" + ",".join(["Homo_sapiens"] + tips) + ");\n")

with open(tmp / "status.tsv", "w") as f:
    f.write("gene\thuman_transcript\tspecies\tstatus\n")
    for g in ("GENEA", "GENEB"):
        for t in ("ENST1", "ENST2"):
            for s in ["Homo_sapiens"] + tips:
                f.write(f"{g}\t{t}\t{s}\t{random.choice(['complete','partial','fragmented'])}\n")

with open(tmp / "copy.tsv", "w") as f:
    f.write("species\tgene\ttotal_copies\tunambiguous_copies\tshared_copies\t"
            "copies_excluding_retro\tfunctional_copies\tretro_copies\n")
    for s in tips:
        for g in ("GENEA", "GENEB"):
            n = random.choice([0, 1, 1, 2, 3, 14])
            f.write(f"{s}\t{g}\t{n}\t{n}\t0\t{n}\t{n}\t0\n")

with open(tmp / "scope.tsv", "w") as f:
    f.write("species\tgene\thuman_transcript\toutcome\tregion_cds_status\t"
            "unrestricted_cds_status\n")
    for g in ("GENEA", "GENEB"):
        for t in ("ENST1", "ENST2"):
            for s in tips:
                f.write(f"{s}\t{g}\t{t}\t{random.choice(['SAME_MODEL','DIFFERENT_LOCUS'])}"
                        f"\tfragmented\tcomplete\n")

with open(tmp / "exons.tsv", "w") as f:
    f.write("gene\ttranscript_id\tis_longest\tchrom\tstrand\tstart\tend\t"
            "transcript_exon_rank\texon_label\texon_label_type\n")
    for ti, t in enumerate(("ENST1", "ENST2")):
        pos = 1000
        for r in range(1, 6):
            f.write(f"GENEA\t{t}\t{int(ti==0)}\t1\t-\t{pos}\t{pos+120}\t{r}\t"
                    f"exon{r}\tlongest\n")
            pos += 4000
PY
    R="${HERE}/R"
    check Rscript "${R}/plot_transcript_status.R" --status "${TMP}/status.tsv" \
        --tree "${TMP}/tree.nwk" --outdir "${TMP}/out" --min-species 5 --format png
    check Rscript "${R}/plot_copy_number.R" --matrix "${TMP}/copy.tsv" \
        --tree "${TMP}/tree.nwk" --outdir "${TMP}/out"
    check Rscript "${R}/plot_scope_comparison.R" --comparison "${TMP}/scope.tsv" \
        --tree "${TMP}/tree.nwk" --outdir "${TMP}/out" --min-species 5
    check Rscript "${R}/plot_exon_models.R" --plot-tsv "${TMP}/exons.tsv" \
        --context-tsv /dev/null --gene GENEA --outdir "${TMP}/out" --format png
else
    echo "  SKIP  R figures (Rscript not found)"
fi

echo
echo "== ${pass} passed, ${fail} failed =="
if [[ ${fail} -ne 0 ]]; then
    # The scratch dir is about to be removed, so keep the log where it can be read.
    cp "${TMP}/log" "${HERE}/tests/smoke_failure.log" 2>/dev/null || true
    echo "see ${HERE}/tests/smoke_failure.log"
    exit 1
fi
