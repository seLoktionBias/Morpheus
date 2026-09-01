# Changelog

## 2.1.0

### Changed

- **`synteny_aware` is renamed `structure_aware`.** The old name was a misnomer:
  the ranking weighs *exon structure* and nothing positional. No scoring term
  reads a coordinate, a flanking gene, or a distance from the home locus —
  `structure` is the longest common subsequence of exon labels, and
  `upstream_gene` / `downstream_gene` / `distance_to_home_locus_bp` appear in the
  output tables as recorded columns only, never in the score.

  The positional constraint on orthology is the *other* axis, the
  `region_restricted` scope. Because of that, `synteny_aware/unrestricted` read
  as a contradiction — "judge by gene order, but ignore where it is" — for a
  combination that is neither of those things. `structure_aware/unrestricted`
  says what it does.

  Scope names are unchanged. Nothing about the computation changed: the same
  candidates get the same scores and the same assignments. Only the label moved.

  Output from earlier versions keeps the old directory name; the README shows the
  three `mv` commands that bring an existing `results/` in line without re-running.


## 2.0.0

Breaking. The pipeline now ends at the codon alignment, and the command line is
explicit rather than convention-driven.

### Removed

- **Selection analysis.** The `hyphy` step, `morpheus/hyphy.py`, the HyPhy Slurm
  array and every `HYPHY_*` setting are gone, along with HyPhy itself from
  `environment.yml`. Model choice, branch labelling and multiple-testing
  correction carry the scientific claim and belong to whoever is making it. Each
  `gene_files_selected/` directory still holds an alignment and a pruned tree
  with identical taxa, so any codon-model tool can be pointed at it directly.
- `results/06_selection/`.

### Changed

- **Results directories renumbered `01`–`06`, contiguous.** `07_plots/` is now
  `06_plots/`. An existing `results/` from 1.x keeps its old `06_selection/` and
  `07_plots/`; both can be deleted.
- Slurm scripts renumbered `00`–`03`; `04_collect.sbatch` is now
  `03_collect.sbatch`. All four now share the `_common.sh` preamble.
- `.not_run_yet` markers are swept again at the end of a run, so a directory
  filled during that run no longer keeps a marker saying it is empty.

### Added

- **Explicit flags** on `morpheus run`: `--gene`, `--gene_list`, `--mode`,
  `--search`, `--tree`, `--bat_dir`, `--ref_dir`, `--output`, `--path_file`,
  `--skip_plot`. Bad values are refused rather than coerced to a default.
- **`--mode slurm`** submits the whole dependency chain in one command and prints
  the job IDs. Possible now that the HyPhy array — whose size was not knowable
  until mid-run — is gone. `--dry-run` works without `sbatch` present.
- **`--search region|similarity|both`** selects the scope and ranking together;
  `region` and `similarity` each run about half the work of `both`.
- **`--skip_plot`** skips every figure, and drops R from the requirements.
- **`--per_gene yes|no`, default `yes`.** Each gene now runs on its own: one
  Slurm job per gene, or one gene at a time locally so a long list cannot swamp
  a workstation. Each gene writes a complete results tree under
  `<output>/genes/<GENE>/`, and `morpheus merge-genes` concatenates them into
  `<output>/combined/` for the whole-list tables and figures. Nothing is shared
  but the flattened GTF cache, which is read-only after the reference job, so
  concurrent gene jobs cannot corrupt a common table and one gene failing leaves
  the rest intact. `--per_gene no` keeps the previous whole-list behaviour,
  which reads each query genome once rather than once per gene.
- **`morpheus merge-genes`** and `scripts/merge_and_plot.sh`, shared by the local
  and Slurm paths so the merge happens one way, not two. Concatenation is by
  column name: a gene missing a column gets `NA` there rather than the next
  column's value sliding into its place.
- **`--env_mode inherit|conda|mamba|venv`, `--env_name` and `--venv_path`**,
  matching Pensieve's model, so a job can resolve its own environment without
  relying on `conda activate` working non-interactively.
- **`--time SPEC`** — hours, `HH:MM:SS` or `D-HH:MM:SS` — applied to all four
  jobs. Per-job defaults stay in the sbatch files, since the jobs differ.
- **Slurm site options**: `--slurm_partition`, `--slurm_account`,
  `--slurm_module`, `--array_throttle` and `--slurm_extra`. There was previously
  no way to name a queue, so `--mode slurm` could not be used at all on a cluster
  whose default partition is not the right one. Each is warned about if passed
  without `--mode slurm`.
- `paths.txt` accepts key synonyms, comments, quotes and `~`, and **reports an
  unrecognised key as an error** instead of ignoring the typo.
- `--path_file` treats the file's directory as the project, so `gene_list.txt`
  and the tree resolve next to it.
- **`install.sh` writes a `morpheus` launcher into the environment**, so
  `mamba activate Morpheus` is all that is needed to get the command on PATH —
  in an interactive session and in a job script alike. Previously the checkout
  had to be added to PATH by hand, or called by full path. The launcher is a
  shim, not a symlink: it records the checkout's absolute path and, if that
  checkout later moves, says so instead of failing with a bare ENOENT.
  `--no-launcher` opts out.
- `morpheus env` reports every resolved input path with `ok` / `MISSING`.
- `examples/paths.txt` and `examples/gene_list.txt`.
- **Worked examples in `morpheus --help` and `morpheus run --help`**: a numbered
  first-run walkthrough including what `paths.txt` contains, and complete
  commands for the local, cluster, no-paths.txt and resume-after-failure cases.
  `tests/check_help.py` holds them to the parser — every flag shown must be
  accepted, every flag accepted must be shown.
- `tests/check_portability.py`: greps every shell file for constructs that work
  on one bash and not another — `${#arr[@]-default}`, `mapfile`/`readarray`,
  `${var^^}`, associative arrays. `bash -n` catches none of these, because the
  first is a runtime expansion error and the rest are missing builtins.
- `tests/check_numbering.py` and `tests/check_slurm.py`: contiguous numbering,
  shared preamble, valid job names and log paths, no Slurm job asking for a step
  the pipeline does not have, and `--mode slurm` submitting exactly the scripts
  that exist), launcher generation, help/parser agreement, and flag rejection.
  Smoke test now 60 checks, including nine on Slurm environment failure modes
  and the per-gene merge's column handling.

### Fixed

- **`${#ARR[@]-0}` is invalid** — the `#` length operator cannot take a
  `-default`. macOS bash 3.2 tolerated it; bash on a RHEL-family cluster
  rejects it, so `00_reference` died on `bad substitution` and every dependent
  job sat in `DependencyNeverSatisfied`. Replaced with `${#ARR[@]}`, which is
  safe under `set -u` for an already-declared empty array.
- **`--per_gene` was never propagated to the submitted jobs**, so each job
  defaulted back to `yes` and re-entered the per-gene branch — under
  `--per_gene no` the reference job tried to orchestrate the whole gene list
  instead of flattening the GTF. The two near-identical `--export` blocks that
  let this drift are now one `build_exports` function, and `slurm/_common.sh`
  always passes `--per_gene no`: whether to split by gene is decided once, at
  submit time.
- `04_collect` chose its behaviour by testing whether `genes/` existed, so a
  directory left by an earlier `--per_gene yes` run would hijack a later
  `--per_gene no` one. It reads `MORPHEUS_PER_GENE` now.
- The smoke test used a fixed scratch directory, so two concurrent runs deleted
  each other's and the loser reported failures unrelated to the code.

- **A mamba-only install was reported as "no conda installation found".**
  `activate_env` searched for a `conda` *binary* to locate the environment. A
  `--backend=mamba` install has no reason to provide one and micromamba never
  does. It now looks for the environment directory itself, and knows about
  `MAMBA_ROOT_PREFIX` and `~/micromamba`.
- **`PLOT_MIN_SPECIES` was a fixed 50** while `MIN_SPECIES_FRACTION` was a
  fraction, so any run over fewer than 50 genomes produced an empty transcript
  status figure — or, with the fix below not yet in place, a failed run. It is
  now derived from `MIN_SPECIES_FRACTION` of the species actually searched
  (rounded up, minimum 1), matching the threshold that already gated
  `gene_files_selected`. At the full 103 genomes that is 52, so existing results
  are unaffected; setting `PLOT_MIN_SPECIES` explicitly still overrides it.
- The merge step counted species with a `find` that also matched
  `06_plots/transcript_status_<scope>_matrix.tsv` — a plotting artefact with
  species as rows and transcripts as columns. It reported 104 species for a
  4-species run, and would have invented a scope named `<scope>_matrix`.
- **A figure that would not draw failed the whole run.** Under `--per_gene` that
  surfaced as the *gene* having failed, when in fact every table was written and
  only a threshold (`PLOT_MIN_SPECIES`) was too high for the species subset.
  Figure failures now warn and are counted; the run still succeeds.
- `mapfile` is bash 4+, and macOS ships bash 3.2 — local runs died immediately
  on `mapfile: command not found`. Replaced with portable while-read loops.
- A cleanup `trap` in the environment-mode test could overwrite the tests' own
  exit status, so a fully passing run could report a failure.

- **Slurm jobs failed with no usable error.** `slurm/_common.sh` ran
  `module load miniforge` and `conda activate Morpheus`, both with `|| true`,
  and both are wrong by default: sites name modules differently, and
  `conda activate` does nothing useful in a non-interactive batch shell. The
  environment half-activated, the job died several steps later complaining about
  a missing file, and the real cause appeared nowhere — surfacing as
  `DependencyNeverSatisfied` on every downstream job in the chain.

  Nothing is guessed now. `--env_mode inherit` is the default and takes the
  submitting shell's PATH; `conda`/`mamba` run the work through
  `<manager> run -n`, which works without an interactive shell; `--slurm_module`
  is opt-in and errors if `module` is unavailable. Every job verifies its tools
  before doing anything and, if they are missing, fails at the top of the log
  naming them and the flags to fix it.
- `conda run` buffers output until the process exits, which would leave a
  twelve-hour job's log empty while it ran. `--no-capture-output` is added when
  the manager supports it, with a note when it does not.

- `morpheus --help` printed `####` and `set -euo pipefail` at the user: it was
  built with `sed -n '3,18p'` over the script's own comment block, and that line
  range drifted as comments were added above it. Both help screens are now
  explicit text, and a test fails if shell source ever leaks again.
- The HyPhy Slurm array read `results/05_genes/manifest.tsv`, a path the
  per-policy split had removed; every task would have failed. Moot now that the
  array is gone, but it was live in 1.0.0.

## 1.0.0

First release.
