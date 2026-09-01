# Changelog

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
- `paths.txt` accepts key synonyms, comments, quotes and `~`, and **reports an
  unrecognised key as an error** instead of ignoring the typo.
- `--path_file` treats the file's directory as the project, so `gene_list.txt`
  and the tree resolve next to it.
- `morpheus env` reports every resolved input path with `ok` / `MISSING`.
- `examples/paths.txt` and `examples/gene_list.txt`.
- `tests/check_numbering.py` and `tests/check_slurm.py`: contiguous numbering,
  shared preamble, valid job names and log paths, no Slurm job asking for a step
  the pipeline does not have, and `--mode slurm` submitting exactly the scripts
  that exist. Smoke test now 49 checks.

### Fixed

- The HyPhy Slurm array read `results/05_genes/manifest.tsv`, a path the
  per-policy split had removed; every task would have failed. Moot now that the
  array is gone, but it was live in 1.0.0.

## 1.0.0

First release.
