# Morpheus

Comparative transcript recovery, orthologous transcript assignment, and gene copy
number across TOGA2-annotated genomes.

Give it a list of human genes, a species tree, and a directory of TOGA2 results —
one per genome — and it recovers each gene's transcripts in every genome, decides
which query transcript corresponds to which human transcript, counts gene copies,
and assembles ready-to-analyse codon alignments with matching pruned trees.

Morpheus stops at the alignment. It deliberately runs no selection tests: which
model, which branches, and which correction to apply are decisions that belong to
whoever is asking the question, and the output is laid out so that any codon-model
tool can be pointed straight at it.

Everything runs locally or on Slurm from the same command. The pipeline logic uses
only the Python standard library; the external tools are BLAST, MAFFT and R.

---

## Contents

- [Why the answers are not obvious](#why-the-answers-are-not-obvious)
- [Install](#install)
- [Quick start](#quick-start)
- [Inputs](#inputs)
- [Two search scopes](#two-search-scopes)
- [Two ranking policies](#two-ranking-policies)
- [Method, step by step](#method-step-by-step)
- [Output layout](#output-layout)
- [Figures](#figures)
- [Running on Slurm](#running-on-slurm)
- [Settings](#settings)
- [Design notes](#design-notes)

---

## Why the answers are not obvious

Three things make this harder than "BLAST the gene and take the best hit", and
most of Morpheus exists to handle them honestly.

**TOGA2's gene labels are not evidence.** In a family with many paralogs a
projection of one gene is routinely filed under a sibling's name. Morpheus never
takes the label at face value: a gene's locus is anchored on where its longest
human isoform actually projected, membership of that locus is decided by exon
overlap, and every candidate is screened against one protein per human gene
genome-wide.

**Sequence similarity cannot separate a tandem array.** Human IFITM1, IFITM2 and
IFITM3 lie within ~12 kb and share most of their coding sequence; every bat IFITM
scores highest against IFITM3 whichever locus it came from. Treating "best hit is
a different gene" as disqualifying deletes two of the three family members. So
inside a paralog family, position decides; only an out-of-family match is
disqualifying.

**The most similar sequence is not always in the expected place, and the thing in
the expected place is not always transcribed.** In *Phyllostomus discolor* the
OAS1 transcript carrying the C-terminal CaaX motif is a retrogene on a different
scaffold; the syntenic locus has no CaaX-bearing model at all. Restricted to the
locus, the animal reads as having lost the domain. Unrestricted, it plainly has
not — but a processed pseudogene may have no promoter, enhancer or TSS and never
be transcribed. Both readings are defensible, so Morpheus computes both and
reports where they disagree, rather than choosing on your behalf.

---

## Install

```bash
git clone https://github.com/seLoktionBias/Morpheus.git
cd Morpheus
bash install.sh
```

This creates a conda environment named `Morpheus`, verifies every external tool,
runs a built-in self test, and installs a `morpheus` launcher into the
environment. **Activating the environment is then all you need** — no PATH
editing, and no full paths in job scripts:

```bash
conda activate Morpheus     # or: mamba activate Morpheus
morpheus env                # confirm what was resolved
```

The launcher is a small script in the environment's `bin/` that hands off to this
checkout. It records the checkout's absolute path, so if you later move or delete
the checkout it says exactly that rather than failing with a bare "No such file
or directory". Re-run `install.sh` from the new location to repoint it, or run
several checkouts side by side with `--env=NAME` per checkout.

Backends, for when a solve is heavy or you have no conda at all:

```bash
bash install.sh --backend=mamba
bash install.sh --backend=conda
bash install.sh --backend=staged    # small solves, for memory-limited login nodes
bash install.sh --backend=current   # no environment; check the active interpreter
bash install.sh --check-only        # verify tools, create nothing
bash install.sh --no-launcher       # do not add `morpheus` to the environment
```

With `--backend=current` or `--no-launcher` there is no environment to install
into, so put the checkout on PATH yourself:

```bash
export PATH="/path/to/Morpheus/bin:$PATH"
```

On a managed cluster, load your site's conda module first — the installer never
calls `module load` for you:

```bash
ml miniforge && bash install.sh
```

`make help` lists the same operations as Makefile targets. `make test` runs the
smoke test (62 checks: the optimal matcher against brute force, the codon table,
ORF reporting, Newick pruning, every CLI subcommand, flag rejection, `paths.txt`
parsing, help/parser agreement, contiguous directory numbering, Slurm script consistency, all shell
syntax, and every figure script against synthetic input).

### Requirements

`python3` ≥ 3.9 (no third-party packages), `blastx`, `blastp`, `makeblastdb`,
`mafft`, and `R` with `ape`, `ggplot2`, `paletteer`. With `--skip_plot`, R is not
needed at all.

---

## Quick start

One gene, everything else taken from `paths.txt` in the current directory:

```bash
morpheus run --gene OAS1
```

A list of genes, on the cluster, submitted as one dependency chain:

```bash
morpheus run --gene_list gene_list.txt --path_file paths.txt --mode slurm
```

Nothing but the tables, and only the syntenic locus:

```bash
morpheus run --gene_list gene_list.txt --search region --skip_plot
```

### Every flag

| Flag | Meaning |
|---|---|
| `--gene NAME` | a single gene symbol; `--gene_list` is then not needed |
| `--gene_list FILE` | a file with one gene symbol per line; `--gene` is then not needed |
| `--mode local\|slurm` | run here, or submit the whole chain (default `local`) |
| `--search region\|similarity\|both` | how much to search (default `both`, see below) |
| `--tree FILE` | species tree in Newick, full path including the filename |
| `--bat_dir DIR` | TOGA2 annotation directory, one sub-directory per genome |
| `--ref_dir DIR` | human genome and Ensembl annotation directory |
| `--output DIR` | working/output directory (default: the current directory) |
| `--path_file FILE` | supply `--bat_dir`, `--ref_dir` and `--output` from a file instead |
| `--skip_plot` | write no figures, only tables — saves time and disk |
| `--per_gene yes\|no` | `yes` (default): one Slurm job per gene, or one gene at a time locally. `no`: one pass over the whole list |
| `--from STEP` | resume from a step |
| `--only STEP` | run one step (repeatable) |
| `--species NAME...` | restrict the search to these `Genus_species` |
| `--overwrite` | redo per-species searches already on disk |
| `--dry-run` | print the plan and stop; works in `--mode slurm` too |
| `--env_mode MODE` | `inherit` (default), `conda`, `mamba` or `venv` — see [Running on Slurm](#running-on-slurm) |
| `--env_name NAME` | environment for `--env_mode conda\|mamba` (default `Morpheus`) |
| `--venv_path DIR` | virtualenv for `--env_mode venv` |

Site options for `--mode slurm` (ignored, with a warning, in `local` mode):

| Flag | Meaning |
|---|---|
| `--time SPEC` | wall clock for every job: hours (`12`), `HH:MM:SS`, or `D-HH:MM:SS`. Omit to keep each job's own default |
| `--slurm_partition NAME` | queue to submit to, e.g. `compute` |
| `--slurm_account NAME` | account to charge |
| `--slurm_module NAME` | site module providing conda (default `miniforge`; `""` if conda is already on PATH) |
| `--array_throttle N` | at most N search tasks running at once |
| `--slurm_extra "..."` | passed verbatim to every `sbatch`, e.g. `"--qos=short"` |

Give either `--gene` or `--gene_list`, not both. Give either `--path_file` or the
`--bat_dir` / `--ref_dir` / `--output` trio; a flag always overrides the file, so
you can keep one `paths.txt` and redirect a single run with `--output`.

**What `--search` selects.** Scope (*where* to look) and ranking (*how* to choose
between what you find) are separate axes, and this one flag deliberately couples
them into the two combinations that answer a real question:

| `--search` | Scope | Ranking | The question it answers |
|---|---|---|---|
| `region` | the gene's own syntenic locus | exon structure counted | *what does this locus produce?* — the conservative orthology answer |
| `similarity` | the gene anywhere, paralogs still excluded | sequence identity | *what is the closest match this animal has?* — finds a functional copy elsewhere when the locus has diverged |
| `both` | both scopes | both rankings | all four combinations, side by side (default) |

`region` and `similarity` each run roughly half the work of `both`. `both` also
writes the scope comparison, which needs both scopes to exist; with `region` or
`similarity` that step and its figure are skipped automatically.

### Partial runs

```bash
morpheus run --dry-run           # print the plan and stop
morpheus run --from assign       # resume from a step
morpheus run --only copy-number  # a single step
morpheus env                     # what did it resolve? what is missing?
```

Every step reads only what earlier steps wrote, so re-running one never
invalidates the others. Expensive work is cached against a content hash, not a
timestamp, so re-running a step to change downstream logic does not repeat a
BLAST that has not actually changed.

Other entry points:

```bash
morpheus step bat-search --help  # the Python CLI directly
morpheus plot exon-models --gene MX1
morpheus version
```

`morpheus env` is the first thing to run when something will not start: it prints
the resolved tool paths and every input path with `ok` or `MISSING` beside it.

---

## Inputs

### `paths.txt`

Morpheus needs three directories. Put them in a `paths.txt` and you never type
them again. A copy ready to edit is in [`examples/paths.txt`](examples/paths.txt).

```
# Morpheus paths.txt - one key=value per line.
# Blank lines and #-comments are ignored; trailing comments are stripped;
# surrounding quotes are tolerated; a leading ~ is expanded; keys are
# case-insensitive.

human_genome_dir=/data/reference/human
bat_annotation_dir=/data/toga2/bat1k
primary_working_dir=/scratch/me/isg_project
```

That is the whole required file. The keys, in full:

| Key | Also accepted as | Required | What it is |
|---|---|---|---|
| `human_genome_dir` | `human_dir`, `ref_dir`, `reference_dir` | yes | human genome + Ensembl annotation |
| `bat_annotation_dir` | `bat_dir`, `toga_dir`, `toga2_dir`, `annotation_dir` | yes | TOGA2 output, one sub-directory per genome |
| `primary_working_dir` | `working_dir`, `output`, `output_dir`, `outdir` | yes | where `results/` is created |
| `tree` | `species_tree`, `tree_file` | no | species tree; defaults to `bat1k_tree.nwk` in the project directory |
| `gene_list` | `genes`, `gene_list_file` | no | defaults to `gene_list.txt` in the project directory |

A key that is not on this list is a **typo, and is reported as an error** rather
than ignored — a silently dropped `hooman_genome_dir` would otherwise send the
whole run at a default directory that does not exist, and you would find out an
hour later.

By default Morpheus reads `paths.txt` from the directory you run in. `--path_file
/elsewhere/paths.txt` reads it from there instead, and treats that directory as
the project, so `gene_list.txt` and the tree are looked for next to it.

Every key has a command-line equivalent that overrides it:

```bash
# no paths.txt at all
morpheus run --gene OAS1 \
  --ref_dir /data/reference/human \
  --bat_dir /data/toga2/bat1k \
  --output  /scratch/me/isg_project \
  --tree    /data/reference/bat1k_tree.nwk

# shared paths.txt, but this run writes somewhere else
morpheus run --gene_list genes.txt --path_file /shared/paths.txt --output ./run2
```

### The files themselves

| Input | Where | Notes |
|---|---|---|
| `gene_list.txt` | project root | one gene symbol per line |
| `bat1k_tree.nwk` | project root | species tree; tips are `Genus_species` |
| `Homo_sapiens.*.gtf.gz` | `human_genome_dir` | Ensembl GTF, any release |
| `Homo_sapiens.*.cds.all.fa` | `human_genome_dir` | Ensembl CDS FASTA |
| `longest_pc_tx.tsv` | `human_genome_dir` | `gene_id, transcript_id, cds_length, gene_symbol` |
| one directory per genome | `bat_annotation_dir` | TOGA2 output |

Genome directory names follow the Bat1K convention
`Genus_species__common_name__HLxxx__assembly`. The prefix is matched against the
tree tips longest-first, so trinomial tips such as
`Rhinolophus_perniger_lanosus` resolve correctly instead of being truncated to a
binomial. A genome with no matching tip is named in a warning and excluded
rather than silently dropped.

Each TOGA2 directory is read for `query_annotation.bed.gz`,
**`processed_pseudogenes.bed.gz`**, **`fragmented_annotation.bed.gz`**,
`nucleotide.fa.gz`, `query_genes.bed.gz`, `loss_summary.tsv.gz` and
`orthology_classification.tsv.gz`. All three BED files matter: the OAS1 retrogene
described above lives in `processed_pseudogenes.bed`, and reading only the main
annotation reports that animal as having lost the domain. Projections are
deduplicated across the three and tagged with an `annotation_source` column.

Nothing is ever written to either input directory, and compressed files are read
in place — nothing is unpacked.

---

## Two search scopes

`--search` chooses between them; `both` is the default. Two questions need two
search scopes, and they genuinely disagree.

| Scope | Searches | Answers |
|---|---|---|
| `region_restricted` | only the locus of the projected longest isoform | what does this gene's own syntenic locus produce? |
| `unrestricted` | the gene anywhere in the genome, paralogs still excluded | does the animal make this protein at all? |

Transcript status uses the region-restricted scope: a paralog elsewhere cannot
stand in for a transcript the locus itself has lost. The sequence sets use
the unrestricted scope: if the animal makes the protein from a different genomic
address, that is the sequence to analyse. The unrestricted scope is still
delimited *by gene* — outside the locus a candidate is kept only when the
whole-proteome screen says the gene it most resembles is the target itself.

`scope_comparison__<policy>.tsv` records, for every species × gene × transcript,
whether the two scopes chose the same model:

| Outcome | Meaning |
|---|---|
| `SAME_MODEL` | both scopes chose the same projection |
| `DIFFERENT_LOCUS` | the unrestricted scope found a better model elsewhere |
| `ONLY_UNRESTRICTED` | nothing acceptable in the locus; found outside it |
| `ONLY_REGION_RESTRICTED` | the locus model lost its slot once rivals were allowed |
| `NEITHER` | no acceptable model in either scope |

---

## Two ranking policies

The scope decides **where** a candidate may come from. The ranking decides **how**
it is scored once found. They are independent, and Morpheus builds both rankings
side by side, keeping them separate through every downstream step.

Every candidate is scored on the same five terms:

```
identity 0.45 · coverage 0.25 · structure 0.12 · positives 0.10 · length 0.08
```

The two policies differ in exactly one rule:

| Policy | Question | Treatment of retrocopies |
|---|---|---|
| `sequence_similarity` | which query sequence most resembles this transcript? | the exon-structure term is **dropped**; remaining weights renormalise |
| `structure_aware` | which query model has the architecture of the real gene? | the structure term **applies**, so an intronless copy pays for having no structure to match |

A retrocopy is a reverse-transcribed mRNA with no introns, so its exon-structure
term is necessarily zero. Scoring it there penalises a processed copy for being
processed, which is circular if the question is sequence resemblance — but it is
the right penalty if the question is which model is the functional gene, because
a processed pseudogene may never be transcribed at all.

> **The ranking is not positional.** No scoring term reads a coordinate, a
> flanking gene, or a distance from the home locus — `structure` is the longest
> common subsequence of *exon labels*, nothing more. The positional constraint on
> orthology lives on the other axis, in the `region_restricted` scope. Keeping
> the two apart is the whole point of having a 2×2 rather than one setting.
>
> This policy was called `synteny_aware` before v2.1.0. The name promised
> positional reasoning the ranking never did, and `structure_aware` reads
> correctly next to `unrestricted` where the old name read as a contradiction.

Measured on the 26-gene, 103-species reference run, the two rankings differ on
**6.6%** of assignments — and in most of those the *same* query models are
recovered and merely re-paired to different human isoforms, because the matching
is one-to-one. Genuinely different models are chosen in 42 of 116 disagreeing
gene × species groups. The differences concentrate in OAS1, IFITM3, IFITM1 and
IFITM2.

`DEFAULT_POLICY` (default `structure_aware`) drives the summary; both trees are
always built.

### Migrating output from an earlier version

Results produced before v2.1.0 carry the old directory name, the old file names,
and - easy to miss - the old name inside `manifest.tsv`, which records an
absolute path to every FASTA and tree it wrote. Renaming the directories alone
leaves that manifest pointing at paths that no longer exist.

```bash
morpheus migrate-names --results /path/to/results          # dry run
morpheus migrate-names --results /path/to/results --apply
```

It renames directories deepest-first, renames files, then rewrites the old name
inside text files only - binary files are left untouched. A rename that would
overwrite something is skipped and reported, never forced, and running it twice
is a no-op. Run it wherever the results live: migrating a local copy does not
migrate the one on the cluster.

---

## Method, step by step

### 1. Human reference — `cds-table`, `human-reference`

The Ensembl GTF is flattened once into a CDS exon table and cached. For each
gene the longest protein-coding isoform defines a **home locus**, and the gene's
transcript set is every protein-coding transcript in that region, on the same
strand, sharing at least one CDS exon with the longest isoform.

Transcripts with an identical CDS are collapsed to one representative —
everything downstream is sequence-based, so two transcripts with the same CDS are
one analysis. Both equivalence classes are recorded: `identical_model_transcripts`
(same CDS *and* same exon coordinates) and
`identical_cds_other_model_transcripts` (same CDS, different exon model, e.g.
alternative first exons encoding the same bases).

Exons get labels — `exon1..exonN` from the longest isoform, `human_novel_exonN`
for isoform-specific exons — and these are what get projected into each genome.
Flanking genes are recorded in transcription order, and a genome-wide
Ensembl-ID-to-symbol map is written so query gene neighbours read as symbols
rather than accessions.

### 2. Paralog families — `families`

Each gene's longest human protein is BLASTed against one protein per human gene.
Genes within 25% of the self-hit bitscore join the family. Derived, not
hard-coded, so it adapts to any gene list.

Membership is **symmetrised** before use. BLAST membership is directional — OAS1
reaches OAS3 above threshold while OAS3 does not reach OAS1 — and left as it
comes, "is this a sibling?" gets different answers depending on which gene is
asked, which made the screen discard OAS3's own locus in several species.

### 3. Region-restricted search — `bat-search`

Per genome, per gene, three candidate pools:

| Pool | What it is | Used for |
|---|---|---|
| `IN_REGION` | every projection sharing **at least one exon** with the projected longest isoform, whatever name TOGA2 gave it | transcripts + copy number |
| `OFF_REGION` | projections of this gene's human transcripts that landed elsewhere | transcripts + copy number |
| `FAMILY` | projections filed under a paralog's name, anywhere | copy number only |

Membership is decided by **exon overlap with the anchor**, not by falling inside
its span. Small genes frequently sit inside a large gene's intron and TOGA2
sometimes files them under the host's name; a span test sweeps them in as
transcripts of the host, an exon test cannot. Rejected projections are recorded
in `all_species_excluded_nested.tsv`.

Query exons inherit a human label only on an exact projected-coordinate match;
anything else becomes `bat_novel_exonN`, so a shifted splice site is never
quietly reported as `exon3`.

### 4. Paralog screen — `screen`

Every candidate is BLASTed against the one-protein-per-human-gene database:

| Verdict | Meaning |
|---|---|
| `CONSISTENT_WITH_TARGET_GENE` | best genome-wide match is the target |
| `CONSISTENT_BY_SYNTENY_WITHIN_FAMILY` | a sibling scores higher, but it sits in the target's own locus |
| `BELONGS_TO_ANOTHER_FAMILY_MEMBER` | a sibling scores higher and it is outside the locus |
| `LOOKS_LIKE_PARALOG` | a gene from *outside* the family scores better |
| `NO_SIMILARITY_TO_TARGET_GENE` / `NO_HIT` | no usable hit |

The BLAST is cached against a content hash of the query, so re-running to change
verdict logic does not repeat it.

### 5. Transcript assignment — `assign`

Two passes, in order.

**Pass 1 — projection identity.** TOGA2 names every projection
`<human_transcript>#<gene>#<chain>`, so a projection already records which human
transcript it was built from. That is not an estimate, it is the construction.
Each human transcript takes the projection built from it; a projection from a
collapsed transcript resolves to its representative.

**Pass 2 — similarity matching for the remainder.** Many human transcripts are
never projected at all. For those, the full pairwise matrix is solved as a
**maximum-weight bipartite matching** (Hungarian), so one query transcript serves
at most one human transcript and vice versa.

The matching is *provably independent of the order* rows and columns are
presented in — verified against brute force on 400 shuffled matrices in the smoke
test. A greedy "best hit per transcript" would both let two isoforms claim the
same model and let a mediocre early pick block a better global pairing.

Pair scores combine, with weights renormalised when a term is dropped:

```
0.45  identity over the longer sequence   (all non-overlapping HSPs, not the best one)
0.10  positives
0.25  alignment coverage
0.12  exon-structure agreement            (ordered LCS over exon labels)
0.08  CDS length similarity
```

Identity is computed **globally**. A single best HSP is a local view: two models
can have near-identical best HSPs while one carries far more indels and loses a
terminal domain, and a local figure ranks them the wrong way round.

Every row records `assignment_basis` (`projection_identity` /
`similarity_matching` / `none`), the highest-identity candidate regardless of
which slot it won (`sequence_best_candidate`), and a note whenever the assigned
model is not the sequence-best one.

### 6. Copy number — `copy-number`

Region-aware, label-independent, and non-overlapping:

1. take all three candidate pools;
2. keep only candidates the screen attaches to the target gene;
3. cluster into loci once per *paralog family*, joining only where projections
   actually **overlap**;
4. award each locus to exactly one gene — the gene whose own syntenic locus it
   is, else the strongest sequence evidence;
5. rejoin neighbouring loci awarded to the *same* gene, walking every gene's loci
   together in coordinate order so a merge never spans a sibling's copy;
6. record the loci a gene claimed but lost, as `shared_copies`.

Steps 3 and 5 are ordered deliberately. Merging on proximity *before* the award
collapses a tandem array into one locus — OAS1/OAS2/OAS3 sit a few kb apart, so a
10 kb tolerance swallowed all three into a single "copy" awarded to OAS3, leaving
the other two reporting zero copies while their transcripts were plainly
complete.

Step 6 exists because a sole award is sometimes a fiction: one bat IFITM locus is
the orthologous position of human IFITM2 *and* IFITM3. Counts are reported three
ways — `unambiguous_copies` (these do partition the genome), `shared_copies`, and
`total_copies`.

Each copy is classified `functional`, `coding_but_disrupted`,
`retro_or_processed`, or `lost_or_uncertain`.

### 7. Deliverables — `deliverables`

Per policy, two trees:

```
05_genes/<policy>/all_gene_files/       every transcript that recovered any model
05_genes/<policy>/gene_files_selected/  those in >= MIN_SPECIES_FRACTION of query species
```

Each transcript directory holds:

```
<GENE>__<TRANSCRIPT>.cds.fa        CDS multifasta, one sequence per species
<GENE>__<TRANSCRIPT>.tree.nwk      species tree pruned to exactly those species
<GENE>__<TRANSCRIPT>.members.tsv   provenance for every sequence
<GENE>__<TRANSCRIPT>.codon.aln.fa  codon alignment            (align step)
```

The alignment and its pruned tree are named alike and hold exactly the same taxa,
so a codon-model tool can be pointed at the pair directly:

```bash
cd results/05_genes/structure_aware/gene_files_selected/OAS1/OAS1__ENST00000202917
hyphy absrel --alignment OAS1__ENST00000202917.codon.aln.fa \
             --tree      OAS1__ENST00000202917.tree.nwk
```

FASTA headers are the bare tree tip label — `Homo_sapiens`, `Myotis_myotis` — so
the multifasta and the Newick agree character for character. The projection each
sequence came from is in `members.tsv`, not stuffed into the header.

Only the selected set goes on to alignment: a transcript recovered in three
species out of a hundred cannot support a codon model, and aligning it pads the
results with noise. Everything is still kept and can be inspected.

### 8. Alignment — `align`

CDS are translated, aligned as protein with MAFFT L-INS-i, and back-translated.
Aligning amino acids cannot open a gap that is not a multiple of three, which is
what codon models require.

TOGA2 models are frequently frameshifted or carry premature stops. Rather than
discarding them, each sequence is trimmed to whole codons and internal stops are
masked to `NNN`, so one disrupted codon does not truncate the protein and wreck
the alignment. Every edit is recorded in `alignment_preparation.tsv`.

---

## Output layout

With `--per_gene yes` (the default) each gene gets its own complete results tree,
and they are merged afterwards:

```
<output>/
├── cache/                       flattened Ensembl GTF, written once, shared
├── logs/
├── genes/
│   ├── MX1/results/             a full 01..06 tree, exactly as below
│   └── OAS1/results/
└── combined/                    every gene's tables concatenated
    ├── copy_number_matrix.tsv
    ├── transcript_status_<scope>__<policy>.tsv
    ├── transcript_assignments_<scope>__<policy>.tsv
    └── plots/                   the whole-list figures
```

Isolation is the point: nothing is shared but the read-only GTF cache, so
concurrent gene jobs cannot corrupt a common table and one gene failing leaves
every other gene's results intact and re-runnable on its own:

```bash
morpheus run --gene OAS1 --output <output>/genes/OAS1
```

The cost is real: each gene job reads every query genome's annotation, so *N*
genes read them *N* times. `--per_gene no` runs one pass over the whole list
instead, reading each genome once and writing a single tree:

```
results/
├── cache/                       flattened Ensembl CDS exon table (reusable)
├── logs/                        one log per step
├── 01_human_reference/          gene_context (with flanking genes), transcripts,
│                                exon labels, human_isoform_exons,
│                                ensembl_gene_id_to_symbol, CDS/protein FASTA,
│                                gene_families
├── 02_bat_search/               all_species_candidates, all_species_loci,
│                                all_species_candidate_exons,
│                                all_species_excluded_nested,
│                                all_species_candidate_cds.fa, paralog_screen,
│                                per_species/
├── 03_transcript_assignment/    transcript_assignments_<scope>__<policy>.tsv
│                                pairwise_similarity_<scope>__<policy>.tsv
│                                scope_comparison__<policy>.tsv
│                                candidates_excluded_<scope>__<policy>.tsv
│                                one_to_one_violations_<scope>__<policy>.tsv
├── 04_copy_number/              copy_number_matrix, copy_number_wide, gene_copies,
│                                contested_loci, overlapping_copy_loci
├── 05_genes/<policy>/           all_gene_files/, gene_files_selected/,
│                                manifest.tsv, transcript_status_<scope>.tsv
├── 06_plots/                    copy-number figure (policy-independent)
│   └── <policy>/                status and scope-comparison figures
├── SUMMARY__<policy>.md         start here
└── SUMMARY__<policy>.tsv
```

The six step directories are numbered `01`–`06` with no gaps. Each is created up
front, so a step that has not run yet leaves a visible empty directory carrying a
`.not_run_yet` marker rather than a hole in the numbering; the marker is removed
as soon as the step produces something. `--skip_plot` omits `06_plots/` entirely.

Scratch BLAST databases live under `02_bat_search/*_work/` and
`03_transcript_assignment/blast_work/` and can be deleted at any time.

---

## Figures

Produced by `morpheus run --only plots`:

- **`<policy>/transcript_status_region_restricted`** and
  **`<policy>/transcript_status_unrestricted`** — species tree beside a matrix of
  one column per human transcript, grouped by gene, coloured by CDS status
  (complete / partial / fragmented / pseudogenized / not_found). Comparing the
  pair shows where confining the search to the syntenic locus changes the answer.
  Columns recovered in fewer than `--min-species` species (default 50) are
  dropped; `--min-species 0` keeps everything.
- **`<policy>/transcript_scope_comparison`** — same layout, coloured by whether
  the two scopes agreed, with a black outline where the unrestricted search
  recovers a complete ORF the syntenic locus does not.
- **`gene_copy_number_phylogeny`** — one column per gene, one colour per copy
  number with no lumped top category, from
  `paletteer_d("colorBlindness::Blue2DarkOrange12Steps")`. Past 12 levels the
  ramp is rebuilt blue-to-dark-orange through that palette's saturated anchors,
  dropping its near-white middle, so a long scale stays legible.
  `--metric total_copies | copies_excluding_retro | functional_copies`.

**Optional, not run by the pipeline** — `plot_exon_models.R` draws UCSC-style
exon models for one gene. Exons keep true genomic width while introns are
compressed by one shared transform per gene, so equivalent exons line up
vertically. Every gene is drawn 5′→3′ left to right whatever strand it is on;
only the reference row carries the full `exon1..exonN` labels, the rest label
only their novel exons.

```bash
# human models for one gene
morpheus plot exon-models --gene MX1

# query models for one species and gene
morpheus plot exon-models \
  --plot-tsv results/02_bat_search/all_species_candidate_exons.tsv \
  --context-tsv results/02_bat_search/all_species_loci.tsv \
  --gene OAS1 --species Phyllostomus_discolor

# every gene at once
morpheus plot exon-models --all-genes
```

Tunables: `--intron-mode log|sqrt|linear|none`, `--intron-factor`,
`--min-exon-width`, `--label-all-rows`, `--no-flip-minus`, `--no-exon-labels`,
`--format pdf|png|both`.

---

## Running on Slurm

```bash
morpheus run --gene_list gene_list.txt --path_file paths.txt --mode slurm \
    --slurm_partition compute
```

That is the whole thing. It resolves the inputs, counts the genes, and submits a
dependent chain. **One job per gene**, by default:

```
morpheus_reference   flatten the GTF once, into the shared cache   32G  6h
morpheus_gene        array 1-N, one task per gene                  16G  12h
morpheus_collect     merge every gene, draw the figures            16G  4h
```

The collect job depends on the array with `afterany`, not `afterok`: one gene
failing must not stop the other twenty being merged and plotted. `merge-genes`
names any gene that produced nothing, so a partial result is never mistaken for
a complete one.

`--per_gene no` submits the other shape instead — one pass over the whole list,
parallelised over *genomes* rather than genes:

```
morpheus_reference     steps 1-3          32G   6h
morpheus_search        array 1-N          8G    2h    one task per genome
morpheus_post_search   screen -> align    32G   12h
morpheus_collect       figures, summary   16G   4h
```

It prints the job IDs and a ready-made `squeue -j` line. Logs land in
`<output>/logs/`. `--dry-run` shows what would be submitted without submitting
it, and works on a machine with no `sbatch` at all — useful for checking a
command on your laptop before running it on the cluster.

`--search` and `--skip_plot` are carried into every job, so one submitted chain
runs the same analysis end to end.

On a cluster that needs a named queue:

```bash
MORPHEUS=/data/home/me/software/Morpheus/bin/morpheus

"$MORPHEUS" run \
    --gene_list "${wdr}/gene_list.txt" \
    --path_file "${wdr}/paths.txt" \
    --mode slurm \
    --slurm_partition compute \
    --slurm_module miniforge \
    --array_throttle 20
```

Per-job `--mem` and `--time` live in the `slurm/*.sbatch` files, because the four
jobs are not alike — the reference build wants 32G for the GTF parse, an array
task wants 8G. `--slurm_extra "--mem=64G"` overrides all four at once if you need
it to.

Note the shape of that command: **one submission for the whole gene list**, not a
loop over genes. The reference build, the paralog screen and the per-genome
search are shared across every gene in the list, so submitting per gene would
redo the entire 103-genome search once per gene. Pass the list and let the array
parallelise over genomes.

### Telling a job where its tools are

This is where cluster runs actually fail, so Morpheus does not guess. `--env_mode`
follows the same three options Pensieve uses:

| `--env_mode` | What the job does | Use it when |
|---|---|---|
| `inherit` *(default)* | nothing — uses the PATH of the shell that submitted the job, carried in by `--export=ALL` | you activate the environment before submitting. Simplest, and what most people want |
| `conda` / `mamba` | runs the work through `<manager> run -n <name>` | the job needs to activate for itself. **`conda activate` does not work in a non-interactive batch shell** — `<manager> run` does, which is the whole reason this mode exists |
| `venv` | sources `--venv_path/bin/activate` | you installed with `--backend=venv` |

```bash
# the usual way: activate, then submit
ml miniforge && mamba activate Morpheus
morpheus run --gene_list gene_list.txt --path_file paths.txt --mode slurm

# or let each job resolve the environment itself
morpheus run --gene_list gene_list.txt --path_file paths.txt --mode slurm \
    --env_mode mamba --env_name Morpheus --slurm_module miniforge
```

Whichever you pick, every job checks up front that `python3`, `mafft`, `blastx`
and `Rscript` are actually reachable, and **fails immediately, naming what is
missing and what to pass instead**:

```
[Morpheus] ERROR: not on PATH in this job: mafft blastx Rscript
  env_mode  inherit
  env_name  Morpheus
  module    <none requested>

  Either activate the environment before submitting and keep the
  default --env_mode inherit, or submit with:
    --env_mode mamba --env_name Morpheus --slurm_module <your conda module>
```

That check exists because of a specific failure: an environment that
half-activates lets the job run on and die several steps later complaining about
a missing file, with the real cause nowhere in the log. In a dependent chain
that shows up as `DependencyNeverSatisfied` on every downstream job and no
explanation anywhere. **When a chain stalls like that, read the `.err` of the
*first* job — the others never ran.**

```bash
squeue -u $USER                  # which are held
cat logs/morpheus_reference_*.err
```

### Getting `morpheus` on PATH in a job

`install.sh` puts a `morpheus` launcher inside the environment, so activating it
is enough — in an interactive session and in a job script alike:

```bash
ml miniforge
mamba activate Morpheus
morpheus run --gene_list gene_list.txt --path_file paths.txt --mode slurm
```

The four submitted jobs do not depend on your PATH: `morpheus run --mode slurm`
passes `MORPHEUS_HOME` to each of them, and each re-activates the environment
through `slurm/_common.sh`. Only submitting the sbatch scripts *by hand* requires
you to export `MORPHEUS_HOME` yourself.

If you installed with `--no-launcher` or `--backend=current`, use the checkout
directly instead:

```bash
export PATH=/path/to/Morpheus/bin:$PATH    # or call /path/to/Morpheus/bin/morpheus
```

### Which module to load

Written for Apocrita, where conda comes from `ml miniforge`. For a site that
names it differently:

```bash
export MORPHEUS_MODULE=<your conda module>   # or "" if conda is already on PATH
```

### Submitting by hand

The scripts work standalone if you would rather control the chain yourself.
They read `paths.txt` from the submission directory as usual.

```bash
export MORPHEUS_HOME=/path/to/Morpheus
cd /path/to/project

ref=$(sbatch --parsable $MORPHEUS_HOME/slurm/00_reference.sbatch)
n=$(find "$(awk -F= '/bat_annotation_dir/{print $2}' paths.txt)" \
      -maxdepth 1 -mindepth 1 -type d ! -name '.*' | wc -l)
arr=$(sbatch --parsable --dependency=afterok:$ref --array=1-$n \
      $MORPHEUS_HOME/slurm/01_bat_search_array.sbatch)
post=$(sbatch --parsable --dependency=afterok:$arr \
       $MORPHEUS_HOME/slurm/02_after_search.sbatch)
sbatch --dependency=afterok:$post $MORPHEUS_HOME/slurm/03_collect.sbatch
```

The search is the only embarrassingly parallel step — each genome is read
independently — so the array turns *N* × 2 minutes into 2 minutes. Array tasks
pass `--no-merge` and a single `merge-search` job combines them afterwards;
merging inside each task would race and leave partial combined tables.

`02_after_search` is the long pole, not the array: it BLASTs every candidate
against the whole human proteome. That result is cached against a content hash of
the query, so a rerun that changes only downstream logic does not repeat it.

## Settings

Override by exporting before running, or edit `config/config.sh`.

| Variable | Default | Meaning |
|---|---|---|
| `THREADS` | all cores | BLAST/MAFFT threads |
| `MIN_SPECIES_FRACTION` | 0.5 | share of query species needed to enter `gene_files_selected` |
| `PLOT_MIN_SPECIES` | `MIN_SPECIES_FRACTION` of the species searched | species needed for a transcript column to be plotted. Derived, not fixed: a hard 50 produced an empty figure whenever fewer than 50 genomes were involved |
| `MIN_ASSIGNMENT_SCORE` | 0.30 | below this a human transcript stays unassigned |
| `MIN_ASSIGNMENT_PIDENT` | 40 | protein identity floor for any assignment |
| `POLICIES` | both | ranking policies to build (normally set by `--search`) |
| `DEFAULT_POLICY` | `structure_aware` | which ranking the summary reports on |
| `FAMILY_MIN_FRACTION` | 0.25 | paralog family membership threshold |
| `COPY_CLUSTER_SLOP` | 10000 | bp gap tolerated when rejoining one gene's fragments |
| `MORPHEUS_ENV` | `Morpheus` | conda environment name |
| `MORPHEUS_MODULE` | `miniforge` | site module to load on a cluster |

```bash
THREADS=16 MIN_SPECIES_FRACTION=0.75 morpheus run --gene_list genes.txt
```

---

## Design notes

A few choices that are not obvious from the code.

**No selection analysis.** Morpheus produces alignments and trees and stops.
Choosing a model, deciding which branches are foreground, and picking a multiple-
testing correction are the parts of a selection analysis that carry the actual
claim, and they belong to whoever is making it — not to a pipeline that ran
overnight with defaults nobody revisited. Each `gene_files_selected/` directory
holds an alignment and a pruned tree with identical taxa and identical labels, so
pointing any codon-model tool at it is one command.

**No third-party Python.** The Hungarian solver, the Newick parser, the BED12 and
GTF readers and every table parser are standard library. A failed conda solve
costs you the external tools, never the pipeline logic, and the smoke test runs
anywhere Python does.

**Caching on content, not timestamps.** Re-running the search rewrites its FASTA
with an identical body; a newer mtime is not a reason to redo hours of BLAST.

**Constraints are verified, not assumed.** `one_to_one_violations_*.tsv` and
`overlapping_copy_loci.tsv` should always be empty, and are written every run so
the guarantee can be checked rather than trusted.

**Nothing is dropped silently.** Candidates excluded by scope, projections nested
in another gene's intron, loci contested between family members, transcripts
below the species threshold — each has its own table.

**Absence of a feature is not evidence against membership** when the absence is a
property of the category. A retrocopy has no introns; an IFITM in a tandem array
scores highest against its sibling. Both were, at one point, silently deleting
real biology here, and both are now handled explicitly.
