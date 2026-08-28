# Morpheus

Comparative transcript recovery, gene copy number, and selection analysis across
TOGA2-annotated genomes.

Give it a list of human genes, a species tree, and a directory of TOGA2 results —
one per genome — and it recovers each gene's transcripts in every genome, decides
which query transcript corresponds to which human transcript, counts gene copies,
assembles ready-to-analyse alignments with matching pruned trees, and runs
selection tests that need no branch labelling.

Everything runs locally or on SLURM. The pipeline logic uses only the Python
standard library; the external tools are BLAST, MAFFT, HyPhy and R.

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
- [Running on SLURM](#running-on-slurm)
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
and runs a built-in self test. Then:

```bash
conda activate Morpheus
export PATH="$PWD/bin:$PATH"
morpheus env          # confirm what was resolved
```

Backends, for when a solve is heavy or you have no conda at all:

```bash
bash install.sh --backend=mamba
bash install.sh --backend=conda
bash install.sh --backend=staged    # small solves, for memory-limited login nodes
bash install.sh --backend=current   # no environment; check the active interpreter
bash install.sh --check-only        # verify tools, create nothing
```

On a managed cluster, load your site's conda module first — the installer never
calls `module load` for you:

```bash
ml miniforge && bash install.sh
```

`make help` lists the same operations as Makefile targets. `make test` runs the
smoke test (35 checks: the optimal matcher against brute force, the codon table,
ORF reporting, Newick pruning, every CLI subcommand, all shell syntax, and every
figure script against synthetic input).

### Requirements

`python3` ≥ 3.9 (no third-party packages), `blastx`, `blastp`, `makeblastdb`,
`mafft`, `hyphy`, and `R` with `ape`, `ggplot2`, `paletteer`.

---

## Quick start

Work from your project directory — the one holding `paths.txt`, `gene_list.txt`
and the species tree. Results are written to `results/` there.

```bash
morpheus run                     # every step, in order
morpheus run --dry-run           # print the plan and stop
morpheus run --from assign       # resume from a step
morpheus run --only copy-number  # a single step
morpheus run --species Myotis_myotis Desmodus_rotundus   # restrict the search
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

---

## Inputs

`paths.txt` in the project directory supplies the three locations:

```
human_genome_dir=/path/to/human_genome
bat_annotation_dir=/path/to/toga2_results
primary_working_dir=/path/to/project
```

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

Two questions need two search scopes, and they genuinely disagree.

| Scope | Searches | Answers |
|---|---|---|
| `region_restricted` | only the locus of the projected longest isoform | what does this gene's own syntenic locus produce? |
| `unrestricted` | the gene anywhere in the genome, paralogs still excluded | does the animal make this protein at all? |

Transcript status uses the region-restricted scope: a paralog elsewhere cannot
stand in for a transcript the locus itself has lost. Alignment and selection use
the unrestricted scope: if the animal makes the protein from a different genomic
address, that is the sequence selection acts on. The unrestricted scope is still
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

Within a scope, there are two defensible ways to rank a candidate. Morpheus
builds **both**, side by side, and keeps them separate through every downstream
step.

| Policy | Question | Treatment of retrocopies |
|---|---|---|
| `sequence_similarity` | which query sequence most resembles this transcript? | the exon-structure term is **dropped**; remaining weights renormalise |
| `synteny_aware` | which query model is the *ortholog*? | the structure term **applies**, so an intronless copy pays for having no structure to match |

A retrocopy is a reverse-transcribed mRNA with no introns, so its exon-structure
term is necessarily zero. Scoring it there penalises a processed copy for being
processed, which is circular if the question is sequence resemblance — but it is
exactly the right penalty if the question is orthology, because orthology is
positional and a processed pseudogene may never be transcribed.

They differ on **7.4%** of assignments in the reference dataset, concentrated
precisely where you would expect: OAS1, IFITM3, IFITM1, IFITM2.

`DEFAULT_POLICY` (default `synteny_aware`) drives the summary; both trees are
always built.

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
hyphy/                             HyPhy JSON and logs        (hyphy step)
```

FASTA headers are the bare tree tip label — `Homo_sapiens`, `Myotis_myotis` — so
the multifasta and the Newick agree character for character. The projection each
sequence came from is in `members.tsv`, not stuffed into the header.

Only the selected set goes on to alignment and selection: a transcript recovered
in three species out of a hundred cannot support a codon model, and aligning it
pads the results with noise. Everything is still kept and can be inspected.

### 8. Alignment — `align`

CDS are translated, aligned as protein with MAFFT L-INS-i, and back-translated.
Aligning amino acids cannot open a gap that is not a multiple of three, which is
what codon models require.

TOGA2 models are frequently frameshifted or carry premature stops. Rather than
discarding them, each sequence is trimmed to whole codons and internal stops are
masked to `NNN`, so one disrupted codon does not truncate the protein and wreck
the alignment. Every edit is recorded in `alignment_preparation.tsv`.

### 9. Selection — `hyphy`

aBSREL, BUSTED, MEME and FEL, all on the unlabelled tree — nothing depends on
choosing foreground lineages in advance. The tree is re-pruned to exactly the
taxa present in each alignment before running.

HyPhy scales poorly inside one analysis, so parallelism is *across* analyses:
`--jobs` runs several at once with `--hyphy-threads` cores each. `--longest-only`
analyses one representative transcript per gene. Resume only reuses a JSON that
actually parses, so an interrupted run cannot leave a truncated file that later
passes for finished.

The JSON is flattened into three tidy tables: `hyphy_gene_level.tsv`,
`hyphy_absrel_branches.tsv`, `hyphy_selected_sites.tsv`.

---

## Output layout

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
├── 06_selection/<policy>/       hyphy_gene_level, hyphy_absrel_branches,
│                                hyphy_selected_sites
├── 07_plots/                    copy-number figure (policy-independent)
│   └── <policy>/                status and scope-comparison figures
├── SUMMARY__<policy>.md         start here
└── SUMMARY__<policy>.tsv
```

Every step directory is created up front, so a step that has not run yet leaves a
visible empty directory with a `.not_run_yet` marker rather than a hole in the
numbering.

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

## Running on SLURM

Written for Apocrita (`ml miniforge`); set `MORPHEUS_MODULE` for a site that
names its conda module differently, or to the empty string if conda is already
on `PATH`.

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

m=$(awk -F'\t' 'NR>1 && $16=="1"' results/05_genes/synteny_aware/manifest.tsv | wc -l)
hy=$(sbatch --parsable --dependency=afterok:$post --array=1-$m%20 \
     $MORPHEUS_HOME/slurm/03_hyphy_array.sbatch)

sbatch --dependency=afterany:$hy $MORPHEUS_HOME/slurm/04_collect.sbatch
```

The search is the only embarrassingly parallel step — each genome is read
independently — so the array turns *N* × 2 minutes into 2 minutes. Array tasks
pass `--no-merge` and a single `merge-search` job combines them afterwards;
merging inside each task would race and leave partial combined tables.

`04_collect` uses `afterany`, not `afterok`: one analysis timing out should not
stop the rest from being collected. Failures are recorded in
`hyphy_gene_level.tsv`.

---

## Settings

Override by exporting before running, or edit `config/config.sh`.

| Variable | Default | Meaning |
|---|---|---|
| `THREADS` | all cores | BLAST/MAFFT threads |
| `MIN_SPECIES_FRACTION` | 0.5 | share of query species needed to enter `gene_files_selected` |
| `PLOT_MIN_SPECIES` | 50 | species needed for a transcript column to be plotted |
| `MIN_ASSIGNMENT_SCORE` | 0.30 | below this a human transcript stays unassigned |
| `MIN_ASSIGNMENT_PIDENT` | 40 | protein identity floor for any assignment |
| `POLICIES` | both | ranking policies to build |
| `DEFAULT_POLICY` | `synteny_aware` | which ranking the summary reports on |
| `FAMILY_MIN_FRACTION` | 0.25 | paralog family membership threshold |
| `COPY_CLUSTER_SLOP` | 10000 | bp gap tolerated when rejoining one gene's fragments |
| `HYPHY_METHODS` | `absrel busted meme fel` | which analyses to run |
| `HYPHY_JOBS` / `HYPHY_THREADS` | cores/2, 2 | concurrent analyses × cores each |
| `HYPHY_TIMEOUT` | none | seconds per analysis |
| `HYPHY_LONGEST_ONLY` | 0 | 1 = one transcript per gene |
| `MORPHEUS_ENV` | `Morpheus` | conda environment name |
| `MORPHEUS_MODULE` | `miniforge` | site module to load on a cluster |

```bash
THREADS=4 HYPHY_METHODS="absrel busted" morpheus run --from hyphy
```

---

## Design notes

A few choices that are not obvious from the code.

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
