"""Command-line interface: `morpheus <step> [options]`.

Every step is independent and reads only files written by earlier steps, so any
part of the analysis can be re-run without repeating the rest.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from . import (align, bat_loci, compare, copy_number, deliverables, family,
               human_ref, merge_genes, pairwise, screen, summary)
from .common import log, read_lines, species_from_dirname


def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--threads", type=int, default=int(os.environ.get("BATCOMP_THREADS", "4")))


def cmd_cds_table(a) -> None:
    human_ref.build_cds_exon_table(a.gtf, a.out)


def cmd_human_reference(a) -> None:
    human_ref.build(a.genes, a.cds_exons, a.longest_pc_tx, a.cds, a.outdir)


def cmd_families(a) -> None:
    family.build(a.reference, a.outdir, min_bitscore_fraction=a.min_fraction,
                 threads=a.threads)


def cmd_bat_search(a) -> None:
    from . import newick
    fam = family.load(a.families) if a.families and Path(a.families).exists() else {}
    tips = newick.read_file(a.tree).leaf_names() if a.tree else None

    run_dirs = sorted(d for d in Path(a.annotations).iterdir()
                      if d.is_dir() and not d.name.startswith("."))
    resolved = [(d, species_from_dirname(d.name, tips)) for d in run_dirs]

    if tips:
        unmatched = sorted({sp for _, sp in resolved if sp not in set(tips)})
        if unmatched:
            log(f"WARNING: {len(unmatched)} genome(s) have no matching tree tip and "
                f"will be excluded downstream: {', '.join(unmatched)}")

    if a.species:
        wanted = set(a.species)
        resolved = [(d, sp) for d, sp in resolved if sp in wanted]
    if not resolved:
        log("no species directories matched")
        return

    outdir = Path(a.outdir) / "per_species"
    for i, (d, sp) in enumerate(resolved, 1):
        done = outdir / f"{sp}.candidates.tsv"
        if done.exists() and not a.overwrite:
            log(f"[{i}/{len(resolved)}] {sp}: already done (use --overwrite to redo)")
            continue
        log(f"[{i}/{len(resolved)}] {sp}")
        bat_loci.run_species(d, a.reference, outdir, locus_slop=a.locus_slop,
                             family_map=fam, tree_tips=tips)
    # In a SLURM array every task holds only its own genome, so merging there
    # would race and each task would overwrite the combined tables with a
    # partial set. The array passes --no-merge and one `merge-search` job
    # collects the lot afterwards.
    if a.no_merge:
        log("per-species output written; run `morpheus merge-search` to combine")
        return
    bat_loci.merge(outdir, a.outdir)


def cmd_merge_search(a) -> None:
    bat_loci.merge(Path(a.outdir) / "per_species", a.outdir)


def cmd_screen(a) -> None:
    screen.run(a.candidate_fasta, a.candidates, a.reference, a.outdir,
               threads=a.threads, families_tsv=a.families)


def cmd_assign(a) -> None:
    pairwise.run_pairwise_blast(a.candidates, a.candidate_fasta, a.reference,
                                a.outdir, threads=a.threads)
    for scope in a.scopes:
        for policy in a.policies:
            pairwise.assign(a.candidates, Path(a.outdir) / "pairwise_blastx.tsv",
                            a.reference, a.outdir, screen_tsv=a.screen,
                            min_score=a.min_score, min_pident=a.min_pident,
                            scope=scope, identity_first=not a.no_identity_first,
                            policy=policy)


def cmd_compare(a) -> None:
    compare.build(a.region_restricted, a.unrestricted, a.outdir,
                  sequences_fasta=a.candidate_fasta, policy=a.policy)


def cmd_copy_number(a) -> None:
    copy_number.run(a.candidates, a.screen, a.loci, a.candidate_fasta, a.outdir,
                    cluster_slop=a.cluster_slop, families_tsv=a.families)


def cmd_deliverables(a) -> None:
    outdir = Path(a.outdir)
    if a.policy:
        outdir = outdir / a.policy
    deliverables.build(a.assignments, a.candidate_fasta, a.reference, a.tree,
                       outdir, min_species_fraction=a.min_species_fraction,
                       scope=a.scope, write_directories=not a.status_only,
                       total_species=a.total_species)


def cmd_align(a) -> None:
    align.align_all(a.manifest, outdir=a.outdir, threads=a.threads,
                    only_passing=not a.include_failing)


def cmd_merge_genes(a) -> None:
    merge_genes.merge(a.genes_root, a.outdir, policies=a.policies,
                      scopes=a.scopes)


def cmd_summary(a) -> None:
    summary.build(a.results, a.outdir, policy=a.policy)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="morpheus",
        description="Region-aware TOGA2 transcript recovery, one-to-one human/bat "
                    "transcript assignment, and gene copy number.")
    sub = p.add_subparsers(dest="step", required=True)

    s = sub.add_parser("cds-table", help="flatten an Ensembl GTF into a CDS exon table (once)")
    s.add_argument("--gtf", required=True)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_cds_table)

    s = sub.add_parser("human-reference", help="build the human reference for the gene list")
    s.add_argument("--genes", required=True)
    s.add_argument("--cds-exons", required=True)
    s.add_argument("--longest-pc-tx", required=True)
    s.add_argument("--cds", required=True, help="Ensembl CDS FASTA")
    s.add_argument("--outdir", required=True)
    s.set_defaults(func=cmd_human_reference)

    s = sub.add_parser("families", help="define each gene's paralog family from the human proteome")
    s.add_argument("--reference", required=True)
    s.add_argument("--outdir", required=True)
    s.add_argument("--min-fraction", type=float, default=0.25,
                   help="minimum bitscore as a fraction of the self-hit (default 0.25)")
    _add_common(s)
    s.set_defaults(func=cmd_families)

    s = sub.add_parser("bat-search", help="region-restricted transcript search in every bat genome")
    s.add_argument("--annotations", required=True, help="directory of TOGA2 run directories")
    s.add_argument("--reference", required=True)
    s.add_argument("--families", default=None)
    s.add_argument("--tree", default=None,
                   help="species tree; used to resolve directory names to tip labels")
    s.add_argument("--outdir", required=True)
    s.add_argument("--species", nargs="*", default=None, help="restrict to these tip labels")
    s.add_argument("--locus-slop", type=int, default=0)
    s.add_argument("--overwrite", action="store_true")
    s.add_argument("--no-merge", action="store_true",
                   help="write per-species output only; for SLURM array tasks")
    _add_common(s)
    s.set_defaults(func=cmd_bat_search)

    s = sub.add_parser("merge-search",
                       help="combine per-species search output into the "
                            "all_species_* tables")
    s.add_argument("--outdir", required=True)
    s.set_defaults(func=cmd_merge_search)

    s = sub.add_parser("screen", help="screen every candidate against the whole human proteome")
    s.add_argument("--candidates", required=True)
    s.add_argument("--candidate-fasta", required=True)
    s.add_argument("--reference", required=True)
    s.add_argument("--families", default=None,
                   help="gene_families.tsv (default: <reference>/gene_families.tsv)")
    s.add_argument("--outdir", required=True)
    _add_common(s)
    s.set_defaults(func=cmd_screen)

    s = sub.add_parser("assign", help="pairwise BLAST and one-to-one transcript assignment")
    s.add_argument("--candidates", required=True)
    s.add_argument("--candidate-fasta", required=True)
    s.add_argument("--reference", required=True)
    s.add_argument("--screen", default=None)
    s.add_argument("--outdir", required=True)
    s.add_argument("--min-score", type=float, default=0.30)
    s.add_argument("--min-pident", type=float, default=40.0)
    s.add_argument("--no-identity-first", action="store_true",
                   help="skip the projection-identity pass and assign purely "
                        "by similarity (for comparison)")
    s.add_argument("--scopes", nargs="+", default=list(pairwise.SCOPES),
                   choices=list(pairwise.SCOPES),
                   help="search scopes to assign under (default: both)")
    s.add_argument("--policies", nargs="+", default=list(pairwise.POLICIES),
                   choices=list(pairwise.POLICIES),
                   help="ranking policies to assign under (default: both)")
    _add_common(s)
    s.set_defaults(func=cmd_assign)

    s = sub.add_parser("compare", help="do the two search scopes pick the same transcript?")
    s.add_argument("--region-restricted", required=True)
    s.add_argument("--unrestricted", required=True)
    s.add_argument("--candidate-fasta", required=True)
    s.add_argument("--outdir", required=True)
    s.add_argument("--policy", default="structure_aware",
                   help="label for this comparison's output files")
    s.set_defaults(func=cmd_compare)

    s = sub.add_parser("copy-number", help="region-aware gene copy number")
    s.add_argument("--candidates", required=True)
    s.add_argument("--candidate-fasta", required=True)
    s.add_argument("--screen", required=True)
    s.add_argument("--loci", required=True)
    s.add_argument("--outdir", required=True)
    s.add_argument("--cluster-slop", type=int, default=10000)
    s.add_argument("--families", default=None,
                   help="gene_families.tsv, so loci are built family-wide")
    s.set_defaults(func=cmd_copy_number)

    s = sub.add_parser("deliverables", help="per-gene/per-transcript FASTA + pruned tree")
    s.add_argument("--assignments", required=True)
    s.add_argument("--candidate-fasta", required=True)
    s.add_argument("--reference", required=True)
    s.add_argument("--tree", required=True)
    s.add_argument("--outdir", required=True)
    s.add_argument("--min-species-fraction", type=float, default=0.5,
                   help="fraction of query species a transcript must be "
                        "recovered in to join gene_files_selected (default 0.5)")
    s.add_argument("--total-species", type=int, default=None,
                   help="number of query genomes searched; inferred if omitted")
    s.add_argument("--scope", default="unrestricted",
                   help="label for the status table written by this run")
    s.add_argument("--policy", default=None,
                   help="write into this policy subdirectory "
                        "(sequence_similarity | structure_aware)")
    s.add_argument("--status-only", action="store_true",
                   help="write the status table but not the sequence directories")
    s.set_defaults(func=cmd_deliverables)

    s = sub.add_parser("align", help="codon-aware alignment of every transcript FASTA")
    s.add_argument("--manifest", required=True)
    s.add_argument("--outdir", default=None)
    s.add_argument("--include-failing", action="store_true")
    _add_common(s)
    s.set_defaults(func=cmd_align)

    s = sub.add_parser("merge-genes",
                       help="concatenate per-gene runs into whole-list tables")
    s.add_argument("--genes-root", required=True,
                   help="directory holding one sub-directory per gene")
    s.add_argument("--outdir", required=True)
    s.add_argument("--policies", nargs="+", default=list(pairwise.POLICIES))
    s.add_argument("--scopes", nargs="+", default=list(pairwise.SCOPES))
    s.set_defaults(func=cmd_merge_genes)

    s = sub.add_parser("summary", help="join every stage into SUMMARY.tsv / SUMMARY.md")
    s.add_argument("--results", required=True, help="the results/ directory")
    s.add_argument("--outdir", default=None)
    s.add_argument("--policy", default="structure_aware",
                   help="which ranking the summary reports on")
    s.set_defaults(func=cmd_summary)

    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
