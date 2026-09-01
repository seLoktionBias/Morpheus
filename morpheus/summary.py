"""Final step - one place to look first.

The analysis writes a lot of tables, and the point of this step is that nobody
should have to open them to find out what happened. It joins every stage into a
single row per gene and writes both a TSV and a short readable report.
"""
from __future__ import annotations

import statistics
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from .common import log, read_tsv, to_float, to_int, write_tsv

SUMMARY_FIELDS = [
    "gene", "chrom", "human_transcripts", "transcripts_with_alignment",
    "species_with_any_transcript", "assigned_transcript_slots",
    "assigned_region_restricted", "assigned_unrestricted",
    "transcripts_complete_only_off_locus", "scope_disagreements",
    "copies_median", "copies_min", "copies_max",
    "species_with_extra_copy", "species_with_no_copy", "retro_copies",
]


def _median(values: Sequence[float], default=0):
    return statistics.median(values) if values else default


def build(results_dir, outdir=None, policy: str = "structure_aware") -> Path:
    results = Path(results_dir)
    outdir = Path(outdir) if outdir else results

    reference = read_tsv(results / "01_human_reference" / "gene_context.tsv")
    assign_dir = results / "03_transcript_assignment"
    region = read_tsv(assign_dir /
                      f"transcript_assignments_region_restricted__{policy}.tsv")
    assignments = read_tsv(assign_dir /
                           f"transcript_assignments_unrestricted__{policy}.tsv")
    comparison = read_tsv(assign_dir / f"scope_comparison__{policy}.tsv")
    copies = read_tsv(results / "04_copy_number" / "copy_number_matrix.tsv")
    copy_loci = read_tsv(results / "04_copy_number" / "gene_copies.tsv")
    manifest = read_tsv(results / "05_genes" / policy / "manifest.tsv")

    genes = [r["gene"] for r in reference if r.get("status") == "OK"]
    chrom = {r["gene"]: r.get("chrom", "NA") for r in reference}

    n_human_tx: Dict[str, int] = defaultdict(int)
    for r in read_tsv(results / "01_human_reference" / "transcripts.tsv"):
        n_human_tx[r["gene"]] += 1

    assigned: Dict[str, int] = defaultdict(int)
    assigned_region: Dict[str, int] = defaultdict(int)
    species_any: Dict[str, set] = defaultdict(set)
    for r in assignments:
        if r.get("assignment_status") != "ASSIGNED":
            continue
        assigned[r["gene"]] += 1
        species_any[r["gene"]].add(r["species"])
    for r in region:
        if r.get("assignment_status") == "ASSIGNED":
            assigned_region[r["gene"]] += 1

    disagree: Dict[str, int] = defaultdict(int)
    rescued: Dict[str, int] = defaultdict(int)
    for r in comparison:
        if r.get("outcome") in ("DIFFERENT_LOCUS", "ONLY_UNRESTRICTED",
                                "ONLY_REGION_RESTRICTED"):
            disagree[r["gene"]] += 1
        if (r.get("unrestricted_cds_status") == "complete"
                and r.get("region_cds_status") != "complete"
                and r.get("outcome") in ("DIFFERENT_LOCUS", "ONLY_UNRESTRICTED")):
            rescued[r["gene"]] += 1

    aligned: Dict[str, int] = defaultdict(int)
    for r in manifest:
        if r.get("passes_species_threshold") == "1":
            aligned[r["gene"]] += 1

    copies_by_gene: Dict[str, List[int]] = defaultdict(list)
    for r in copies:
        copies_by_gene[r["gene"]].append(to_int(r.get("total_copies")))
    retro: Dict[str, int] = defaultdict(int)
    for r in copy_loci:
        if r.get("copy_class") == "retro_or_processed":
            retro[r["gene"]] += 1

    rows = []
    for gene in sorted(genes):
        cn = copies_by_gene.get(gene, [])
        rows.append({
            "gene": gene, "chrom": chrom.get(gene, "NA"),
            "human_transcripts": n_human_tx.get(gene, 0),
            "transcripts_with_alignment": aligned.get(gene, 0),
            "species_with_any_transcript": len(species_any.get(gene, ())),
            "assigned_transcript_slots": assigned.get(gene, 0),
            "assigned_region_restricted": assigned_region.get(gene, 0),
            "assigned_unrestricted": assigned.get(gene, 0),
            "transcripts_complete_only_off_locus": rescued.get(gene, 0),
            "scope_disagreements": disagree.get(gene, 0),
            "copies_median": _median(cn),
            "copies_min": min(cn) if cn else 0,
            "copies_max": max(cn) if cn else 0,
            "species_with_extra_copy": sum(1 for x in cn if x > 1),
            "species_with_no_copy": sum(1 for x in cn if x == 0),
            "retro_copies": retro.get(gene, 0),
        })

    log(f"summary built from the '{policy}' ranking")
    tsv = write_tsv(outdir / f"SUMMARY__{policy}.tsv", rows, SUMMARY_FIELDS)
    md = outdir / f"SUMMARY__{policy}.md"
    _write_report(rows, results, md, policy)
    log(f"summary written to {tsv.name} and {md.name}")
    return tsv


def _write_report(rows: List[dict], results: Path, path: Path,
                  policy: str = "structure_aware") -> None:
    def table(headers, keys, sort_key=None, reverse=False):
        data = sorted(rows, key=sort_key, reverse=reverse) if sort_key else rows
        out = ["| " + " | ".join(headers) + " |",
               "|" + "|".join("---" for _ in headers) + "|"]
        for r in data:
            out.append("| " + " | ".join(str(r[k]) for k in keys) + " |")
        return "\n".join(out)

    total_rescued = sum(r["transcripts_complete_only_off_locus"] for r in rows)
    total_disagree = sum(r["scope_disagreements"] for r in rows)

    lines = [
        f"# Analysis summary - {policy} ranking",
        "",
        f"{len(rows)} genes, ranked under the **{policy}** policy. The other "
        "ranking is reported in its own summary alongside this one; they differ "
        "only where a processed copy competes with a syntenic model.",
        "",
        f"Full tables are under `{results.name}/`; this file is the index.",
        "",
        "## Transcript recovery",
        "",
        table(["Gene", "Chr", "Human tx", "Tx aligned", "Species with data",
               "Assigned slots"],
              ["gene", "chrom", "human_transcripts", "transcripts_with_alignment",
               "species_with_any_transcript", "assigned_transcript_slots"]),
        "",
        "## Does the search scope change the answer?",
        "",
        "Transcript status uses the gene's own syntenic locus; the sequence "
        "sets use the gene wherever it occurs, paralogs still excluded.",
        "",
        f"- **{total_disagree}** transcript/species cells differ between the two "
        "scopes.",
        f"- **{total_rescued}** are complete ORFs recoverable *only* by leaving "
        "the gene's own locus - the animal makes the protein, from a different "
        "genomic address.",
        "",
        table(["Gene", "Assigned (locus)", "Assigned (anywhere)",
               "Complete only off-locus", "Scope disagreements"],
              ["gene", "assigned_region_restricted", "assigned_unrestricted",
               "transcripts_complete_only_off_locus", "scope_disagreements"],
              sort_key=lambda r: -r["scope_disagreements"]),
        "",
        "## Copy number",
        "",
        table(["Gene", "Median", "Min", "Max", "Species >1 copy",
               "Species 0 copies", "Retro copies"],
              ["gene", "copies_median", "copies_min", "copies_max",
               "species_with_extra_copy", "species_with_no_copy", "retro_copies"],
              sort_key=lambda r: -r["copies_max"]),
        "",
        "## Where things are",
        "",
        "| What | Path |",
        "|---|---|",
        f"| Per-transcript FASTA, tree, alignment | `{results.name}/05_genes/{policy}/all_gene_files/<GENE>/<GENE>__<TX>/` |",
        f"| Subset analysed (>= min species) | `{results.name}/05_genes/{policy}/gene_files_selected/` |",
        f"| Transcript assignments | `{results.name}/03_transcript_assignment/transcript_assignments_<scope>__<policy>.tsv` |",
        f"| Sequences, sequence-similarity ranking | `{results.name}/05_genes/sequence_similarity/` |",
        f"| Sequences, structure-aware ranking | `{results.name}/05_genes/structure_aware/` |",
        f"| Region-restricted vs unrestricted | `{results.name}/03_transcript_assignment/scope_comparison.tsv` |",
        f"| Every pairwise similarity considered | `{results.name}/03_transcript_assignment/pairwise_similarity_<scope>.tsv` |",
        f"| Copy number per species | `{results.name}/04_copy_number/copy_number_matrix.tsv` |",
        f"| Each copy locus, with evidence | `{results.name}/04_copy_number/gene_copies.tsv` |",
        f"| Paralog screen verdicts | `{results.name}/02_bat_search/paralog_screen.tsv` |",
        f"| Figures | `{results.name}/06_plots/` |",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")
