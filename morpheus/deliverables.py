"""Step 5 - assemble the per-gene, per-transcript deliverables.

Two trees of per-transcript directories:

    all_gene_files/       every human transcript that recovered any model
    gene_files_selected/  only those present in at least a set fraction of the
                          query species (50% by default)

Alignment and selection run on the selected set alone. A transcript recovered in
three species out of a hundred cannot support a codon model or a branch test,
and aligning it wastes time and pads the results with noise; but it is still
evidence, so it is kept and can be inspected.

Each directory holds:

    <GENE>__<TRANSCRIPT>/
        <GENE>__<TRANSCRIPT>.cds.fa       CDS multifasta, one sequence per species
        <GENE>__<TRANSCRIPT>.tree.nwk     species tree pruned to exactly those species
        <GENE>__<TRANSCRIPT>.members.tsv  provenance for every sequence

FASTA headers are the bare tree tip label (Homo_sapiens, Myotis_myotis, ...), so
the multifasta and the Newick file agree exactly, character for character. The
TOGA2 projection each sequence came from is recorded in members.tsv rather than
being stuffed into the header.

This runs once per search scope. The unrestricted scope writes the sequence
directories, because alignment and selection should use the model the animal
actually has wherever it sits. Both scopes write a transcript-status table, and
they are plotted side by side: the region-restricted one says what the gene's own
syntenic locus produces, the unrestricted one says what the animal produces.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .common import (die, fasta_iter, log, orf_report, read_tsv, sanitize,
                     to_float, to_int, write_fasta, write_tsv)
from . import newick

MANIFEST_FIELDS = ["gene", "human_transcript", "directory", "n_sequences",
                   "n_query_species", "total_query_species",
                   "fraction_of_species", "min_fraction_required",
                   "human_present", "tree_tips",
                   "n_complete", "n_partial", "n_fragmented", "n_pseudogenized",
                   "passes_species_threshold", "fasta", "newick",
                   "selected_directory", "selected_fasta", "selected_newick"]

MEMBER_FIELDS = ["header", "species", "gene", "human_transcript", "source",
                 "candidate_id", "projection", "toga_transcript",
                 "toga_gene_label", "toga_status", "orthology_class", "pool",
                 "chrom", "start", "end", "strand",
                 "upstream_gene", "downstream_gene",
                 "distance_to_home_locus_bp", "screen_verdict",
                 "pident", "alignment_coverage", "similarity_score",
                 "sequence_length", "cds_status"]

STATUS_FIELDS = ["gene", "human_transcript", "species", "status", "toga_status",
                 "pool", "sequence_length", "similarity_score"]


# A stop this near the end of the supplied CDS is the gene's real stop codon,
# not a premature one: TOGA2 projections routinely carry a few trailing bases.
TERMINAL_STOP_FRACTION = 0.95


def cds_status(seq: str) -> str:
    """Classify a CDS for the transcript-status plot.

    complete       ATG start and an in-frame stop at the end of the ORF
    fragmented     no in-frame stop at all, or no ATG start
    partial        genuine premature stop, but past half the codons
    pseudogenized  premature stop in the first half
    """
    if not seq:
        return "not_found"
    r = orf_report(seq)
    if not r["n_codons"]:
        return "not_found"

    if r["stop_codon_index"] is None:
        return "fragmented"
    if r["orf_fraction"] >= TERMINAL_STOP_FRACTION:
        return "complete" if r["has_start"] else "fragmented"
    if not r["has_start"]:
        return "fragmented"
    return "partial" if r["orf_fraction"] > 0.5 else "pseudogenized"


def build(assignments_tsv, candidate_fasta, reference_dir, tree_path, outdir,
          min_species_fraction: float = 0.5, human_label: str = "Homo_sapiens",
          scope: str = "unrestricted", write_directories: bool = True,
          total_species: Optional[int] = None) -> Path:
    import shutil
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    all_dir = outdir / "all_gene_files"
    sel_dir = outdir / "gene_files_selected"

    assignments = [r for r in read_tsv(assignments_tsv)
                   if r.get("assignment_status") == "ASSIGNED"]
    bat_seqs = {h.split()[0]: s for h, s in fasta_iter(candidate_fasta)}
    human_seqs = {h.split()[0]: s for h, s in
                  fasta_iter(Path(reference_dir) / "human_transcript_cds.fa")}
    ref_tx = {r["transcript_id"]: r
              for r in read_tsv(Path(reference_dir) / "transcripts.tsv")}

    tree_tips = set(newick.read_file(tree_path).leaf_names())

    # "How many species were searched", not "how many happened to have a hit" -
    # a threshold measured against the surviving species would drift upwards as
    # recovery got worse, which is exactly backwards.
    n_total = total_species or len({r["species"] for r in read_tsv(assignments_tsv)})
    if n_total <= 0:
        die("cannot determine the total number of query species")

    by_transcript: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    for r in assignments:
        by_transcript[(r["gene"], r["human_transcript"])].append(r)

    manifest: List[dict] = []
    status_rows: List[dict] = []
    skipped_species = set()

    for (gene, tx), rows in sorted(by_transcript.items()):
        # sanitize() collapses runs of underscores, so build the "__"
        # separator after sanitizing each half.
        label = f"{sanitize(gene)}__{sanitize(tx)}"
        gene_dir = all_dir / gene / label
        records: List[Tuple[str, str]] = []
        members: List[dict] = []
        counts: Dict[str, int] = defaultdict(int)

        human_seq = human_seqs.get(tx, "")
        if human_seq:
            st = cds_status(human_seq)
            records.append((human_label, human_seq))
            counts[st] += 1
            members.append({
                "header": human_label, "species": human_label, "gene": gene,
                "human_transcript": tx, "source": "human_reference",
                "candidate_id": "NA", "projection": "NA",
                "toga_transcript": tx, "toga_gene_label": gene,
                "toga_status": "NA", "orthology_class": "NA", "pool": "NA",
                "chrom": ref_tx.get(tx, {}).get("chrom", "NA"),
                "start": ref_tx.get(tx, {}).get("start", "NA"),
                "end": ref_tx.get(tx, {}).get("end", "NA"),
                "strand": ref_tx.get(tx, {}).get("strand", "NA"),
                "upstream_gene": "NA", "downstream_gene": "NA",
                "distance_to_home_locus_bp": 0, "screen_verdict": "NA",
                "pident": 100.0, "alignment_coverage": 1.0,
                "similarity_score": "NA",
                "sequence_length": len(human_seq), "cds_status": st,
            })
            status_rows.append({"gene": gene, "human_transcript": tx,
                                "species": human_label, "status": st,
                                "toga_status": "NA", "pool": "REFERENCE",
                                "sequence_length": len(human_seq),
                                "similarity_score": "NA"})

        # one sequence per species; the assignment step already guarantees that
        seen_species: Dict[str, dict] = {}
        for r in sorted(rows, key=lambda x: (x["species"],
                                             -to_float(x.get("similarity_score")))):
            sp = r["species"]
            if sp not in tree_tips:
                skipped_species.add(sp)
                continue
            if sp in seen_species:
                continue
            seq = bat_seqs.get(r["candidate_id"], "")
            if not seq:
                continue
            seen_species[sp] = r
            st = cds_status(seq)
            counts[st] += 1
            records.append((sp, seq))
            members.append({
                "header": sp, "species": sp, "gene": gene,
                "human_transcript": tx, "source": "toga2_projection",
                **{k: r.get(k, "NA") for k in
                   ("candidate_id", "projection", "toga_transcript",
                    "toga_gene_label", "toga_status", "orthology_class", "pool",
                    "chrom", "start", "end", "strand",
                    "upstream_gene", "downstream_gene",
                    "distance_to_home_locus_bp", "screen_verdict", "pident",
                    "alignment_coverage", "similarity_score")},
                "sequence_length": len(seq), "cds_status": st,
            })
            status_rows.append({"gene": gene, "human_transcript": tx,
                                "species": sp, "status": st,
                                "toga_status": r.get("toga_status", "NA"),
                                "pool": r.get("pool", "NA"),
                                "sequence_length": len(seq),
                                "similarity_score": r.get("similarity_score", "NA")})

        n_query = len(seen_species)
        fraction = n_query / max(1, n_total)
        passes = fraction >= min_species_fraction and bool(human_seq)

        if not write_directories:
            continue

        gene_dir.mkdir(parents=True, exist_ok=True)
        fasta_path = gene_dir / f"{label}.cds.fa"
        write_fasta(fasta_path, records)
        write_tsv(gene_dir / f"{label}.members.tsv", members, MEMBER_FIELDS)

        tree_out = gene_dir / f"{label}.tree.nwk"
        tips = newick.prune_to_file(tree_path, [n for n, _ in records], tree_out)
        if tips is None:
            tree_out = None
            log(f"  {label}: fewer than 2 tips - no pruned tree written")

        # The selected set is a real copy, not a link: it travels to a cluster
        # or a collaborator intact, and the alignments written into it later
        # stay with the set they belong to.
        sel_gene_dir = sel_fasta = sel_tree = None
        if passes:
            sel_gene_dir = sel_dir / gene / label
            if sel_gene_dir.exists():
                shutil.rmtree(sel_gene_dir)
            sel_gene_dir.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(gene_dir, sel_gene_dir)
            sel_fasta = sel_gene_dir / fasta_path.name
            sel_tree = (sel_gene_dir / Path(tree_out).name) if tree_out else None

        manifest.append({
            "gene": gene, "human_transcript": tx, "directory": str(gene_dir),
            "n_sequences": len(records), "n_query_species": n_query,
            "total_query_species": n_total,
            "fraction_of_species": round(fraction, 4),
            "min_fraction_required": min_species_fraction,
            "human_present": int(bool(human_seq)),
            "tree_tips": len(tips) if tips else 0,
            "n_complete": counts["complete"], "n_partial": counts["partial"],
            "n_fragmented": counts["fragmented"],
            "n_pseudogenized": counts["pseudogenized"],
            "passes_species_threshold": int(passes),
            "fasta": str(fasta_path),
            "newick": str(tree_out) if tree_out else "NA",
            "selected_directory": str(sel_gene_dir) if sel_gene_dir else "NA",
            "selected_fasta": str(sel_fasta) if sel_fasta else "NA",
            "selected_newick": str(sel_tree) if sel_tree else "NA",
        })

    if skipped_species:
        log(f"WARNING: {len(skipped_species)} species absent from the tree were "
            f"dropped: {', '.join(sorted(skipped_species))}")

    status_path = outdir / f"transcript_status_{scope}.tsv"
    write_tsv(status_path, status_rows, STATUS_FIELDS)
    if not write_directories:
        n_sp = len({r["species"] for r in status_rows}) - 1
        log(f"[{scope}] status table only: {len(status_rows)} rows across "
            f"{n_sp} query species -> {status_path.name}")
        return status_path

    write_tsv(outdir / "manifest.tsv", manifest, MANIFEST_FIELDS)
    n_pass = sum(r["passes_species_threshold"] for r in manifest)
    log(f"[{scope}] {len(manifest)} transcript directories under "
        f"{all_dir.name}/; {n_pass} reached {min_species_fraction:.0%} of "
        f"{n_total} query species and were copied to {sel_dir.name}/")
    return outdir / "manifest.tsv"
