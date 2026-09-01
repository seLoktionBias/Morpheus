"""Step 3b - does the search scope change which transcript you recover?

Two assignments exist for every human transcript in every species: one confined
to the gene's own syntenic locus, one allowed anywhere the gene occurs. Most of
the time they choose the same model and the distinction does not matter. Where
they differ it matters a great deal, and that is what this table records.

The motivating case is OAS1 in Phyllostomus discolor. The transcript carrying
the C-terminal CaaX domain is degraded at the syntenic locus; a copy on another
scaffold carries it at ~98% identity to human. Restricted to the locus, the
animal reads as having lost the domain. Unrestricted, it clearly has not. Both
statements are true of different questions, so both are reported and the
disagreement is made explicit rather than resolved silently.

Outcomes per (species, gene, human transcript):

  SAME_MODEL              both scopes chose the same projection
  DIFFERENT_LOCUS         the unrestricted scope found a better model elsewhere
  ONLY_UNRESTRICTED       nothing acceptable in the locus; found outside it
  ONLY_REGION_RESTRICTED  the locus model lost its slot once rivals were allowed
  NEITHER                 no acceptable model in either scope
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

from .common import log, read_tsv, to_float, to_int, write_tsv

COMPARE_FIELDS = [
    "species", "gene", "human_transcript", "outcome",
    "region_candidate", "region_projection", "region_chrom", "region_start",
    "region_end", "region_pident", "region_score", "region_cds_status",
    "unrestricted_candidate", "unrestricted_projection", "unrestricted_chrom",
    "unrestricted_start", "unrestricted_end", "unrestricted_pident",
    "unrestricted_score", "unrestricted_cds_status", "unrestricted_pool",
    "distance_to_home_locus_bp", "same_chrom",
    "pident_gain", "score_gain", "length_gain_bp",
    "region_upstream_gene", "region_downstream_gene",
    "unrestricted_upstream_gene", "unrestricted_downstream_gene",
]

SUMMARY_FIELDS = ["gene", "human_transcript", "n_species", "n_same_model",
                  "n_different_locus", "n_only_unrestricted",
                  "n_only_region_restricted", "n_neither",
                  "n_scope_changes", "median_pident_gain"]


def _index(rows: List[dict]) -> Dict[Tuple[str, str, str], dict]:
    return {(r["species"], r["gene"], r["human_transcript"]): r
            for r in rows if r.get("assignment_status") == "ASSIGNED"}


def build(region_tsv, unrestricted_tsv, outdir, sequences_fasta=None,
          policy: str = "structure_aware") -> Path:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    region_rows = read_tsv(region_tsv)
    unres_rows = read_tsv(unrestricted_tsv)
    region = _index(region_rows)
    unres = _index(unres_rows)

    status: Dict[str, str] = {}
    if sequences_fasta and Path(sequences_fasta).exists():
        from .common import fasta_iter
        from .deliverables import cds_status
        status = {h.split()[0]: cds_status(s) for h, s in fasta_iter(sequences_fasta)}

    keys = sorted(set(region) | set(unres)
                  | {(r["species"], r["gene"], r["human_transcript"])
                     for r in region_rows + unres_rows})

    rows: List[dict] = []
    for key in keys:
        species, gene, tx = key
        a, b = region.get(key), unres.get(key)
        if a is None and b is None:
            outcome = "NEITHER"
        elif b is None:
            outcome = "ONLY_REGION_RESTRICTED"
        elif a is None:
            outcome = "ONLY_UNRESTRICTED"
        elif a["candidate_id"] == b["candidate_id"]:
            outcome = "SAME_MODEL"
        else:
            outcome = "DIFFERENT_LOCUS"

        def field(rec, name, default="NA"):
            return rec.get(name, default) if rec else default

        same_chrom = ("NA" if not (a and b)
                      else int(a.get("chrom") == b.get("chrom")))
        rows.append({
            "species": species, "gene": gene, "human_transcript": tx,
            "outcome": outcome,
            "region_candidate": field(a, "candidate_id"),
            "region_projection": field(a, "projection"),
            "region_chrom": field(a, "chrom"), "region_start": field(a, "start"),
            "region_end": field(a, "end"), "region_pident": field(a, "pident"),
            "region_score": field(a, "similarity_score"),
            "region_cds_status": status.get(field(a, "candidate_id"), "not_found"),
            "unrestricted_candidate": field(b, "candidate_id"),
            "unrestricted_projection": field(b, "projection"),
            "unrestricted_chrom": field(b, "chrom"),
            "unrestricted_start": field(b, "start"),
            "unrestricted_end": field(b, "end"),
            "unrestricted_pident": field(b, "pident"),
            "unrestricted_score": field(b, "similarity_score"),
            "unrestricted_cds_status": status.get(field(b, "candidate_id"), "not_found"),
            "unrestricted_pool": field(b, "pool"),
            "distance_to_home_locus_bp": field(b, "distance_to_home_locus_bp"),
            "same_chrom": same_chrom,
            "pident_gain": (round(to_float(b["pident"]) - to_float(a["pident"]), 3)
                            if a and b else "NA"),
            "score_gain": (round(to_float(b["similarity_score"])
                                 - to_float(a["similarity_score"]), 4)
                           if a and b else "NA"),
            "length_gain_bp": (to_int(b.get("seq_length")) - to_int(a.get("seq_length"))
                               if a and b else "NA"),
            "region_upstream_gene": field(a, "upstream_gene"),
            "region_downstream_gene": field(a, "downstream_gene"),
            "unrestricted_upstream_gene": field(b, "upstream_gene"),
            "unrestricted_downstream_gene": field(b, "downstream_gene"),
        })

    write_tsv(outdir / f"scope_comparison__{policy}.tsv", rows, COMPARE_FIELDS)

    # per-transcript rollup, so a gene with a systematic difference stands out
    per_tx: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    for r in rows:
        per_tx[(r["gene"], r["human_transcript"])].append(r)

    import statistics
    summary = []
    for (gene, tx), group in sorted(per_tx.items()):
        counts = defaultdict(int)
        for r in group:
            counts[r["outcome"]] += 1
        gains = [to_float(r["pident_gain"]) for r in group
                 if r["outcome"] == "DIFFERENT_LOCUS" and r["pident_gain"] != "NA"]
        changed = (counts["DIFFERENT_LOCUS"] + counts["ONLY_UNRESTRICTED"]
                   + counts["ONLY_REGION_RESTRICTED"])
        summary.append({
            "gene": gene, "human_transcript": tx, "n_species": len(group),
            "n_same_model": counts["SAME_MODEL"],
            "n_different_locus": counts["DIFFERENT_LOCUS"],
            "n_only_unrestricted": counts["ONLY_UNRESTRICTED"],
            "n_only_region_restricted": counts["ONLY_REGION_RESTRICTED"],
            "n_neither": counts["NEITHER"],
            "n_scope_changes": changed,
            "median_pident_gain": round(statistics.median(gains), 3) if gains else "NA",
        })
    write_tsv(outdir / f"scope_comparison_summary__{policy}.tsv", summary,
              SUMMARY_FIELDS)

    totals = defaultdict(int)
    for r in rows:
        totals[r["outcome"]] += 1
    log(f"[{policy}] scope comparison: " +
        ", ".join(f"{k}={v}" for k, v in sorted(totals.items())))
    n_recovered = sum(1 for r in rows
                      if r["outcome"] in ("DIFFERENT_LOCUS", "ONLY_UNRESTRICTED")
                      and r["unrestricted_cds_status"] == "complete"
                      and r["region_cds_status"] != "complete")
    log(f"  {n_recovered} transcript(s) are complete only when the search leaves "
        f"the gene's own locus")
    return outdir / f"scope_comparison__{policy}.tsv"
