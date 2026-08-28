"""Step 4 - region-aware gene copy number variation.

Counting copies from TOGA2 labels alone would be wrong in both directions: a
paralog misfiled under the target's name inflates the count, and a real copy
misfiled under a sibling's name deflates it. So a copy is defined by evidence,
not by a name:

  1. gather every projection that could belong to the gene - its own home
     locus, off-region projections of its human transcripts, and projections
     TOGA2 filed under any paralog (the FAMILY pool);
  2. keep only those the whole-proteome screen attaches to the target gene;
  3. cluster what survives into **non-overlapping** genomic loci, across the
     whole paralog family at once, joining projections only where they actually
     overlap;
  4. give each locus to exactly one gene;
  5. only then rejoin neighbouring loci that went to the *same* gene, so a
     fragmented projection is not counted twice;
  6. record, for every gene, the loci it also has a claim on but did not win.

Step 6 exists because a 1:1 award is sometimes a fiction. A single bat IFITM
locus is the orthologous position of human IFITM2 *and* IFITM3; awarding it to
one and reporting the other as zero copies contradicts the transcript status
table, which quite correctly recovers transcripts of both from that locus. So a
gene's count is reported three ways: `unambiguous_copies` (loci that are its
alone), `shared_copies` (loci it shares with a sibling), and `total_copies`, the
sum. Only the unambiguous set is guaranteed to partition the genome; the shared
count says "this gene is present here too", which is what makes the copy number
and the transcript status agree.

Step 3 has to span the family: clustering per gene would let one stretch of a
tandem array count as a copy of IFITM1 *and* IFITM2 *and* IFITM3.

Step 3 also has to join on overlap alone, and step 5 has to come after the
award. Merging on proximity first collapses a tandem array into a single locus -
the OAS1/OAS2/OAS3 genes sit a few kb apart, so a 10 kb tolerance swallows all
three into one "copy" that is then awarded to one gene, and the other two report
zero copies while their transcripts are plainly present. Distinct neighbouring
genes do not overlap; the fragments of one gene do.

Each copy is then classified so the count can be read at different strictness:
the home copy, additional functional copies, retro/processed copies, and copies
TOGA2 considers lost.
"""
from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .common import (log, orf_report, read_tsv, sanitize, to_float, to_int,
                     write_tsv)

# TOGA2 loss_summary codes that indicate a still-coding projection.
INTACT_STATUSES = {"FI", "I", "PI"}
LOST_STATUSES = {"L", "UL", "M", "N"}

COPY_FIELDS = ["species", "gene", "copy_id", "chrom", "start", "end", "strand",
               "is_home_copy", "copy_class", "attribution", "awarded_to",
               "n_projections", "pools",
               "claimed_by_genes",
               "toga_gene_labels", "best_projection", "best_toga_status",
               "best_orthology_class", "max_cds_bp", "n_exons_best",
               "all_retro", "screen_verdicts", "best_human_gene",
               "cds_integrity", "upstream_gene", "downstream_gene"]

MATRIX_FIELDS = ["species", "gene", "total_copies", "unambiguous_copies",
                 "shared_copies", "shared_with", "copies_excluding_retro",
                 "functional_copies", "home_copy_present",
                 "additional_functional_copies",
                 "retro_copies", "lost_copies", "copy_loci"]


# A stop this close to the end of the supplied CDS is the gene's real stop
# codon; TOGA2 projections routinely carry a few bases past it.
TERMINAL_STOP_FRACTION = 0.95


def _cds_integrity(seq: str) -> str:
    """Coarse CDS integrity of a candidate coding sequence."""
    if not seq:
        return "no_sequence"
    r = orf_report(seq)
    if not r["n_codons"]:
        return "no_sequence"
    if r["stop_codon_index"] is None:
        # runs to the end of the projection without an in-frame stop
        return "truncated_orf"
    if r["orf_fraction"] >= TERMINAL_STOP_FRACTION:
        return "complete_orf" if r["has_start"] else "truncated_orf"
    return "disrupted_orf"


def _cluster(rows: Sequence[dict], slop: int) -> List[List[dict]]:
    """Cluster candidate rows into non-overlapping loci by chromosome.

    Strand is deliberately ignored: an inverted duplicate is still one copy.
    """
    by_chrom: Dict[str, List[dict]] = defaultdict(list)
    for r in rows:
        by_chrom[r["chrom"]].append(r)
    clusters: List[List[dict]] = []
    for chrom, members in by_chrom.items():
        members = sorted(members, key=lambda r: (to_int(r["start"]), to_int(r["end"])))
        cur: List[dict] = []
        hi = 0
        for r in members:
            s, e = to_int(r["start"]), to_int(r["end"])
            if cur and s - slop <= hi:
                cur.append(r)
                hi = max(hi, e)
            else:
                if cur:
                    clusters.append(cur)
                cur, hi = [r], e
        if cur:
            clusters.append(cur)
    return clusters


def _family_of(gene: str, families: Dict[str, set]) -> frozenset:
    """The set of target genes that share a paralog family with `gene`."""
    return frozenset(families.get(gene.upper(), {gene.upper()}))


def run(candidates_tsv, screen_tsv, loci_tsv, sequences_fasta, outdir,
        cluster_slop: int = 10000, min_cds_bp: int = 150,
        families_tsv=None) -> Path:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    from .common import fasta_iter
    seqs = {h.split()[0]: s for h, s in fasta_iter(sequences_fasta)}
    screen = {r["candidate_id"]: r for r in read_tsv(screen_tsv)}
    cands = read_tsv(candidates_tsv, required=["candidate_id", "species", "gene"])

    # Group the target genes into families, so loci can be built family-wide and
    # a shared tandem array is not counted once per member.
    families: Dict[str, set] = {}
    if families_tsv and Path(families_tsv).exists():
        from . import family as family_mod
        families = family_mod.load(families_tsv)
    target_genes = sorted({r["gene"] for r in cands})

    # BLAST family membership is not symmetric: OAS3 can be in OAS1's family
    # without OAS1 being in OAS3's. Left directed, OAS1 and OAS3 get separate
    # components, their loci are clustered twice, and the same stretch of
    # sequence is counted as a copy of both. Symmetrise first, then take
    # connected components, so every locus is clustered exactly once.
    adjacency: Dict[str, set] = {g: {g} for g in target_genes}
    for g in target_genes:
        fam = families.get(g.upper(), set())
        for other in target_genes:
            if other.upper() in fam:
                adjacency[g].add(other)
                adjacency[other].add(g)

    comp_of: Dict[str, frozenset] = {}
    unassigned = set(target_genes)
    while unassigned:
        seed = unassigned.pop()
        component, frontier = {seed}, [seed]
        while frontier:
            for nxt in adjacency[frontier.pop()]:
                if nxt not in component:
                    component.add(nxt)
                    frontier.append(nxt)
        component = frozenset(component)
        for g in component:
            comp_of[g] = component
        unassigned -= component
    for comp in sorted({tuple(sorted(c)) for c in comp_of.values()}):
        if len(comp) > 1:
            log(f"  copy loci clustered jointly for: {', '.join(comp)}")

    # home locus per (species, gene), from the region-restricted search
    home: Dict[Tuple[str, str], Tuple[str, int, int]] = {}
    for r in read_tsv(loci_tsv):
        if r.get("is_home_locus") == "1":
            home[(r["species"], r["gene"])] = (r["chrom"], to_int(r["start"]),
                                               to_int(r["end"]))

    # ---- keep only candidates that really look like the target gene -------
    kept: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    n_dropped_screen = n_dropped_short = 0
    for c in cands:
        s = screen.get(c["candidate_id"])
        verdict = s["screen_verdict"] if s else "NA"
        if verdict not in {"CONSISTENT_WITH_TARGET_GENE",
                           "CONSISTENT_BY_SYNTENY_WITHIN_FAMILY"}:
            n_dropped_screen += 1
            continue
        if to_int(c.get("cds_bp")) < min_cds_bp:
            n_dropped_short += 1
            continue
        kept[(c["species"], c["gene"])].append({**c, **{
            "screen_verdict": verdict,
            "best_human_gene": s["best_human_gene"] if s else "NA"}})

    log(f"copy number: {n_dropped_screen} candidate(s) dropped by the paralog "
        f"screen, {n_dropped_short} dropped as too short (<{min_cds_bp} bp)")

    copy_rows: List[dict] = []
    matrix_rows: List[dict] = []
    species_seen = {sp for sp, _ in kept}
    genes_seen = {g for _, g in kept}

    # Cluster once per (species, paralog family), so loci never overlap and each
    # is awarded to a single gene.
    by_family: Dict[Tuple[str, frozenset], List[dict]] = defaultdict(list)
    for (species, gene), rows in kept.items():
        by_family[(species, comp_of.get(gene, frozenset({gene})))].extend(rows)

    owned: Dict[Tuple[str, str], List[List[dict]]] = defaultdict(list)
    shared: Dict[Tuple[str, str], List[tuple]] = defaultdict(list)
    contested_rows: List[dict] = []
    for (species, fam), rows in sorted(by_family.items(), key=lambda kv: (kv[0][0], sorted(kv[0][1]))):
        # Overlap only. Proximity merging happens after the award, below.
        for members in _cluster(rows, 0):
            chrom = members[0]["chrom"]
            lo = min(to_int(m["start"]) for m in members)
            hi = max(to_int(m["end"]) for m in members)
            claimants = sorted({m["gene"] for m in members})

            # The gene whose own syntenic locus this is wins; otherwise the gene
            # with the strongest sequence evidence at this locus.
            syntenic = [g for g in claimants
                        if home.get((species, g))
                        and home[(species, g)][0] == chrom
                        and hi > home[(species, g)][1] and lo < home[(species, g)][2]]
            if len(syntenic) == 1:
                owner, basis = syntenic[0], "syntenic_home_locus"
            else:
                pool = syntenic or claimants

                def evidence(g):
                    hits = [to_float(screen.get(m["candidate_id"], {})
                                     .get("target_gene_bitscore", 0))
                            for m in members if m["gene"] == g]
                    return (max(hits) if hits else 0.0, g)

                owner = max(pool, key=evidence)
                basis = ("competing_syntenic_loci" if len(syntenic) > 1
                         else "best_sequence_evidence")

            if len(claimants) > 1:
                contested_rows.append({
                    "species": species, "chrom": chrom, "start": lo, "end": hi,
                    "claimed_by": ",".join(claimants), "awarded_to": owner,
                    "basis": basis, "n_projections": len(members)})
                # The losers keep a recorded claim on this locus, so a gene that
                # is genuinely present here is not reported as having no copy.
                for other in claimants:
                    if other != owner:
                        shared[(species, other)].append(
                            (chrom, lo, hi, owner,
                             [m for m in members if m["gene"] == other] or members))

            owned[(species, owner)].append([m for m in members])

    # Now that every locus belongs to one gene, rejoin neighbouring loci of the
    # same gene: those are fragments of one copy, not separate copies.
    #
    # The walk is done per species over every gene's loci at once, in coordinate
    # order, and two fragments are joined only when they are *adjacent in that
    # order*. Merging a gene's loci in isolation would let one gene's copy span
    # a neighbour's copy sitting between them in an interleaved array, which
    # reintroduces exactly the overlap this design exists to prevent.
    per_species: Dict[str, List[Tuple[str, str, int, int, List[dict]]]] = defaultdict(list)
    for (species, gene), clusters in owned.items():
        for c in clusters:
            per_species[species].append(
                (gene, c[0]["chrom"], min(to_int(m["start"]) for m in c),
                 max(to_int(m["end"]) for m in c), c))

    owned = defaultdict(list)
    for species, loci in per_species.items():
        loci.sort(key=lambda x: (x[1], x[2], x[3]))
        cur = None
        for gene, chrom, lo, hi, members in loci:
            same_run = (cur is not None and cur[0] == gene and cur[1] == chrom
                        and lo - cur[3] <= cluster_slop)
            if same_run:
                cur = (gene, chrom, cur[2], max(cur[3], hi), cur[4] + members)
            else:
                if cur is not None:
                    owned[(species, cur[0])].append(cur[4])
                cur = (gene, chrom, lo, hi, list(members))
        if cur is not None:
            owned[(species, cur[0])].append(cur[4])

    for (species, gene) in sorted(set(list(kept.keys())) | set(owned.keys())):
        home_span = home.get((species, gene))
        clusters = owned.get((species, gene), [])

        copies = []
        for n, members in enumerate(sorted(clusters,
                                           key=lambda m: (m[0]["chrom"],
                                                          to_int(m[0]["start"]))), 1):
            chrom = members[0]["chrom"]
            lo = min(to_int(m["start"]) for m in members)
            hi = max(to_int(m["end"]) for m in members)
            best = max(members, key=lambda m: (to_int(m["cds_bp"]),
                                               to_int(m["seq_length"])))
            all_retro = all(m.get("is_retro") == "1" for m in members)
            statuses = {m.get("toga_status", "NA") for m in members}
            any_intact = bool(statuses & INTACT_STATUSES)
            integrity = _cds_integrity(seqs.get(best["candidate_id"], ""))

            is_home = bool(home_span and home_span[0] == chrom
                           and hi > home_span[1] and lo < home_span[2])

            if all_retro:
                copy_class = "retro_or_processed"
            elif any_intact and integrity in {"complete_orf", "truncated_orf"}:
                copy_class = "functional"
            elif any_intact:
                copy_class = "coding_but_disrupted"
            else:
                copy_class = "lost_or_uncertain"

            copies.append({
                "species": species, "gene": gene,
                "copy_id": sanitize(f"{species}|{gene}|copy{n}"),
                "chrom": chrom, "start": lo, "end": hi,
                "strand": best.get("strand", "."),
                "is_home_copy": int(is_home), "copy_class": copy_class,
                "attribution": "sole", "awarded_to": gene,
                "n_projections": len(members),
                "pools": ",".join(sorted({m["pool"] for m in members})),
                "claimed_by_genes": ",".join(sorted({m["gene"] for m in members})),
                "toga_gene_labels": ",".join(sorted({m["toga_gene_label"] for m in members})),
                "best_projection": best["projection"],
                "best_toga_status": best.get("toga_status", "NA"),
                "best_orthology_class": best.get("orthology_class", "NA"),
                "max_cds_bp": to_int(best["cds_bp"]),
                "n_exons_best": to_int(best.get("exon_count")),
                "all_retro": int(all_retro),
                "screen_verdicts": ",".join(sorted({m["screen_verdict"] for m in members})),
                "best_human_gene": best.get("best_human_gene", "NA"),
                "cds_integrity": integrity,
                "upstream_gene": best.get("upstream_gene", "NA"),
                "downstream_gene": best.get("downstream_gene", "NA"),
            })

        # loci this gene shares with a sibling paralog
        shared_here = []
        for n, (chrom, lo, hi, owner, members) in enumerate(
                sorted(shared.get((species, gene), []), key=lambda x: (x[0], x[1])), 1):
            best = max(members, key=lambda m: (to_int(m["cds_bp"]),
                                               to_int(m["seq_length"])))
            shared_here.append({
                "species": species, "gene": gene,
                "copy_id": sanitize(f"{species}|{gene}|shared{n}"),
                "chrom": chrom, "start": lo, "end": hi,
                "strand": best.get("strand", "."),
                "is_home_copy": 0, "copy_class": "shared_with_paralog",
                "attribution": "shared", "awarded_to": owner,
                "n_projections": len(members),
                "pools": ",".join(sorted({m["pool"] for m in members})),
                "claimed_by_genes": ",".join(sorted({m["gene"] for m in members})),
                "toga_gene_labels": ",".join(sorted({m["toga_gene_label"] for m in members})),
                "best_projection": best["projection"],
                "best_toga_status": best.get("toga_status", "NA"),
                "best_orthology_class": best.get("orthology_class", "NA"),
                "max_cds_bp": to_int(best["cds_bp"]),
                "n_exons_best": to_int(best.get("exon_count")),
                "all_retro": 0,
                "screen_verdicts": ",".join(sorted({m["screen_verdict"] for m in members})),
                "best_human_gene": best.get("best_human_gene", "NA"),
                "cds_integrity": _cds_integrity(seqs.get(best["candidate_id"], "")),
                "upstream_gene": best.get("upstream_gene", "NA"),
                "downstream_gene": best.get("downstream_gene", "NA"),
            })

        copy_rows.extend(copies)
        copy_rows.extend(shared_here)
        functional = [c for c in copies if c["copy_class"] == "functional"]
        matrix_rows.append({
            "species": species, "gene": gene,
            "total_copies": len(copies) + len(shared_here),
            "unambiguous_copies": len(copies),
            "shared_copies": len(shared_here),
            "shared_with": ",".join(sorted({c["awarded_to"] for c in shared_here})) or "NONE",
            # Retro/processed copies are real copies, but they are not paralog
            # expansion and one species can carry many; kept separate so a
            # figure can use either count.
            "copies_excluding_retro":
                sum(1 for c in copies + shared_here
                    if c["copy_class"] != "retro_or_processed"),
            "functional_copies": len(functional),
            "home_copy_present": int(any(c["is_home_copy"] for c in copies)),
            "additional_functional_copies":
                sum(1 for c in functional if not c["is_home_copy"]),
            "retro_copies": sum(1 for c in copies if c["copy_class"] == "retro_or_processed"),
            "lost_copies": sum(1 for c in copies
                               if c["copy_class"] in {"lost_or_uncertain",
                                                      "coding_but_disrupted"}),
            "copy_loci": ";".join(f"{c['chrom']}:{c['start']}-{c['end']}" for c in copies),
        })

    # Species x gene with no surviving candidate at all still need a zero row.
    all_genes = sorted(genes_seen)
    present = {(r["species"], r["gene"]) for r in matrix_rows}
    for species in sorted(species_seen):
        for gene in all_genes:
            if (species, gene) not in present:
                matrix_rows.append({"species": species, "gene": gene,
                                    "total_copies": 0,
                                    "unambiguous_copies": 0,
                                    "shared_copies": 0, "shared_with": "NONE",
                                    "copies_excluding_retro": 0,
                                    "functional_copies": 0,
                                    "home_copy_present": 0,
                                    "additional_functional_copies": 0,
                                    "retro_copies": 0, "lost_copies": 0,
                                    "copy_loci": "NONE"})
    matrix_rows.sort(key=lambda r: (r["gene"], r["species"]))

    write_tsv(outdir / "gene_copies.tsv", copy_rows, COPY_FIELDS)
    write_tsv(outdir / "contested_loci.tsv", contested_rows,
              ["species", "chrom", "start", "end", "claimed_by", "awarded_to",
               "basis", "n_projections"])
    if contested_rows:
        log(f"copy number: {len(contested_rows)} locus/loci were claimed by more "
            f"than one family member and awarded to one - see contested_loci.tsv")
    _check_non_overlapping(copy_rows, outdir)
    write_tsv(outdir / "copy_number_matrix.tsv", matrix_rows, MATRIX_FIELDS)
    _write_wide(matrix_rows, outdir / "copy_number_wide.tsv")

    log(f"copy number: {len(copy_rows)} copies across {len(species_seen)} species "
        f"and {len(all_genes)} genes")
    return outdir / "copy_number_matrix.tsv"


def _check_non_overlapping(copy_rows: List[dict], outdir: Path) -> None:
    """Copy loci within a species must not overlap, whichever gene owns them."""
    by_species: Dict[str, List[dict]] = defaultdict(list)
    for r in copy_rows:
        # Shared rows are a second attribution of a locus already counted under
        # its owner; only the awarded set is required to partition the genome.
        if r.get("attribution") == "sole":
            by_species[r["species"]].append(r)
    overlaps = []
    for species, rows in by_species.items():
        rows = sorted(rows, key=lambda r: (r["chrom"], to_int(r["start"])))
        for a, b in zip(rows, rows[1:]):
            if a["chrom"] == b["chrom"] and to_int(b["start"]) < to_int(a["end"]):
                overlaps.append({"species": species, "chrom": a["chrom"],
                                 "copy_a": a["copy_id"], "gene_a": a["gene"],
                                 "copy_b": b["copy_id"], "gene_b": b["gene"],
                                 "a_end": a["end"], "b_start": b["start"]})
    write_tsv(outdir / "overlapping_copy_loci.tsv", overlaps,
              ["species", "chrom", "copy_a", "gene_a", "copy_b", "gene_b",
               "a_end", "b_start"])
    if overlaps:
        log(f"WARNING: {len(overlaps)} pair(s) of overlapping copy loci remain")
    else:
        log("copy loci are non-overlapping within every species")


def _write_wide(rows: Sequence[dict], path) -> None:
    """species x gene matrix of functional copy counts, for quick inspection."""
    genes = sorted({r["gene"] for r in rows})
    by_species: Dict[str, Dict[str, int]] = defaultdict(dict)
    for r in rows:
        by_species[r["species"]][r["gene"]] = r["functional_copies"]
    out = []
    for species in sorted(by_species):
        row = {"species": species}
        row.update({g: by_species[species].get(g, 0) for g in genes})
        out.append(row)
    write_tsv(path, out, ["species"] + genes)
