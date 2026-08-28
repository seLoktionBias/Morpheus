"""Step 2 - region-restricted transcript search in each bat genome.

For one gene in one bat genome we build two candidate pools.

  IN_REGION  Every TOGA2 projection sharing at least one exon with the
             projected longest human isoform, regardless of the gene name TOGA2
             attached to it. This is the region-restricted search: it recovers
             transcripts TOGA2 named after a paralog, and it refuses to trust
             the name alone.

             The test is exon overlap, not span overlap. A smaller gene often
             sits inside a large gene's intron, and TOGA2 sometimes files it
             under the host gene's name; anything keyed on spans sweeps it in as
             a transcript of the host. Requiring a shared exon means an intronic
             neighbour cannot qualify however deeply it is nested.

  OFF_REGION Projections of this gene's human transcripts that landed in a
             different genomic locus. These are kept as candidates because the
             most similar CDS is not always in the expected region - the OAS1
             projection in Phyllostomus discolor is the motivating case.

  FAMILY     Projections TOGA2 labelled with any *paralog* of the gene, anywhere
             in the genome. These are not transcript candidates; they exist so
             the copy-number step can find a real gene copy that TOGA2 filed
             under a sibling gene's name. They are screened against the whole
             human proteome before being counted.

The home locus is anchored on the projection of the human longest isoform, so
"expected region" is defined by the data rather than by TOGA2's labelling.

`query_annotation.bed` is not the whole annotation. TOGA2 also writes
`processed_pseudogenes.bed` and `fragmented_annotation.bed`, and a functionally
important copy can sit in either. In Phyllostomus discolor the OAS1 transcript
carrying the C-terminal CaaX motif is a retrocopy on a different scaffold, filed
as a processed pseudogene; the syntenic locus has no CaaX-bearing model at all.
Reading only the main annotation would report that animal as having lost the
domain. All three files are therefore read, deduplicated, and tagged with the
file each projection came from.
"""
from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .common import (Bed12, ProjectionName, clean_nt, die, fasta_iter, log,
                     merge_intervals, opn, parse_bed12, parse_bed6,
                     parse_projection, read_tsv, reciprocal_overlap, sanitize,
                     species_from_dirname, strip_version, write_fasta,
                     write_tsv)

# TOGA2 annotation files, in order of precedence when the same projection
# appears in more than one.
ANNOTATION_FILES = [
    ("query_annotation", ("query_annotation.bed.gz", "query_annotation.bed")),
    ("processed_pseudogene", ("processed_pseudogenes.bed.gz",
                              "processed_pseudogenes.bed")),
    ("fragmented", ("fragmented_annotation.bed.gz", "fragmented_annotation.bed")),
]

CANDIDATE_FIELDS = [
    "species", "gene", "candidate_id", "projection", "toga_transcript",
    "toga_gene_label", "chain", "is_retro", "annotation_source",
    "toga_status", "orthology_class",
    "chrom", "start", "end", "strand", "exon_count", "cds_bp", "seq_length",
    "locus_id", "pool", "distance_to_home_locus_bp",
    "exon_structure", "matched_longest_exons", "n_matched_longest_exons",
    "n_bat_novel_exons", "upstream_gene", "downstream_gene",
]

# One row per candidate exon, for the exon-model figures.
EXON_FIELDS = ["species", "gene", "candidate_id", "projection",
               "toga_transcript", "toga_gene_label", "pool", "chrom", "strand",
               "start", "end", "length_bp", "transcript_exon_rank",
               "exon_label", "exon_label_type", "is_anchor",
               "previous_gene", "next_gene", "gene_order"]

NESTED_FIELDS = ["species", "gene", "projection", "toga_gene_label", "chrom",
                 "start", "end", "strand", "exon_count", "anchor_projection",
                 "anchor_start", "anchor_end", "reason"]

LOCUS_FIELDS = [
    "species", "gene", "locus_id", "chrom", "start", "end", "strand",
    "is_home_locus", "n_projections", "n_distinct_human_transcripts",
    "anchor_projection", "anchor_start", "anchor_end", "anchor_exon_count",
    "n_in_region_after_exon_test", "n_dropped_no_shared_exon",
    "max_cds_bp", "any_retro", "toga_gene_labels",
    "upstream_gene", "downstream_gene",
]


# --------------------------------------------------------------------------
# TOGA2 side tables
# --------------------------------------------------------------------------


def read_loss_summary(path) -> Dict[str, str]:
    """{projection_label: status} for PROJECTION-level rows."""
    out: Dict[str, str] = {}
    if not Path(path).exists():
        return out
    with opn(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[0] == "PROJECTION":
                out[f[1]] = f[2]
    return out


def read_orthology(path) -> Dict[str, str]:
    """{query_projection_label: orthology_class}."""
    out: Dict[str, str] = {}
    if not Path(path).exists():
        return out
    with opn(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            qi, ci = header.index("q_transcript"), header.index("orthology_class")
        except ValueError:
            return out
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(qi, ci):
                out[f[qi]] = f[ci]
    return out


# --------------------------------------------------------------------------
# locus clustering
# --------------------------------------------------------------------------


def cluster_loci(beds: Sequence[Bed12], slop: int = 0
                 ) -> List[Tuple[str, str, int, int, List[Bed12]]]:
    """Group projections into non-overlapping loci by chromosome and strand."""
    by_key: Dict[Tuple[str, str], List[Bed12]] = defaultdict(list)
    for b in beds:
        by_key[(b.chrom, b.strand)].append(b)

    loci = []
    for (chrom, strand), members in by_key.items():
        members = sorted(members, key=lambda b: (b.start, b.end))
        cur: List[Bed12] = []
        cur_lo = cur_hi = 0
        for b in members:
            if cur and b.start - slop <= cur_hi:
                cur.append(b)
                cur_hi = max(cur_hi, b.end)
            else:
                if cur:
                    loci.append((chrom, strand, cur_lo, cur_hi, cur))
                cur, cur_lo, cur_hi = [b], b.start, b.end
        if cur:
            loci.append((chrom, strand, cur_lo, cur_hi, cur))
    loci.sort(key=lambda x: (x[0], x[2]))
    return loci


_ENSG = re.compile(r"^(ENSG\d+)")


def _readable_gene(name: str, ensg_map: Dict[str, str]) -> str:
    """TOGA2 names gene regions by Ensembl ID; show the symbol where known."""
    m = _ENSG.match(str(name))
    return ensg_map.get(m.group(1), name) if m else name


def neighbouring_genes(gene_bed: Sequence[dict], chrom: str, lo: int, hi: int,
                       strand: str, ensg_map: Optional[Dict[str, str]] = None
                       ) -> Tuple[str, str]:
    """Nearest annotated TOGA2 gene region on each side, in transcription order."""
    same = [r for r in gene_bed if r["chrom"] == chrom]
    if not same:
        return "NA", "NA"
    left = [r for r in same if r["end"] <= lo]
    right = [r for r in same if r["start"] >= hi]
    prev_gene = max(left, key=lambda r: r["end"])["name"] if left else "NA"
    next_gene = min(right, key=lambda r: r["start"])["name"] if right else "NA"
    if ensg_map:
        prev_gene = _readable_gene(prev_gene, ensg_map)
        next_gene = _readable_gene(next_gene, ensg_map)
    return (prev_gene, next_gene) if strand != "-" else (next_gene, prev_gene)


# --------------------------------------------------------------------------
# exon-label projection
# --------------------------------------------------------------------------


def _shares_exon(a: Sequence[Tuple[int, int]],
                 b: Sequence[Tuple[int, int]], min_bp: int = 1) -> bool:
    """Do two exon block lists overlap in at least one exon?"""
    for s1, e1 in a:
        for s2, e2 in b:
            if min(e1, e2) - max(s1, s2) >= min_bp:
                return True
    return False


def project_exon_labels(candidates: Sequence[Bed12],
                        structure_by_tx: Dict[str, List[str]]
                        ) -> Dict[Tuple[int, int], str]:
    """Map projected bat exon coordinates onto human exon labels.

    A bat exon only inherits a human label when its projected coordinates match
    a human exon projection exactly. Everything else becomes bat_novel_exonN
    later, so that a shifted splice site is never silently called "exon3".
    """
    coord_to_label: Dict[Tuple[int, int], str] = {}
    for bed in candidates:
        labels = structure_by_tx.get(strip_version(parse_projection(bed.name).transcript))
        if not labels:
            continue
        blocks = bed.blocks_transcript_order
        for i, label in enumerate(labels):
            if i >= len(blocks):
                break
            coord_to_label.setdefault(blocks[i], label)
    return coord_to_label


# --------------------------------------------------------------------------
# per-species driver
# --------------------------------------------------------------------------


def _toga_file(run_dir: Path, *names: str) -> Optional[Path]:
    for n in names:
        p = run_dir / n
        if p.exists() and p.stat().st_size > 0:
            return p
    return None


def run_species(run_dir, reference_dir, outdir, locus_slop: int = 0,
                family_map: Optional[Dict[str, set]] = None,
                tree_tips: Optional[Sequence[str]] = None) -> Path:
    """Extract candidates for every gene of interest from one TOGA2 run."""
    run_dir = Path(run_dir)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    species = species_from_dirname(run_dir.name, tree_tips)

    fa_path = _toga_file(run_dir, "nucleotide.fa.gz", "nucleotide.fa")
    if fa_path is None:
        die(f"{run_dir.name}: no nucleotide.fa(.gz)")

    gene_bed = parse_bed6(_toga_file(run_dir, "query_genes.bed.gz", "query_genes.bed") or "")
    status_by_proj = read_loss_summary(_toga_file(run_dir, "loss_summary.tsv.gz",
                                                  "loss_summary.tsv") or "")
    orth_by_proj = read_orthology(_toga_file(run_dir, "orthology_classification.tsv.gz",
                                             "orthology_classification.tsv") or "")

    ensg_map = {r["gene_id"]: r["gene"] for r in
                read_tsv(Path(reference_dir) / "ensembl_gene_id_to_symbol.tsv")}
    ref_tx = read_tsv(Path(reference_dir) / "transcripts.tsv",
                      required=["gene", "transcript_id", "exon_structure"])
    ref_ctx = {r["gene"].upper(): r
               for r in read_tsv(Path(reference_dir) / "gene_context.tsv")
               if r.get("status") == "OK"}

    seeds_by_gene: Dict[str, set] = defaultdict(set)
    structure_by_gene: Dict[str, Dict[str, List[str]]] = defaultdict(dict)
    for r in ref_tx:
        g = r["gene"].upper()
        tid = strip_version(r["transcript_id"])
        seeds_by_gene[g].add(tid)
        structure_by_gene[g][tid] = [x for x in r["exon_structure"].split(",") if x]
        # Transcripts collapsed into this representative are still seeds, so
        # their TOGA2 projections are found when anchoring the locus.
        for extra in r.get("identical_model_transcripts", "NONE").split(","):
            extra = strip_version(extra)
            if extra and extra != "NONE":
                seeds_by_gene[g].add(extra)
                # identical exon model, so the label chain transfers
                structure_by_gene[g][extra] = structure_by_gene[g][tid]
        for extra in r.get("identical_cds_other_model_transcripts", "NONE").split(","):
            extra = strip_version(extra)
            if extra and extra != "NONE":
                # same CDS but a different exon model: seed it so its
                # projection is found, but do not lend it a label chain that
                # does not describe its exons
                seeds_by_gene[g].add(extra)

    # Read every annotation file TOGA2 wrote, not just the main one. The same
    # projection can appear in several; the first file to list it wins, and its
    # origin is recorded so a retrocopy is never mistaken for a locus model.
    all_beds: List[Bed12] = []
    source_of: Dict[Tuple[str, str, int, int], str] = {}
    counts: Dict[str, int] = {}
    for source, names in ANNOTATION_FILES:
        path = _toga_file(run_dir, *names)
        if path is None:
            continue
        added = 0
        for b in parse_bed12(path):
            key = (b.name, b.chrom, b.start, b.end)
            if key in source_of:
                continue
            source_of[key] = source
            all_beds.append(b)
            added += 1
        counts[source] = added
    if not all_beds:
        die(f"{run_dir.name}: no usable annotation BED found")
    log(f"{species}: " + ", ".join(f"{k}={v}" for k, v in counts.items()))

    # Index projections by versionless human transcript id and by the gene name
    # TOGA2 attached, once per species.
    beds_by_seed: Dict[str, List[Bed12]] = defaultdict(list)
    beds_by_toga_gene: Dict[str, List[Bed12]] = defaultdict(list)
    for b in all_beds:
        proj = parse_projection(b.name)
        beds_by_seed[proj.transcript_base].append(b)
        beds_by_toga_gene[proj.gene.upper()].append(b)
    family_map = family_map or {}

    candidate_rows: List[dict] = []
    locus_rows: List[dict] = []
    nested_rows: List[dict] = []
    exon_rows: List[dict] = []
    # One projection can be a candidate for more than one gene (a paralog's
    # projection is IN_REGION for its own gene and FAMILY for its sibling), so
    # a FASTA key maps to a list of candidate ids.
    wanted_seq_keys: Dict[str, List[str]] = defaultdict(list)

    for gene_upper, ctx in sorted(ref_ctx.items()):
        gene = ctx["gene"]
        seeds = seeds_by_gene.get(gene_upper, set())
        if not seeds:
            continue

        seed_beds = [b for s in seeds for b in beds_by_seed.get(s, [])]
        if not seed_beds:
            locus_rows.append({"species": species, "gene": gene,
                               "locus_id": "NA", "is_home_locus": 0,
                               "n_projections": 0,
                               "anchor_projection": "NO_PROJECTION_FOUND"})
            continue

        loci = cluster_loci(seed_beds, slop=locus_slop)
        longest_base = strip_version(ctx["longest_transcript"])

        # Home locus = where the longest human isoform projected. If the longest
        # isoform is absent, fall back to the locus carrying the most distinct
        # human transcripts, then the largest CDS.
        def locus_rank(locus):
            _, _, _, _, members = locus
            has_longest = any(parse_projection(m.name).transcript_base == longest_base
                              for m in members)
            n_tx = len({parse_projection(m.name).transcript_base for m in members})
            return (has_longest, n_tx, max(m.cds_bp for m in members))

        home = max(loci, key=locus_rank)
        home_chrom, home_strand, home_lo, home_hi, home_members = home
        # The anchor is the projection of the human longest isoform where one
        # reached this locus, else the largest projection in it.
        home_anchor = max(
            home_members,
            key=lambda m: (parse_projection(m.name).transcript_base == longest_base,
                           m.cds_bp, m.block_count))

        for idx, (chrom, strand, lo, hi, members) in enumerate(loci, 1):
            is_home = (chrom, strand, lo, hi) == (home_chrom, home_strand, home_lo, home_hi)
            up, down = neighbouring_genes(gene_bed, chrom, lo, hi, strand, ensg_map)
            anchor = max(members, key=lambda m: (
                parse_projection(m.name).transcript_base == longest_base, m.cds_bp))
            locus_rows.append({
                "species": species, "gene": gene,
                "locus_id": f"{species}|{gene}|locus{idx}",
                "chrom": chrom, "start": lo, "end": hi, "strand": strand,
                "is_home_locus": int(is_home), "n_projections": len(members),
                "n_distinct_human_transcripts":
                    len({parse_projection(m.name).transcript_base for m in members}),
                "anchor_projection": anchor.name,
                "anchor_start": anchor.start, "anchor_end": anchor.end,
                "anchor_exon_count": anchor.block_count,
                "max_cds_bp": max(m.cds_bp for m in members),
                "any_retro": int(any(parse_projection(m.name).is_retro for m in members)),
                "toga_gene_labels": ",".join(sorted(
                    {parse_projection(m.name).gene for m in members})),
                "upstream_gene": up, "downstream_gene": down,
            })

        # ---- IN_REGION: shares an exon with the projected longest isoform --
        # `home_anchor` is the projection of the human longest isoform, so the
        # locus is defined by that transcript's exons rather than by a window.
        anchor_blocks = home_anchor.blocks
        in_region = [b for b in all_beds
                     if b.chrom == home_chrom and b.strand == home_strand
                     and b.end > home_lo and b.start < home_hi
                     and _shares_exon(b.blocks, anchor_blocks)]
        in_region_names = {b.name for b in in_region}

        # Projections that fall inside the locus span but share no exon with the
        # anchor: typically a smaller gene nested in one of its introns.
        span_only = [b for b in all_beds
                     if b.chrom == home_chrom and b.strand == home_strand
                     and b.end > home_lo and b.start < home_hi
                     and b.name not in in_region_names]
        for b in span_only:
            proj = parse_projection(b.name)
            nested_rows.append({
                "species": species, "gene": gene,
                "projection": proj.raw, "toga_gene_label": proj.gene,
                "chrom": b.chrom, "start": b.start, "end": b.end,
                "strand": b.strand, "exon_count": b.block_count,
                "anchor_projection": home_anchor.name,
                "anchor_start": home_anchor.start, "anchor_end": home_anchor.end,
                "reason": "inside_the_locus_span_but_shares_no_exon_with_the_anchor",
            })

        # ---- OFF_REGION: seed projections that landed elsewhere ------------
        off_region = [b for b in seed_beds if b.name not in in_region_names]
        off_region_names = {b.name for b in off_region}

        # ---- FAMILY: projections filed under a paralog's name --------------
        family_genes = family_map.get(gene_upper, set())
        family = [b for fg in family_genes for b in beds_by_toga_gene.get(fg, [])
                  if b.name not in in_region_names and b.name not in off_region_names]

        coord_labels = project_exon_labels(in_region, structure_by_gene[gene_upper])
        bat_novel: Dict[Tuple[int, int], str] = {}
        home_up, home_down = neighbouring_genes(gene_bed, home_chrom, home_lo,
                                                home_hi, home_strand, ensg_map)

        for pool, beds in (("IN_REGION", in_region), ("OFF_REGION", off_region),
                           ("FAMILY", family)):
            for b in sorted(beds, key=lambda x: (x.chrom, x.start, x.name)):
                proj = parse_projection(b.name)
                exs_order = b.blocks_transcript_order
                labels = []
                for block in exs_order:
                    lab = coord_labels.get(block)
                    if lab is None:
                        lab = bat_novel.get(block)
                        if lab is None:
                            lab = f"bat_novel_exon{len(bat_novel) + 1}"
                            bat_novel[block] = lab
                    labels.append(lab)
                matched = sorted({int(m.group(1)) for m in
                                  (re.fullmatch(r"exon(\d+)", x) for x in labels) if m})

                if b.chrom != home_chrom:
                    distance = -1                      # different scaffold
                elif b.end <= home_lo:
                    distance = home_lo - b.end
                elif b.start >= home_hi:
                    distance = b.start - home_hi
                else:
                    distance = 0

                source = source_of.get((b.name, b.chrom, b.start, b.end), "NA")
                candidate_id = sanitize(
                    f"{species}|{gene}|{b.chrom}|{b.start}|{proj.raw}", "cand")
                up, down = ((home_up, home_down) if pool == "IN_REGION"
                            else neighbouring_genes(gene_bed, b.chrom, b.start,
                                                    b.end, b.strand, ensg_map))
                candidate_rows.append({
                    "species": species, "gene": gene, "candidate_id": candidate_id,
                    "projection": proj.raw, "toga_transcript": proj.transcript,
                    "toga_gene_label": proj.gene, "chain": proj.chain,
                    "is_retro": int(proj.is_retro or source == "processed_pseudogene"),
                    "annotation_source": source,
                    "toga_status": status_by_proj.get(proj.raw,
                                                      status_by_proj.get(proj.short, "NA")),
                    "orthology_class": orth_by_proj.get(proj.raw,
                                                        orth_by_proj.get(proj.short, "NA")),
                    "chrom": b.chrom, "start": b.start, "end": b.end,
                    "strand": b.strand, "exon_count": b.block_count,
                    "cds_bp": b.cds_bp, "seq_length": 0,
                    "locus_id": f"{species}|{gene}|{b.chrom}:{b.start}-{b.end}",
                    "pool": pool, "distance_to_home_locus_bp": distance,
                    "exon_structure": ",".join(labels),
                    "matched_longest_exons": ",".join(map(str, matched)) or "NONE",
                    "n_matched_longest_exons": len(matched),
                    "n_bat_novel_exons": sum(1 for x in labels
                                             if x.startswith("bat_novel_exon")),
                    "upstream_gene": up, "downstream_gene": down,
                })

                is_anchor = int(pool == "IN_REGION" and b.name == home_anchor.name)
                for rank, ((ex_s, ex_e), lab) in enumerate(zip(exs_order, labels), 1):
                    exon_rows.append({
                        "species": species, "gene": gene,
                        "candidate_id": candidate_id, "projection": proj.raw,
                        "toga_transcript": proj.transcript,
                        "toga_gene_label": proj.gene, "pool": pool,
                        "chrom": b.chrom, "strand": b.strand,
                        "start": ex_s, "end": ex_e,
                        "length_bp": ex_e - ex_s,
                        "transcript_exon_rank": rank,
                        "exon_label": lab,
                        "exon_label_type": ("longest" if re.fullmatch(r"exon\d+", lab)
                                            else "human_novel"
                                            if lab.startswith("human_novel_exon")
                                            else "bat_novel"),
                        "is_anchor": is_anchor,
                        "previous_gene": up, "next_gene": down,
                        "gene_order": f"{up}->{gene}->{down}",
                    })
                wanted_seq_keys[proj.short].append(candidate_id)

    # ---- one pass over the (large) nucleotide FASTA ------------------------
    log(f"{species}: extracting {len(candidate_rows)} candidate CDS from {fa_path.name}")
    seq_by_candidate: Dict[str, str] = {}
    for header, seq in fasta_iter(fa_path):
        for cid in wanted_seq_keys.get(parse_projection(header).short, ()):
            seq_by_candidate[cid] = clean_nt(seq)

    records = []
    for row in candidate_rows:
        seq = seq_by_candidate.get(row["candidate_id"], "")
        row["seq_length"] = len(seq)
        if seq:
            records.append((row["candidate_id"], seq))

    write_tsv(outdir / f"{species}.candidates.tsv", candidate_rows, CANDIDATE_FIELDS)
    write_tsv(outdir / f"{species}.loci.tsv", locus_rows, LOCUS_FIELDS)
    write_tsv(outdir / f"{species}.excluded_nested.tsv", nested_rows, NESTED_FIELDS)
    write_tsv(outdir / f"{species}.candidate_exons.tsv", exon_rows, EXON_FIELDS)
    write_fasta(outdir / f"{species}.candidate_cds.fa", records)

    n_missing = len(candidate_rows) - len(records)
    log(f"{species}: {len(candidate_rows)} candidates across "
        f"{len(ref_ctx)} genes ({n_missing} without CDS sequence)")
    return outdir / f"{species}.candidates.tsv"


def merge(per_species_dir, outdir) -> None:
    """Concatenate the per-species tables and FASTAs into combined files."""
    per_species_dir, outdir = Path(per_species_dir), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for suffix, fields in (("candidates.tsv", CANDIDATE_FIELDS),
                           ("loci.tsv", LOCUS_FIELDS),
                           ("excluded_nested.tsv", NESTED_FIELDS),
                           ("candidate_exons.tsv", EXON_FIELDS)):
        rows = []
        for p in sorted(per_species_dir.glob(f"*.{suffix}")):
            rows.extend(read_tsv(p))
        write_tsv(outdir / f"all_species_{suffix}", rows, fields)
        log(f"merged {len(rows)} rows -> all_species_{suffix}")

    records = []
    for p in sorted(per_species_dir.glob("*.candidate_cds.fa")):
        records.extend(fasta_iter(p))
    write_fasta(outdir / "all_species_candidate_cds.fa", records)
    log(f"merged {len(records)} CDS -> all_species_candidate_cds.fa")
