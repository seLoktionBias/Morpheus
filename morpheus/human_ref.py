"""Step 1 - build the human reference for the genes of interest.

For every gene in the gene list this produces:

  * the genomic region of its longest protein-coding CDS ("home locus");
  * every human protein-coding transcript whose CDS lies in that region and
    shares at least one CDS exon with the longest isoform (region-restricted,
    so a neighbouring paralog does not leak in);
  * an exon-label vocabulary (exon1..exonN from the longest isoform, plus
    human_novel_exonN for isoform-specific exons) used later to project the
    human exon structure into each bat genome;
  * per-transcript CDS and protein FASTAs;
  * a per-exon table (human_isoform_exons.tsv) for the optional gene-model
    figure, so any gene can be drawn on demand.

It also writes a genome-wide "longest protein per gene" FASTA, which the copy
number step uses to tell a real gene copy apart from a TOGA2-misannotated
paralog.
"""
from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .common import (die, fasta_iter, log, opn, read_lines, sanitize,
                     strip_version, write_fasta, write_tsv)

CDS_EXON_FIELDS = ["chrom", "start", "end", "strand", "gene_id", "gene",
                   "transcript_id", "exon_number"]


# --------------------------------------------------------------------------
# Ensembl GTF -> flat CDS exon table (cached, built once)
# --------------------------------------------------------------------------

_ATTR = re.compile(r'(\S+) "([^"]*)"')


def build_cds_exon_table(gtf_path, out_path) -> Path:
    """Flatten the CDS features of an Ensembl GTF into a sorted TSV.

    Only protein_coding transcripts are kept. Coordinates stay 1-based
    inclusive, exactly as in the GTF.
    """
    out_path = Path(out_path)
    if out_path.exists() and out_path.stat().st_size > 0:
        log(f"reusing cached CDS exon table {out_path.name}")
        return out_path

    log(f"parsing CDS features from {Path(gtf_path).name} (this takes a few minutes)")
    rows = []
    kept = 0
    with opn(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            attrs = dict(_ATTR.findall(f[8]))
            if attrs.get("transcript_biotype") != "protein_coding":
                continue
            rows.append((f[0], int(f[3]), int(f[4]), f[6],
                         attrs.get("gene_id", "NA"),
                         attrs.get("gene_name", attrs.get("gene_id", "NA")),
                         strip_version(attrs.get("transcript_id", "NA")),
                         int(attrs.get("exon_number", 0) or 0)))
            kept += 1
            if kept % 500000 == 0:
                log(f"  {kept:,} CDS exons parsed")

    rows.sort(key=lambda r: (r[0], r[1], r[2], r[6], r[7]))
    # Keep the .gz suffix on the temporary name so opn() still compresses it.
    tmp = out_path.with_name("." + out_path.name + ".partial" + out_path.suffix)
    tmp.parent.mkdir(parents=True, exist_ok=True)
    with opn(tmp, "wt") as oh:
        oh.write("\t".join(CDS_EXON_FIELDS) + "\n")
        for r in rows:
            oh.write("\t".join(str(x) for x in r) + "\n")
    tmp.replace(out_path)
    log(f"wrote {kept:,} CDS exons for "
        f"{len({r[6] for r in rows}):,} transcripts -> {out_path.name}")
    return out_path


def load_cds_exon_table(path) -> Dict[str, dict]:
    """Load the flat table into {transcript_id: {...,'exons':[...] }}."""
    tx: Dict[str, dict] = {}
    with opn(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            tid = f[idx["transcript_id"]]
            rec = tx.get(tid)
            if rec is None:
                rec = tx[tid] = {"transcript_id": tid,
                                 "chrom": f[idx["chrom"]],
                                 "strand": f[idx["strand"]],
                                 "gene_id": f[idx["gene_id"]],
                                 "gene": f[idx["gene"]],
                                 "exons": []}
            rec["exons"].append({"start": int(f[idx["start"]]),
                                 "end": int(f[idx["end"]]),
                                 "exon_number": int(f[idx["exon_number"]])})
    for rec in tx.values():
        rec["exons"].sort(key=lambda e: e["exon_number"] or e["start"])
    return tx


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------


def transcript_region(rec: dict) -> Tuple[str, int, int, str]:
    starts = [e["start"] for e in rec["exons"]]
    ends = [e["end"] for e in rec["exons"]]
    return rec["chrom"], min(starts), max(ends), rec["strand"]


def _exons_overlap(a: List[dict], b: List[dict]) -> bool:
    for x in a:
        for y in b:
            if x["end"] >= y["start"] and x["start"] <= y["end"]:
                return True
    return False


def parse_longest_pc_tx(path) -> Dict[str, dict]:
    """longest_pc_tx.tsv: gene_id, transcript_id, cds_length, gene_symbol."""
    by_symbol: Dict[str, dict] = {}
    for line in read_lines(path):
        f = line.split("\t")
        if len(f) < 4:
            continue
        if f[0].lower().startswith("gene") and f[1].lower().startswith("trans"):
            continue  # header
        rec = {"gene_id": f[0].strip(), "transcript_id": strip_version(f[1]),
               "cds_length": f[2].strip(), "gene": f[3].strip()}
        by_symbol[rec["gene"].upper()] = rec
    return by_symbol


# --------------------------------------------------------------------------
# main entry point
# --------------------------------------------------------------------------


def build(gene_list_path, cds_exon_table, longest_pc_tx, human_cds_fasta,
          outdir) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genes = read_lines(gene_list_path)
    longest_by_symbol = parse_longest_pc_tx(longest_pc_tx)
    log(f"loading CDS exon table")
    all_tx = load_cds_exon_table(cds_exon_table)

    log("indexing human CDS FASTA")
    cds_by_tx: Dict[str, str] = {}
    for header, seq in fasta_iter(human_cds_fasta):
        cds_by_tx[strip_version(header.split()[0])] = seq

    # transcripts grouped by chromosome, so the region scan stays cheap
    by_chrom: Dict[str, List[dict]] = defaultdict(list)
    for rec in all_tx.values():
        by_chrom[rec["chrom"]].append(rec)

    # Genome-wide CDS span per gene, so each target gene can report the
    # protein-coding genes flanking it. Synteny is how a projected locus gets
    # checked, so the neighbours belong in the reference, not just in the
    # per-genome tables.
    spans: Dict[str, Dict[str, List[int]]] = defaultdict(dict)
    for rec in all_tx.values():
        chrom, lo, hi, _ = transcript_region(rec)
        cur = spans[chrom].get(rec["gene"])
        if cur is None:
            spans[chrom][rec["gene"]] = [lo, hi]
        else:
            cur[0] = min(cur[0], lo)
            cur[1] = max(cur[1], hi)
    gene_order_by_chrom = {c: sorted(g.items(), key=lambda kv: kv[1])
                           for c, g in spans.items()}

    context_rows, transcript_rows, label_rows = [], [], []
    exon_rows: List[dict] = []
    cds_records, protein_records = [], []

    for gene in genes:
        info = longest_by_symbol.get(gene.upper())
        if not info:
            context_rows.append({"gene": gene, "status": "NOT_IN_LONGEST_PC_TX"})
            log(f"  {gene}: not present in longest_pc_tx.tsv - skipped")
            continue
        longest_id = info["transcript_id"]
        longest = all_tx.get(longest_id)
        if not longest:
            context_rows.append({"gene": gene, "gene_id": info["gene_id"],
                                 "longest_transcript": longest_id,
                                 "status": "LONGEST_TX_HAS_NO_CDS"})
            log(f"  {gene}: longest transcript {longest_id} absent from GTF CDS - skipped")
            continue

        chrom, lo, hi, strand = transcript_region(longest)

        # Region-restricted transcript set: same chromosome and strand, span
        # overlapping the longest isoform, and sharing >=1 CDS exon with it.
        members = []
        for rec in by_chrom[chrom]:
            if rec["strand"] != strand:
                continue
            _, s, e, _ = transcript_region(rec)
            if e < lo or s > hi:
                continue
            if rec["transcript_id"] == longest_id or _exons_overlap(rec["exons"], longest["exons"]):
                members.append(rec)
        if longest_id not in {r["transcript_id"] for r in members}:
            members.append(longest)

        # ---- exon label vocabulary -------------------------------------
        longest_label = {}
        for i, ex in enumerate(longest["exons"], 1):
            longest_label[(ex["start"], ex["end"])] = f"exon{i}"
            label_rows.append({"gene": gene, "gene_id": info["gene_id"],
                               "exon_label": f"exon{i}", "chrom": chrom,
                               "start": ex["start"], "end": ex["end"],
                               "strand": strand, "label_type": "longest_isoform"})

        novel_label: Dict[Tuple[int, int], str] = {}
        structures: Dict[str, List[str]] = {}
        for rec in sorted(members, key=lambda r: r["transcript_id"]):
            labels = []
            for ex in rec["exons"]:
                key = (ex["start"], ex["end"])
                lab = longest_label.get(key)
                if lab is None:
                    lab = novel_label.get(key)
                    if lab is None:
                        lab = f"human_novel_exon{len(novel_label) + 1}"
                        novel_label[key] = lab
                        label_rows.append({
                            "gene": gene, "gene_id": info["gene_id"],
                            "exon_label": lab, "chrom": chrom,
                            "start": ex["start"], "end": ex["end"],
                            "strand": strand,
                            "label_type": "isoform_specific"})
                labels.append(lab)
            structures[rec["transcript_id"]] = labels

        # ---- collapse transcripts with an identical CDS ----------------
        # The key is the CDS sequence alone. Everything downstream - BLAST,
        # codon alignment - is sequence-based, so two transcripts with
        # the same CDS are the same analysis. Keying on exon structure as well
        # would keep such a pair apart, and each copy would then consume a
        # separate query transcript in the one-to-one matching and be aligned
        # and tested twice.
        #
        # Transcripts can reach an identical CDS from different exon
        # coordinates (alternative first exons that encode the same bases), so
        # whether the underlying model is also identical is recorded separately;
        # the exon-label projection must not borrow a structure that does not
        # actually match.
        chain_of = {rec["transcript_id"]: tuple((e["start"], e["end"])
                                                for e in rec["exons"])
                    for rec in members}

        groups: Dict[str, List[str]] = defaultdict(list)
        for rec in members:
            tid = rec["transcript_id"]
            groups[cds_by_tx.get(tid, "")].append(tid)

        representatives: Dict[str, List[str]] = {}
        same_model: Dict[str, List[str]] = {}
        def canonical_exon_count(tid: str) -> int:
            return sum(1 for lab in structures[tid] if re.fullmatch(r"exon\d+", lab))

        for _, tids in groups.items():
            tids = sorted(tids)
            if longest_id in tids:
                rep = longest_id
            else:
                # Among equal-CDS transcripts, keep the one whose exon chain is
                # expressed most in canonical exonN labels. That chain is what
                # gets projected into every query genome, so a representative
                # carrying human_novel labels for exons the longest isoform
                # also has would weaken the structure comparison downstream.
                rep = max(tids, key=lambda t: (canonical_exon_count(t), -len(t), t))
            representatives[rep] = tids
            same_model[rep] = [t for t in tids
                               if chain_of[t] == chain_of[rep]]

        n_written = 0
        for rec in sorted(members, key=lambda r: r["transcript_id"]):
            tid = rec["transcript_id"]
            if tid not in representatives:
                continue
            seq = cds_by_tx.get(tid, "")
            if not seq:
                log(f"  {gene}: no CDS sequence for {tid} - transcript dropped")
                continue
            exons = rec["exons"]
            transcript_rows.append({
                "gene": gene, "gene_id": info["gene_id"], "transcript_id": tid,
                "is_longest": int(tid == longest_id),
                "chrom": chrom, "strand": strand,
                "start": min(e["start"] for e in exons),
                "end": max(e["end"] for e in exons),
                "cds_exon_count": len(exons),
                "cds_bp": sum(e["end"] - e["start"] + 1 for e in exons),
                "cds_sequence_length": len(seq),
                "exon_structure": ",".join(structures[tid]),
                # same CDS *and* same exon coordinates: safe to share this
                # transcript's exon-label chain when projecting into a genome
                "identical_model_transcripts":
                    ",".join(t for t in same_model[tid] if t != tid) or "NONE",
                # same CDS, different exon model: still redundant for sequence
                # analysis, but its structure is genuinely different
                "identical_cds_other_model_transcripts":
                    ",".join(t for t in representatives[tid]
                             if t != tid and t not in same_model[tid]) or "NONE",
            })
            cds_records.append((tid, seq))
            protein_records.append((tid, _protein(seq)))
            for rank, ex in enumerate(exons, 1):
                exon_rows.append({
                    "gene": gene, "gene_id": info["gene_id"],
                    "transcript_id": tid, "is_longest": int(tid == longest_id),
                    "chrom": chrom, "strand": strand,
                    "start": ex["start"], "end": ex["end"],
                    "length_bp": ex["end"] - ex["start"] + 1,
                    "transcript_exon_rank": rank,
                    "gtf_exon_number": ex["exon_number"],
                    "exon_label": structures[tid][rank - 1],
                    "exon_label_type": ("longest_isoform"
                                        if re.fullmatch(r"exon\d+", structures[tid][rank - 1])
                                        else "isoform_specific"),
                })
            n_written += 1

        ordered = gene_order_by_chrom.get(chrom, [])
        idx = next((i for i, (name, _) in enumerate(ordered) if name == gene), None)
        prev_gene = ordered[idx - 1][0] if idx not in (None, 0) else "NA"
        next_gene = (ordered[idx + 1][0]
                     if idx is not None and idx + 1 < len(ordered) else "NA")
        # report in transcription order, not chromosome order
        upstream, downstream = ((next_gene, prev_gene) if strand == "-"
                                else (prev_gene, next_gene))

        context_rows.append({
            "gene": gene, "gene_id": info["gene_id"],
            "previous_gene": prev_gene, "next_gene": next_gene,
            "upstream_gene": upstream, "downstream_gene": downstream,
            "gene_order": f"{upstream}->{gene}->{downstream}",
            "longest_transcript": longest_id, "chrom": chrom,
            "start": lo, "end": hi, "strand": strand,
            "longest_cds_exon_count": len(longest["exons"]),
            "longest_cds_bp": sum(e["end"] - e["start"] + 1 for e in longest["exons"]),
            "n_transcripts": n_written,
            "n_transcripts_before_dedup": len(members),
            "status": "OK",
        })
        log(f"  {gene}: {n_written} representative transcript(s) "
            f"from {len(members)} in region {chrom}:{lo}-{hi}({strand})")

    write_tsv(outdir / "gene_context.tsv", context_rows,
              ["gene", "gene_id", "longest_transcript", "chrom", "start", "end",
               "strand", "longest_cds_exon_count", "longest_cds_bp",
               "previous_gene", "next_gene", "upstream_gene", "downstream_gene",
               "gene_order", "n_transcripts", "n_transcripts_before_dedup",
               "status"])
    write_tsv(outdir / "transcripts.tsv", transcript_rows,
              ["gene", "gene_id", "transcript_id", "is_longest", "chrom", "strand",
               "start", "end", "cds_exon_count", "cds_bp", "cds_sequence_length",
               "exon_structure", "identical_model_transcripts",
               "identical_cds_other_model_transcripts"])
    write_tsv(outdir / "human_isoform_exons.tsv", exon_rows,
              ["gene", "gene_id", "transcript_id", "is_longest", "chrom",
               "strand", "start", "end", "length_bp", "transcript_exon_rank",
               "gtf_exon_number", "exon_label", "exon_label_type"])
    write_tsv(outdir / "exon_label_coordinates.tsv", label_rows,
              ["gene", "gene_id", "exon_label", "chrom", "start", "end",
               "strand", "label_type"])
    write_fasta(outdir / "human_transcript_cds.fa", cds_records)
    write_fasta(outdir / "human_transcript_proteins.fa", protein_records)

    # TOGA2 names query gene regions by Ensembl gene ID, so the neighbouring
    # genes in every locus table come out as ENSG accessions unless they can be
    # translated. Emit the map once, genome-wide.
    ensg = {rec["gene_id"]: rec["gene"] for rec in all_tx.values()
            if rec.get("gene_id") and rec.get("gene")}
    write_tsv(outdir / "ensembl_gene_id_to_symbol.tsv",
              [{"gene_id": k, "gene": v} for k, v in sorted(ensg.items())],
              ["gene_id", "gene"])
    log(f"ENSG -> symbol map: {len(ensg):,} genes")

    _write_genomewide_protein_db(longest_by_symbol, cds_by_tx,
                                 outdir / "human_longest_proteins_all_genes.fa")

    log(f"human reference written to {outdir}")


def _protein(cds: str) -> str:
    from .common import translate
    return translate(cds).rstrip("*")


def _write_genomewide_protein_db(longest_by_symbol, cds_by_tx, path) -> None:
    """One protein per human gene, used to catch TOGA2 paralog misannotation."""
    records = []
    for symbol, info in sorted(longest_by_symbol.items()):
        seq = cds_by_tx.get(info["transcript_id"])
        if not seq:
            continue
        records.append((f"{info['gene']}|{info['transcript_id']}", _protein(seq)))
    write_fasta(path, records)
    log(f"genome-wide paralog-screening protein DB: {len(records):,} genes")
