"""Step 2.5 - genome-wide paralog screen.

TOGA2 attaches a human gene name to every projection, and that name is not
always right: in a family with many paralogs a projection of one gene can be
labelled with a transcript of its neighbour. Neither the transcript assignment
nor the copy number count should take the label at face value.

So every bat candidate CDS is BLASTed once against a database holding one
protein per human protein-coding gene (the longest isoform, genome-wide). The
result gives, per candidate:

  * which human gene it actually looks most like;
  * how much better that gene scores than the gene we are targeting.

A candidate that sits in the target gene's home locus but looks far more like a
paralog is exactly the misannotation the project needs to catch, and a
candidate that looks like the target gene but sits elsewhere is exactly the
off-region copy that must not be silently discarded.

One thing sequence similarity cannot do is separate the members of a tight
tandem array. Human IFITM1, IFITM2 and IFITM3 sit within ~12 kb of each other
and share most of their coding sequence; every bat IFITM scores highest against
IFITM3 whichever locus it came from. Treating "best hit is a different gene" as
disqualifying there deletes two of the three family members. So the verdict
depends on whether the better-scoring gene is *inside the target's own paralog
family*:

  * outside the family - a genuine misannotation, and disqualifying;
  * inside the family, in the home locus - synteny decides, and the candidate is
    kept and flagged rather than dropped;
  * inside the family, outside the locus - most simply read as the sibling's own
    copy, so it is recorded but does not become a transcript or a copy here.
"""
from __future__ import annotations

import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

from .common import die, fasta_iter, log, read_tsv, to_float, write_tsv


def _digest(path) -> str:
    """Content hash of a file, read in chunks so large FASTAs stay cheap."""
    import hashlib
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

SCREEN_FIELDS = ["candidate_id", "species", "gene", "pool",
                 "best_human_gene", "best_bitscore", "best_pident",
                 "target_gene_bitscore", "target_gene_pident",
                 "bitscore_margin_over_target", "looks_like_target_gene",
                 "best_gene_is_in_target_family",
                 "second_best_human_gene", "second_best_bitscore",
                 "screen_verdict"]

# Verdicts that keep a candidate in play for the target gene.
ACCEPTING_VERDICTS = {"CONSISTENT_WITH_TARGET_GENE",
                      "CONSISTENT_BY_SYNTENY_WITHIN_FAMILY"}


def run(candidate_fasta, candidates_tsv, reference_dir, outdir,
        threads: int = 4, evalue: str = "1e-5", top_n: int = 25,
        families_tsv=None) -> Path:
    """BLAST all candidates against the one-protein-per-human-gene database."""
    if shutil.which("blastx") is None or shutil.which("makeblastdb") is None:
        die("blastx/makeblastdb not on PATH. Activate the analysis environment.")

    outdir = Path(outdir)
    work = outdir / "screen_work"
    work.mkdir(parents=True, exist_ok=True)

    db_fa = Path(reference_dir) / "human_longest_proteins_all_genes.fa"
    if not db_fa.exists():
        die(f"missing paralog-screening database: {db_fa}")

    from . import family as family_mod
    fam_path = families_tsv or (Path(reference_dir) / "gene_families.tsv")
    families = family_mod.load(fam_path) if Path(fam_path).exists() else {}
    if not families:
        log("WARNING: no gene_families.tsv - every better-scoring gene will be "
            "treated as out-of-family, which is too strict for tandem arrays")

    db = str(work / "human_all_genes")
    if not Path(db + ".pin").exists() and not Path(db + ".pdb").exists():
        log("building genome-wide human protein BLAST database")
        subprocess.run(["makeblastdb", "-in", str(db_fa), "-dbtype", "prot",
                        "-out", db], check=True, stdout=subprocess.DEVNULL)

    # The BLAST is the expensive part and depends only on the query FASTA and
    # the database, not on how verdicts are decided. Reuse it when the query has
    # not changed since the hits were written.
    hits_path = work / "candidates_vs_all_human_genes.tsv"
    # Compare content, not timestamps: re-running the search rewrites the FASTA
    # with an identical body, and a newer mtime is not a reason to redo hours of
    # BLAST.
    stamp_path = work / "candidates_vs_all_human_genes.query"
    stamp = _digest(candidate_fasta) + " " + _digest(db_fa)
    fresh = (hits_path.exists() and hits_path.stat().st_size > 0
             and stamp_path.exists() and stamp_path.read_text().strip() == stamp)
    if fresh:
        log(f"reusing existing BLAST hits ({hits_path.name}); "
            f"delete it to force a rescan")
    else:
        log("screening every candidate against all human genes")
        subprocess.run(
            ["blastx", "-query", str(candidate_fasta), "-db", db,
             "-out", str(hits_path), "-evalue", evalue, "-seg", "no",
             "-max_target_seqs", str(top_n), "-num_threads", str(threads),
             "-outfmt", "6 qseqid sseqid pident length bitscore evalue"],
            check=True)
        stamp_path.write_text(stamp + "\n")

    # best bitscore per (candidate, human gene)
    best: Dict[str, Dict[str, Tuple[float, float]]] = defaultdict(dict)
    with open(hits_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            gene = f[1].split("|")[0].upper()
            bits, pid = to_float(f[4]), to_float(f[2])
            cur = best[f[0]].get(gene)
            if cur is None or bits > cur[0]:
                best[f[0]][gene] = (bits, pid)

    rows = []
    for c in read_tsv(candidates_tsv, required=["candidate_id", "gene"]):
        cid, target = c["candidate_id"], c["gene"].upper()
        per_gene = best.get(cid, {})
        ranked = sorted(per_gene.items(), key=lambda kv: -kv[1][0])
        top_gene, top_bits, top_pid = ("NA", 0.0, 0.0)
        if ranked:
            top_gene, (top_bits, top_pid) = ranked[0][0], ranked[0][1]
        second_gene, second_bits = ("NA", 0.0)
        if len(ranked) > 1:
            second_gene, second_bits = ranked[1][0], ranked[1][1][0]
        tgt_bits, tgt_pid = per_gene.get(target, (0.0, 0.0))

        looks_like_target = int(top_gene == target)
        in_family = int(top_gene in families.get(target, set()))
        if not ranked:
            verdict = "NO_HIT"
        elif looks_like_target:
            verdict = "CONSISTENT_WITH_TARGET_GENE"
        elif tgt_bits <= 0 and not in_family:
            verdict = "NO_SIMILARITY_TO_TARGET_GENE"
        elif not in_family:
            # a gene from outside the family scores better: real misannotation
            verdict = "LOOKS_LIKE_PARALOG"
        elif c.get("pool") == "IN_REGION":
            # sequence prefers a sibling, but the candidate is in this gene's
            # own locus; within a family synteny is the better evidence
            verdict = "CONSISTENT_BY_SYNTENY_WITHIN_FAMILY"
        else:
            # outside this gene's locus *and* it looks like a sibling: the
            # simplest reading is that it is the sibling's own copy. Recorded,
            # but it is not a transcript or a copy of the target.
            verdict = "BELONGS_TO_ANOTHER_FAMILY_MEMBER"

        rows.append({
            "candidate_id": cid, "species": c["species"], "gene": c["gene"],
            "pool": c.get("pool", "NA"),
            "best_human_gene": top_gene,
            "best_bitscore": round(top_bits, 1),
            "best_pident": round(top_pid, 2),
            "target_gene_bitscore": round(tgt_bits, 1),
            "target_gene_pident": round(tgt_pid, 2),
            "bitscore_margin_over_target": round(top_bits - tgt_bits, 1),
            "looks_like_target_gene": looks_like_target,
            "best_gene_is_in_target_family": in_family,
            "second_best_human_gene": second_gene,
            "second_best_bitscore": round(second_bits, 1),
            "screen_verdict": verdict,
        })

    out = write_tsv(outdir / "paralog_screen.tsv", rows, SCREEN_FIELDS)
    counts: Dict[str, int] = defaultdict(int)
    for r in rows:
        counts[r["screen_verdict"]] += 1
    log("paralog screen: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    return out


def load(path) -> Dict[str, dict]:
    return {r["candidate_id"]: r for r in read_tsv(path)}
