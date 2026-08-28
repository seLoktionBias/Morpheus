"""Define each target gene's paralog family from the human proteome.

Copy number cannot be counted honestly without knowing which other human genes
could plausibly be confused with the target. Rather than hard-coding families
(IFITM1/2/3/5/10, OAS1/2/3/OASL, MX1/MX2 ...), the family is derived: the
target's longest human protein is BLASTed against one protein per human gene,
and every gene scoring within a bitscore fraction of the self-hit joins the
family.

Two things use the result:

  * the bat search widens its net to projections TOGA2 labelled with *any*
    family member, so a real gene copy misannotated as its paralog is still
    found;
  * the copy-number step then keeps only those whose best genome-wide match is
    the target itself, which is what delimits a true copy from the paralog.
"""
from __future__ import annotations

import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set

from .common import die, fasta_iter, log, read_lines, read_tsv, to_float, write_tsv

FAMILY_FIELDS = ["gene", "family_gene", "bitscore", "pident",
                 "bitscore_fraction_of_self", "is_self"]


def build(reference_dir, outdir, min_bitscore_fraction: float = 0.25,
          threads: int = 4) -> Path:
    """BLAST each target's human protein against all human genes."""
    if shutil.which("blastp") is None or shutil.which("makeblastdb") is None:
        die("blastp/makeblastdb not on PATH. Activate the analysis environment.")

    reference_dir, outdir = Path(reference_dir), Path(outdir)
    work = outdir / "family_work"
    work.mkdir(parents=True, exist_ok=True)

    db_fa = reference_dir / "human_longest_proteins_all_genes.fa"
    all_prot = {h.split()[0]: s for h, s in fasta_iter(db_fa)}

    # One query protein per target gene: its longest isoform.
    ctx = [r for r in read_tsv(reference_dir / "gene_context.tsv")
           if r.get("status") == "OK"]
    targets = {r["gene"]: r["longest_transcript"] for r in ctx}

    query_fa = work / "target_proteins.fa"
    records = []
    for gene, tx in sorted(targets.items()):
        key = next((k for k in all_prot if k.split("|")[0].upper() == gene.upper()), None)
        if key is None:
            log(f"  {gene}: absent from the genome-wide protein DB - family = itself only")
            continue
        records.append((f"{gene}|{tx}", all_prot[key]))
    if not records:
        die("no target proteins available for family definition")
    from .common import write_fasta
    write_fasta(query_fa, records)

    db = str(work / "human_all_genes")
    if not Path(db + ".pin").exists() and not Path(db + ".pdb").exists():
        subprocess.run(["makeblastdb", "-in", str(db_fa), "-dbtype", "prot",
                        "-out", db], check=True, stdout=subprocess.DEVNULL)

    hits = work / "target_vs_all_human_genes.tsv"
    subprocess.run(
        ["blastp", "-query", str(query_fa), "-db", db, "-out", str(hits),
         "-evalue", "1e-5", "-seg", "no", "-max_target_seqs", "500",
         "-num_threads", str(threads),
         "-outfmt", "6 qseqid sseqid pident bitscore"],
        check=True)

    best: Dict[str, Dict[str, tuple]] = defaultdict(dict)
    with open(hits) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            q = f[0].split("|")[0].upper()
            s = f[1].split("|")[0].upper()
            bits, pid = to_float(f[3]), to_float(f[2])
            cur = best[q].get(s)
            if cur is None or bits > cur[0]:
                best[q][s] = (bits, pid)

    rows = []
    for gene in sorted(targets):
        g = gene.upper()
        self_bits = best[g].get(g, (0.0, 0.0))[0]
        if self_bits <= 0:
            self_bits = max((b for b, _ in best[g].values()), default=1.0)
        n = 0
        for other, (bits, pid) in sorted(best[g].items(), key=lambda kv: -kv[1][0]):
            frac = bits / self_bits if self_bits else 0.0
            if frac < min_bitscore_fraction and other != g:
                continue
            rows.append({"gene": gene, "family_gene": other,
                         "bitscore": round(bits, 1), "pident": round(pid, 2),
                         "bitscore_fraction_of_self": round(frac, 4),
                         "is_self": int(other == g)})
            n += 1
        log(f"  {gene}: paralog family of {n} human gene(s)")

    return write_tsv(outdir / "gene_families.tsv", rows, FAMILY_FIELDS)


def load(path) -> Dict[str, Set[str]]:
    """{TARGET_GENE_UPPER: {family gene symbols, upper-case}}

    The relation is symmetrised on load. BLAST membership is directional - OAS1
    reaches OAS3 above the threshold while OAS3 does not reach OAS1 - and left
    as it comes, "is this a sibling?" gets different answers depending on which
    gene is asked. That inconsistency reads as a real paralog to the screen, so
    OAS3's own locus in some species was thrown away for looking like OAS2.
    Paralogy is a symmetric relation; make it one before anyone consumes it.
    """
    raw: Dict[str, Set[str]] = defaultdict(set)
    for r in read_tsv(path):
        raw[r["gene"].upper()].add(r["family_gene"].upper())

    fam: Dict[str, Set[str]] = defaultdict(set)
    for gene, members in raw.items():
        fam[gene] |= members | {gene}
        for other in members:
            # only reciprocate between genes that are themselves targets
            if other in raw:
                fam[other].add(gene)
    return fam
