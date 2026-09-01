"""Recombine per-gene runs into the whole-list tables and figures.

When each gene runs on its own -- a separate Slurm job, or one at a time
locally -- each writes a complete, self-contained results/ tree under its own
directory. That isolation is the point: one gene failing cannot corrupt a shared
table, and nothing has to be locked.

The cost is that the figures people actually want are multi-gene: a status
matrix with one column per human transcript across every gene, and a copy-number
matrix with one column per gene. This step concatenates the per-gene tables back
into whole-list ones so those figures can be drawn.

Concatenation is by column name, not position: a per-gene table written by an
older version with a different column order still merges correctly, and a
missing column becomes NA rather than silently shifting every value one place
to the left.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence

from .common import log, read_tsv, write_tsv

# Per-gene path -> merged path, relative to a gene's results/ and to the
# combined directory. {policy} and {scope} are expanded over what exists.
MERGE_TARGETS = [
    ("04_copy_number/copy_number_matrix.tsv", "copy_number_matrix.tsv"),
    ("04_copy_number/gene_copies.tsv",        "gene_copies.tsv"),
    ("02_bat_search/paralog_screen.tsv",      "paralog_screen.tsv"),
    ("01_human_reference/gene_context.tsv",   "gene_context.tsv"),
]

POLICY_TARGETS = [
    ("05_genes/{policy}/manifest.tsv",                    "manifest__{policy}.tsv"),
    ("05_genes/{policy}/transcript_status_{scope}.tsv",   "transcript_status_{scope}__{policy}.tsv"),
]

ASSIGN_TARGETS = [
    ("03_transcript_assignment/transcript_assignments_{scope}__{policy}.tsv",
     "transcript_assignments_{scope}__{policy}.tsv"),
    ("03_transcript_assignment/scope_comparison__{policy}.tsv",
     "scope_comparison__{policy}.tsv"),
]


def _concat(paths: Sequence[Path], out: Path) -> Optional[Path]:
    """Concatenate headed TSVs by column name. Returns None if none existed."""
    rows: List[dict] = []
    fields: List[str] = []
    seen = set()
    for p in paths:
        for r in read_tsv(p):
            for k in r:
                if k not in seen:
                    seen.add(k)
                    fields.append(k)
            rows.append(r)
    if not fields:
        return None
    # A gene whose table lacked a column gets NA there, rather than the next
    # column's value sliding into its place.
    for r in rows:
        for k in fields:
            r.setdefault(k, "NA")
    out.parent.mkdir(parents=True, exist_ok=True)
    return write_tsv(out, rows, fields)


def _gene_dirs(genes_root: Path) -> List[Path]:
    return sorted(d for d in genes_root.iterdir()
                  if d.is_dir() and (d / "results").is_dir())


def merge(genes_root, outdir, policies: Sequence[str],
          scopes: Sequence[str]) -> Path:
    genes_root = Path(genes_root)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    dirs = _gene_dirs(genes_root)
    if not dirs:
        raise SystemExit(f"no per-gene results under {genes_root}")

    written: Dict[str, int] = {}

    def do(rel: str, name: str) -> None:
        paths = [d / "results" / rel for d in dirs]
        present = [p for p in paths if p.is_file()]
        got = _concat(present, outdir / name)
        if got is not None:
            written[name] = len(present)

    for rel, name in MERGE_TARGETS:
        do(rel, name)
    for policy in policies:
        for rel, name in POLICY_TARGETS:
            if "{scope}" in rel:
                for scope in scopes:
                    do(rel.format(policy=policy, scope=scope),
                       name.format(policy=policy, scope=scope))
            else:
                do(rel.format(policy=policy), name.format(policy=policy))
        for rel, name in ASSIGN_TARGETS:
            if "{scope}" in rel:
                for scope in scopes:
                    do(rel.format(policy=policy, scope=scope),
                       name.format(policy=policy, scope=scope))
            else:
                do(rel.format(policy=policy), name.format(policy=policy))

    missing = [d.name for d in dirs
               if not (d / "results" / "04_copy_number" / "copy_number_matrix.tsv").is_file()]
    if missing:
        log(f"WARNING: {len(missing)} gene(s) produced no copy-number table and "
            f"may have failed: {', '.join(sorted(missing))}")

    log(f"merged {len(dirs)} gene directories into {outdir}")
    for name in sorted(written):
        log(f"  {name}  <- {written[name]} genes")
    return outdir
