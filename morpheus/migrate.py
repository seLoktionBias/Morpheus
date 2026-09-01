"""Bring a results tree produced by an older version in line with current names.

v2.1.0 renamed the ranking policy `synteny_aware` to `structure_aware`. The
computation did not change, so an existing results tree is still completely
valid - but its directory names, file names and the paths recorded inside its
tables all carry the old word, and the current pipeline looks for the new one.

Renaming by hand is three `mv` commands plus a `sed`, and the `sed` is the part
people forget: manifest.tsv records absolute paths to every FASTA and tree it
wrote, and a rename that misses those leaves a manifest pointing at directories
that no longer exist.

Nothing here is destructive. A rename that would overwrite something is refused,
not forced, and the default is a dry run.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .common import log

# Renames this version knows how to perform, oldest first.
RENAMES: List[Tuple[str, str]] = [
    ("synteny_aware", "structure_aware"),
]

# Rewriting the inside of a file is only safe for text. Anything else is left
# alone; the figures were checked and carry no policy name of their own.
TEXT_SUFFIXES = {".tsv", ".csv", ".txt", ".md", ".log", ".nwk", ".fa", ".fasta",
                 ".json", ".yml", ".yaml", ".sh", ".py", ".R", ""}
MAX_REWRITE_BYTES = 200 * 1024 * 1024


def _is_text(path: Path) -> bool:
    if path.suffix not in TEXT_SUFFIXES:
        return False
    try:
        if path.stat().st_size > MAX_REWRITE_BYTES:
            return False
        with path.open("rb") as fh:
            return b"\0" not in fh.read(8192)
    except OSError:
        return False


def _plan(root: Path, old: str, new: str):
    """Directories to rename, files to rename, files to rewrite."""
    dirs, files, contents = [], [], []
    for p in root.rglob("*"):
        if old in p.name:
            (dirs if p.is_dir() else files).append(p)
    # Deepest first, so renaming a parent never invalidates a child's path.
    dirs.sort(key=lambda p: len(p.parts), reverse=True)
    files.sort(key=lambda p: len(p.parts), reverse=True)

    for p in root.rglob("*"):
        if not p.is_file() or not _is_text(p):
            continue
        try:
            if old in p.read_text(errors="ignore"):
                contents.append(p)
        except OSError:
            continue
    return dirs, files, sorted(contents)


def migrate(results_dir, apply: bool = False,
            renames: Optional[List[Tuple[str, str]]] = None) -> Dict[str, int]:
    root = Path(results_dir)
    if not root.is_dir():
        raise SystemExit(f"not a directory: {root}")
    renames = renames or RENAMES
    totals = {"directories": 0, "files": 0, "contents": 0, "skipped": 0}

    for old, new in renames:
        dirs, files, contents = _plan(root, old, new)
        if not (dirs or files or contents):
            log(f"nothing to migrate for '{old}' -> '{new}'")
            continue

        log(f"{'applying' if apply else 'would apply'} '{old}' -> '{new}' under {root}")

        for p in dirs + files:
            target = p.with_name(p.name.replace(old, new))
            kind = "dir " if p.is_dir() else "file"
            if target.exists():
                log(f"  SKIP {kind} {p.relative_to(root)}  ({target.name} already exists)")
                totals["skipped"] += 1
                continue
            log(f"  {kind} {p.relative_to(root)}  ->  {target.name}")
            if apply:
                p.rename(target)
            totals["directories" if kind == "dir " else "files"] += 1

        # Re-plan the contents: the paths just moved.
        if apply:
            _, _, contents = _plan(root, old, new)
        for p in contents:
            try:
                text = p.read_text(errors="ignore")
            except OSError:
                continue
            n = text.count(old)
            log(f"  text {p.relative_to(root)}  ({n} occurrence{'s' if n != 1 else ''})")
            if apply:
                p.write_text(text.replace(old, new))
            totals["contents"] += 1

    verb = "renamed" if apply else "would rename"
    log(f"{verb} {totals['directories']} directories and {totals['files']} files; "
        f"{'rewrote' if apply else 'would rewrite'} {totals['contents']} files"
        + (f"; {totals['skipped']} skipped" if totals["skipped"] else ""))
    if not apply and any(totals[k] for k in ("directories", "files", "contents")):
        log("this was a dry run - pass --apply to perform it")
    return totals
