"""Shared I/O and sequence helpers for the Morpheus pipeline.

Nothing in this module knows about genes or species; it is deliberately generic
so that the analysis modules stay readable.
"""
from __future__ import annotations

import csv
import gzip
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

# --------------------------------------------------------------------------
# logging
# --------------------------------------------------------------------------


def log(msg: str) -> None:
    print(f"[morpheus] {msg}", file=sys.stderr, flush=True)


def die(msg: str) -> "NoReturn":  # noqa: F821
    raise SystemExit(f"[morpheus] ERROR: {msg}")


# --------------------------------------------------------------------------
# files
# --------------------------------------------------------------------------


def opn(path, mode: str = "rt"):
    """Open a plain or gzipped file transparently."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_lines(path) -> List[str]:
    out = []
    with opn(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                out.append(line)
    return out


def read_tsv(path, required: Sequence[str] = ()) -> List[dict]:
    """Read a headed TSV into a list of dicts."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        rows = [dict(r) for r in rdr]
        missing = [c for c in required if rows and c not in rows[0]]
        if missing:
            die(f"{path} is missing column(s): {', '.join(missing)}")
        return rows


def write_tsv(path, rows: Iterable[dict], fields: Sequence[str]) -> Path:
    """Write rows as a headed TSV. Always writes the header, even if empty."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(fields), delimiter="\t",
                           extrasaction="ignore", restval="NA")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return path


# --------------------------------------------------------------------------
# FASTA
# --------------------------------------------------------------------------


def fasta_iter(path) -> Iterator[Tuple[str, str]]:
    """Yield (header_without_gt, sequence) pairs."""
    header, chunks = None, []
    with opn(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def read_fasta(path) -> Dict[str, str]:
    return {h.split()[0]: s for h, s in fasta_iter(path)}


def write_fasta(path, records: Iterable[Tuple[str, str]], wrap: int = 60) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            if wrap and wrap > 0:
                for i in range(0, len(seq), wrap):
                    fh.write(seq[i:i + wrap] + "\n")
            else:
                fh.write(seq + "\n")
    return path


# --------------------------------------------------------------------------
# sequence
# --------------------------------------------------------------------------

STOP_CODONS = {"TGA", "TAA", "TAG"}

_CODON_TABLE = {}
for _i, _b1 in enumerate("TCAG"):
    for _b2 in "TCAG":
        for _b3 in "TCAG":
            _CODON_TABLE[_b1 + _b2 + _b3] = (
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
            )[len(_CODON_TABLE)]

_COMPLEMENT = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv",
                            "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def clean_nt(seq: str) -> str:
    """Uppercase, strip whitespace, U->T, non-ACGTN -> N."""
    seq = re.sub(r"\s+", "", seq).upper().replace("U", "T")
    return re.sub(r"[^ACGTN]", "N", seq)


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def translate(seq: str, to_stop: bool = False) -> str:
    """Translate a CDS in frame 1. Incomplete/ambiguous codons become 'X'."""
    seq = clean_nt(seq)
    out = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        aa = _CODON_TABLE.get(seq[i:i + 3], "X")
        if to_stop and aa == "*":
            break
        out.append(aa)
    return "".join(out)


def codons(seq: str) -> List[str]:
    seq = clean_nt(seq)
    return [seq[i:i + 3] for i in range(0, len(seq), 3)]


def orf_report(seq: str) -> dict:
    """Describe the reading frame of a CDS.

    TOGA2 query CDS are not always a whole number of codons and frequently
    carry a few bases past the genuine stop, because the projection follows the
    query's own indels. Judging "complete" by whether the *final* codon is a
    stop therefore mislabels intact genes as truncated. Instead the ORF is taken
    to end at the first in-frame stop, and what matters is how much of the
    supplied sequence that ORF accounts for.
    """
    seq = clean_nt(seq)
    whole = seq[:len(seq) - len(seq) % 3]
    cs = [whole[i:i + 3] for i in range(0, len(whole), 3)]
    if not cs:
        return {"n_codons": 0, "has_start": False, "stop_codon_index": None,
                "orf_codons": 0, "orf_fraction": 0.0, "trailing_bp": len(seq),
                "internal_stops_before_orf_end": 0, "in_frame": len(seq) % 3 == 0}

    stop_index = next((i for i, c in enumerate(cs) if c in STOP_CODONS), None)
    orf_codons = stop_index if stop_index is not None else len(cs)
    return {
        "n_codons": len(cs),
        "has_start": cs[0] == "ATG",
        "stop_codon_index": stop_index,
        "orf_codons": orf_codons,
        # how much of the supplied CDS the ORF accounts for
        "orf_fraction": (stop_index + 1) / len(cs) if stop_index is not None else 1.0,
        "trailing_bp": len(seq) - (stop_index + 1) * 3 if stop_index is not None else len(seq) % 3,
        "in_frame": len(seq) % 3 == 0,
    }


# --------------------------------------------------------------------------
# identifiers
# --------------------------------------------------------------------------


def strip_version(tx_id: str) -> str:
    """ENST00000202917.10 -> ENST00000202917 ; NM_001320151.2 -> NM_001320151"""
    return re.sub(r"\.\d+$", "", str(tx_id).strip())


def sanitize(text: str, fallback: str = "NA") -> str:
    """Make a token safe for filenames and FASTA/Newick headers."""
    text = re.sub(r"\s+", "_", str(text).strip())
    text = re.sub(r"[^A-Za-z0-9_.-]", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or fallback


def species_from_dirname(dirname: str,
                         tree_tips: Optional[Iterable[str]] = None) -> str:
    """Bat1K/TOGA2 run directory name -> tree tip label.

    'Antrozous_pallidus__pallid_bat__HLantPal2__mAntPal2.1.pri'
        -> 'Antrozous_pallidus'

    Most names are binomial, but some tree tips are trinomial
    ('Rhinolophus_perniger_lanosus'). Taking the first two tokens would silently
    drop those species, so when the tree tips are known the longest tip that the
    directory name starts with wins; the binomial guess is only the fallback.
    """
    head = sanitize(str(dirname).rstrip("/").split("/")[-1].split("__")[0])
    if tree_tips:
        matches = [t for t in tree_tips
                   if head == t or head.startswith(t + "_")]
        if matches:
            return max(matches, key=len)
    parts = [p for p in head.split("_") if p]
    return "_".join(parts[:2]) if len(parts) >= 2 else head


def assembly_from_dirname(dirname: str) -> str:
    """Trailing assembly token of a Bat1K run directory, or 'NA'."""
    parts = str(dirname).rstrip("/").split("/")[-1].split("__")
    return sanitize(parts[-1]) if len(parts) >= 2 else "NA"


# --------------------------------------------------------------------------
# TOGA2 projection names
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class ProjectionName:
    """A TOGA2 projection label: <transcript>#<gene>#<chain>[#retro]."""
    raw: str
    transcript: str          # e.g. ENST00000202917.10
    transcript_base: str     # e.g. ENST00000202917
    gene: str                # gene symbol TOGA2 attached (may be wrong!)
    chain: str
    is_retro: bool

    @property
    def short(self) -> str:
        return f"{self.transcript}#{self.gene}#{self.chain}"


def parse_projection(name: str) -> ProjectionName:
    raw = str(name).split()[0]
    parts = raw.split("#")
    tx = parts[0] if parts else raw
    gene = parts[1] if len(parts) > 1 else "NA"
    chain = parts[2] if len(parts) > 2 else "NA"
    retro = any(p.lower() in {"retro", "pp"} for p in parts[3:])
    return ProjectionName(raw, tx, strip_version(tx), gene, chain, retro)


# --------------------------------------------------------------------------
# BED12
# --------------------------------------------------------------------------


@dataclass
class Bed12:
    chrom: str
    start: int          # 0-based, inclusive
    end: int            # 0-based, exclusive
    name: str
    score: str
    strand: str
    thick_start: int
    thick_end: int
    rgb: str
    block_count: int
    block_sizes: List[int]
    block_starts: List[int]

    @property
    def blocks(self) -> List[Tuple[int, int]]:
        return [(self.start + s, self.start + s + z)
                for s, z in zip(self.block_starts, self.block_sizes)]

    @property
    def blocks_transcript_order(self) -> List[Tuple[int, int]]:
        b = self.blocks
        return b if self.strand != "-" else list(reversed(b))

    @property
    def cds_bp(self) -> int:
        return sum(max(0, min(e, self.thick_end) - max(s, self.thick_start))
                   for s, e in self.blocks)

    @property
    def span(self) -> int:
        return self.end - self.start


def parse_bed12(path) -> List[Bed12]:
    out: List[Bed12] = []
    with opn(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            try:
                sizes = [int(x) for x in f[10].rstrip(",").split(",") if x != ""]
                starts = [int(x) for x in f[11].rstrip(",").split(",") if x != ""]
                out.append(Bed12(f[0], int(f[1]), int(f[2]), f[3], f[4], f[5],
                                 int(f[6]), int(f[7]), f[8], int(f[9]), sizes, starts))
            except ValueError as exc:
                die(f"malformed BED12 line {n} of {path}: {exc}")
    return out


def parse_bed6(path) -> List[dict]:
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with opn(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                rows.append({"chrom": f[0], "start": int(f[1]), "end": int(f[2]),
                             "name": f[3] if len(f) > 3 else "NA",
                             "strand": f[5] if len(f) > 5 else "."})
            except ValueError:
                continue
    rows.sort(key=lambda r: (r["chrom"], r["start"], r["end"]))
    return rows


# --------------------------------------------------------------------------
# intervals
# --------------------------------------------------------------------------


def overlap_bp(a: Tuple[int, int], b: Tuple[int, int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def reciprocal_overlap(a: Tuple[int, int], b: Tuple[int, int]) -> float:
    ov = overlap_bp(a, b)
    if ov <= 0:
        return 0.0
    return ov / max(1, min(a[1] - a[0], b[1] - b[0]))


def merge_intervals(intervals: Sequence[Tuple[int, int]], slop: int = 0
                    ) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    out = []
    for s, e in sorted(intervals):
        if out and s - slop <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return [tuple(x) for x in out]


# --------------------------------------------------------------------------
# misc
# --------------------------------------------------------------------------


def to_float(x, default: float = 0.0) -> float:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return default
    return default if v != v else v  # NaN guard


def to_int(x, default: int = 0) -> int:
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return default


def pct(numerator: float, denominator: float) -> float:
    return 100.0 * numerator / denominator if denominator else 0.0
