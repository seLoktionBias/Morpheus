"""Step 6 - codon-aware alignment of each transcript multifasta.

CDS sequences are translated, the proteins are aligned with MAFFT (L-INS-i),
and the nucleotides are threaded back onto the protein alignment. Aligning
amino acids and back-translating is frame-safe: it cannot open a gap that is
not a multiple of three, which is what HyPhy's codon models require.

TOGA2 models are often frameshifted or carry premature stops. Rather than
discarding them, each sequence is trimmed to a whole number of codons and its
internal stop codons are masked to NNN before translation, so a single disrupted
codon does not truncate the whole protein and wreck the alignment. The masking
is recorded per sequence so nothing is lost silently.
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .common import (STOP_CODONS, clean_nt, codons, die, fasta_iter, log,
                     read_tsv, translate, write_fasta, write_tsv)

PREP_FIELDS = ["gene", "human_transcript", "species", "raw_length",
               "used_length", "trimmed_bp", "n_internal_stops_masked",
               "n_ambiguous_codons", "terminal_stop_removed"]


def _require_mafft() -> str:
    exe = shutil.which("mafft")
    if exe is None:
        die("mafft not on PATH. Activate the analysis environment.")
    return exe


def prepare_cds(seq: str) -> Tuple[str, dict]:
    """Trim to whole codons, drop a terminal stop, mask internal stops."""
    seq = clean_nt(seq)
    raw_len = len(seq)
    trimmed = raw_len - raw_len % 3
    seq = seq[:trimmed]
    cs = codons(seq)

    terminal_stop = bool(cs) and cs[-1] in STOP_CODONS
    if terminal_stop:
        cs = cs[:-1]

    masked = 0
    ambiguous = 0
    for i, c in enumerate(cs):
        if c in STOP_CODONS:
            cs[i] = "NNN"
            masked += 1
        elif "N" in c:
            ambiguous += 1

    out = "".join(cs)
    return out, {
        "raw_length": raw_len,
        "used_length": len(out),
        "trimmed_bp": raw_len - trimmed,
        "n_internal_stops_masked": masked,
        "n_ambiguous_codons": ambiguous,
        "terminal_stop_removed": int(terminal_stop),
    }


def backtranslate(protein_alignment: Dict[str, str],
                  nucleotides: Dict[str, str]) -> List[Tuple[str, str]]:
    """Thread codons back onto an aligned protein, gap for gap."""
    out = []
    for name, aligned_aa in protein_alignment.items():
        nt = nucleotides.get(name)
        if nt is None:
            log(f"  back-translation: aligner renamed '{name}' - sequence dropped")
            continue
        pieces, pos = [], 0
        for aa in aligned_aa:
            if aa == "-":
                pieces.append("---")
            else:
                pieces.append(nt[pos:pos + 3].ljust(3, "N"))
                pos += 3
        out.append((name, "".join(pieces)))
    return out


def align_file(fasta_path, out_path, threads: int = 2,
               mafft_args: Optional[Sequence[str]] = None) -> Optional[dict]:
    """Codon-align one transcript multifasta. Returns per-sequence prep stats."""
    mafft = _require_mafft()
    fasta_path, out_path = Path(fasta_path), Path(out_path)

    records = list(fasta_iter(fasta_path))
    if len(records) < 3:
        log(f"  {fasta_path.name}: only {len(records)} sequence(s) - not aligned")
        return None

    prepared: Dict[str, str] = {}
    stats: Dict[str, dict] = {}
    for header, seq in records:
        name = header.split()[0]
        nt, st = prepare_cds(seq)
        if len(nt) < 30:
            log(f"  {fasta_path.name}: {name} shorter than 10 codons - dropped")
            continue
        prepared[name] = nt
        stats[name] = st

    if len(prepared) < 3:
        log(f"  {fasta_path.name}: fewer than 3 usable sequences - not aligned")
        return None

    with tempfile.TemporaryDirectory() as tmp:
        prot_in = Path(tmp) / "prot.fa"
        write_fasta(prot_in, [(n, translate(s)) for n, s in prepared.items()])
        cmd = [mafft, "--anysymbol", "--quiet", "--thread", str(threads)]
        cmd += list(mafft_args) if mafft_args else ["--maxiterate", "1000", "--localpair"]
        cmd.append(str(prot_in))
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            log(f"  {fasta_path.name}: mafft failed - {result.stderr.strip()[:200]}")
            return None
        aln_path = Path(tmp) / "prot.aln"
        aln_path.write_text(result.stdout)
        protein_alignment = {h.split()[0]: s.upper() for h, s in fasta_iter(aln_path)}

    codon_alignment = backtranslate(protein_alignment, prepared)
    # keep input order stable and deterministic
    order = {n: i for i, n in enumerate(prepared)}
    codon_alignment.sort(key=lambda kv: order.get(kv[0], 1 << 30))
    write_fasta(out_path, codon_alignment)
    return stats


def align_all(manifest_tsv, outdir=None, threads: int = 2,
              only_passing: bool = True) -> Path:
    """Align every transcript directory listed in the deliverables manifest."""
    manifest = read_tsv(manifest_tsv, required=["gene", "human_transcript", "fasta"])
    prep_rows: List[dict] = []
    n_ok = n_skip = 0

    for row in manifest:
        if only_passing and row.get("passes_species_threshold") != "1":
            n_skip += 1
            continue
        # align inside the selected set, so its alignments travel with it
        fasta = Path(row.get("selected_fasta") or row["fasta"])
        if str(fasta) == "NA":
            fasta = Path(row["fasta"])
        if not fasta.exists():
            n_skip += 1
            continue
        out_path = fasta.with_name(fasta.name.replace(".cds.fa", ".codon.aln.fa"))
        stats = align_file(fasta, out_path, threads=threads)
        if stats is None:
            n_skip += 1
            continue
        n_ok += 1
        for species, st in stats.items():
            prep_rows.append({"gene": row["gene"],
                              "human_transcript": row["human_transcript"],
                              "species": species, **st})
        log(f"  aligned {row['gene']} {row['human_transcript']} "
            f"({len(stats)} sequences)")

    out = Path(outdir) if outdir else Path(manifest_tsv).parent
    write_tsv(out / "alignment_preparation.tsv", prep_rows, PREP_FIELDS)
    log(f"codon alignments: {n_ok} written, {n_skip} skipped")
    return out / "alignment_preparation.tsv"
