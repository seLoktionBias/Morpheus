"""Step 7 - HyPhy selection analyses that need no branch labelling.

  aBSREL  adaptive branch-site REL: which branches show episodic diversifying
          selection, with no a-priori foreground
  BUSTED  gene-wide test for episodic diversifying selection on any branch
  MEME    sites under episodic diversifying selection
  FEL     sites under pervasive selection, positive or negative

All four run on the unlabelled species tree, so nothing here depends on
choosing foreground lineages in advance.

The JSON HyPhy emits is not something anyone wants to read directly, so each
run is flattened into three tidy tables: one row per branch, one per gene, and
one per selected site.
"""
from __future__ import annotations

import json
import shutil
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .common import die, fasta_iter, log, read_tsv, to_float, to_int, write_tsv
from . import newick

METHODS = ("absrel", "busted", "meme", "fel")

GENE_FIELDS = ["gene", "human_transcript", "method", "n_sequences", "n_sites",
               "test_statistic", "p_value", "significant_at_0.05",
               "n_branches_tested", "n_branches_significant", "runtime_status"]

BRANCH_FIELDS = ["gene", "human_transcript", "species", "branch",
                 "is_terminal", "corrected_p_value", "uncorrected_p_value",
                 "lrt", "significant_at_0.05", "omega_classes",
                 "max_omega", "proportion_max_omega"]

SITE_FIELDS = ["gene", "human_transcript", "method", "codon_site", "alpha",
               "beta", "p_value", "significant_at_0.05", "selection_direction",
               "n_branches_under_selection"]


def _require_hyphy() -> str:
    exe = shutil.which("hyphy")
    if exe is None:
        die("hyphy not on PATH. Activate the analysis environment.")
    return exe


def _tree_for_alignment(alignment: Path, tree_path: Path, out_path: Path
                        ) -> Optional[Path]:
    """Prune the transcript tree to exactly the taxa present in the alignment."""
    names = [h.split()[0] for h, _ in fasta_iter(alignment)]
    tips = newick.prune_to_file(tree_path, names, out_path)
    if tips is None or len(tips) < 3:
        return None
    missing = set(names) - set(tips)
    if missing:
        log(f"  {alignment.name}: {len(missing)} alignment taxa absent from the "
            f"tree: {', '.join(sorted(missing))}")
        return None
    return out_path


def run_one(alignment, tree, outdir, method: str, threads: int = 2,
            timeout: Optional[int] = None) -> Optional[Path]:
    """Run one HyPhy method. Returns the JSON path, or None on failure."""
    hyphy = _require_hyphy()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_json = outdir / f"{method}.json"
    # Resume only from a JSON that actually parses. An interrupted run leaves a
    # non-empty but truncated file behind, and treating that as finished would
    # silently drop the analysis on every later attempt.
    if out_json.exists() and out_json.stat().st_size > 0:
        if _load(out_json) is not None:
            return out_json
        log(f"  {method}: discarding incomplete {out_json.name} and rerunning")
        out_json.unlink()

    cmd = [hyphy, "CPU=%d" % threads, method,
           "--alignment", str(alignment), "--tree", str(tree),
           "--output", str(out_json)]
    log_path = outdir / f"{method}.log"
    try:
        with open(log_path, "w") as lh:
            subprocess.run(cmd, check=True, stdout=lh, stderr=subprocess.STDOUT,
                           timeout=timeout)
    except subprocess.TimeoutExpired:
        log(f"  {method} timed out after {timeout}s for {Path(alignment).name}")
        return None
    except subprocess.CalledProcessError:
        tail = log_path.read_text()[-300:] if log_path.exists() else ""
        log(f"  {method} failed for {Path(alignment).name}: {tail.strip()[:200]}")
        return None
    return out_json if out_json.exists() else None


# --------------------------------------------------------------------------
# JSON -> tidy tables
# --------------------------------------------------------------------------


def _load(path) -> Optional[dict]:
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def _n_sites(doc: dict) -> int:
    info = doc.get("input", {})
    return to_int(info.get("number of sites", 0))


def _n_seqs(doc: dict) -> int:
    return to_int(doc.get("input", {}).get("number of sequences", 0))


def parse_absrel(doc: dict, gene: str, tx: str) -> Tuple[dict, List[dict]]:
    attrs = doc.get("branch attributes", {}).get("0", {})
    tested = doc.get("tested", {}).get("0", {})
    branches: List[dict] = []
    n_sig = 0
    for branch, vals in sorted(attrs.items()):
        if tested and tested.get(branch) != "test":
            continue
        corrected = vals.get("Corrected P-value")
        uncorrected = vals.get("Uncorrected P-value")
        rates = vals.get("Rate Distributions") or []
        omegas = [(to_float(r[0]), to_float(r[1])) for r in rates
                  if isinstance(r, (list, tuple)) and len(r) >= 2]
        max_omega, prop = max(omegas, key=lambda x: x[0]) if omegas else (0.0, 0.0)
        sig = int(corrected is not None and to_float(corrected, 1.0) < 0.05)
        n_sig += sig
        branches.append({
            "gene": gene, "human_transcript": tx,
            "species": branch if "Node" not in branch else "NA",
            "branch": branch,
            "is_terminal": int("Node" not in branch),
            "corrected_p_value": corrected if corrected is not None else "NA",
            "uncorrected_p_value": uncorrected if uncorrected is not None else "NA",
            "lrt": vals.get("LRT", "NA"),
            "significant_at_0.05": sig,
            "omega_classes": len(omegas),
            "max_omega": round(max_omega, 4),
            "proportion_max_omega": round(prop, 4),
        })
    gene_row = {
        "gene": gene, "human_transcript": tx, "method": "absrel",
        "n_sequences": _n_seqs(doc), "n_sites": _n_sites(doc),
        "test_statistic": "NA", "p_value": "NA",
        "significant_at_0.05": int(n_sig > 0),
        "n_branches_tested": len(branches), "n_branches_significant": n_sig,
        "runtime_status": "OK",
    }
    return gene_row, branches


def parse_busted(doc: dict, gene: str, tx: str) -> dict:
    tr = doc.get("test results", {})
    p = tr.get("p-value")
    return {
        "gene": gene, "human_transcript": tx, "method": "busted",
        "n_sequences": _n_seqs(doc), "n_sites": _n_sites(doc),
        "test_statistic": tr.get("LRT", "NA"),
        "p_value": p if p is not None else "NA",
        "significant_at_0.05": int(p is not None and to_float(p, 1.0) < 0.05),
        "n_branches_tested": len(doc.get("tested", {}).get("0", {})),
        "n_branches_significant": "NA", "runtime_status": "OK",
    }


def _site_table(doc: dict) -> Tuple[List[str], List[list]]:
    content = doc.get("MLE", {}).get("content", {}).get("0", [])
    headers = [h[0] for h in doc.get("MLE", {}).get("headers", [])]
    return headers, content


def parse_sites(doc: dict, gene: str, tx: str, method: str) -> Tuple[dict, List[dict]]:
    headers, content = _site_table(doc)
    idx = {h.lower(): i for i, h in enumerate(headers)}

    def col(*names) -> Optional[int]:
        for n in names:
            if n in idx:
                return idx[n]
        return None

    i_alpha = col("alpha", "&alpha;", "synonymous rate")
    i_beta = col("beta", "&beta;", "non-synonymous rate",
                 "&beta;<sup>+</sup>")
    i_p = col("p-value", "p_value", "p")
    i_nb = col("# branches under selection")

    rows: List[dict] = []
    n_sig = 0
    for site, vals in enumerate(content, 1):
        p = to_float(vals[i_p], 1.0) if i_p is not None and i_p < len(vals) else 1.0
        if p >= 0.05:
            continue
        alpha = to_float(vals[i_alpha]) if i_alpha is not None and i_alpha < len(vals) else 0.0
        beta = to_float(vals[i_beta]) if i_beta is not None and i_beta < len(vals) else 0.0
        n_sig += 1
        n_branches = (to_int(vals[i_nb]) if i_nb is not None and i_nb < len(vals)
                      else "NA")
        # MEME only tests for episodic positive selection, so its significant
        # sites are positive by construction; FEL is two-sided.
        direction = "positive" if method == "meme" else (
            "positive" if beta > alpha else "negative" if beta < alpha else "neutral")
        rows.append({
            "gene": gene, "human_transcript": tx, "method": method,
            "codon_site": site, "alpha": round(alpha, 4), "beta": round(beta, 4),
            "p_value": p, "significant_at_0.05": 1,
            "selection_direction": direction,
            "n_branches_under_selection": n_branches,
        })
    gene_row = {
        "gene": gene, "human_transcript": tx, "method": method,
        "n_sequences": _n_seqs(doc), "n_sites": _n_sites(doc),
        "test_statistic": "NA", "p_value": "NA",
        "significant_at_0.05": int(n_sig > 0),
        "n_branches_tested": "NA", "n_branches_significant": "NA",
        "runtime_status": f"{n_sig}_significant_sites",
    }
    return gene_row, rows


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------


def _select_alignments(manifest: List[dict], only_passing: bool,
                       longest_only: bool) -> List[dict]:
    """Which transcript directories to analyse."""
    rows = [r for r in manifest
            if not only_passing or r.get("passes_species_threshold") == "1"]
    if not longest_only:
        return rows
    # One transcript per gene: the one with the most species, then the longest
    # alignment, so the analysed model is the best-supported one.
    best: Dict[str, dict] = {}
    for r in rows:
        cur = best.get(r["gene"])
        key = (to_int(r.get("n_query_species")), to_int(r.get("n_sequences")))
        if cur is None or key > (to_int(cur.get("n_query_species")),
                                 to_int(cur.get("n_sequences"))):
            best[r["gene"]] = r
    return [best[g] for g in sorted(best)]


def run_all(manifest_tsv, outdir, methods: Sequence[str] = METHODS,
            threads: int = 2, timeout: Optional[int] = None,
            only_passing: bool = True, jobs: int = 1,
            longest_only: bool = False) -> Path:
    """Run the requested methods over every selected alignment.

    HyPhy scales poorly past a handful of cores, so the useful parallelism is
    across analyses rather than inside one. `jobs` analyses run at once, each
    given `threads` cores.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = read_tsv(manifest_tsv, required=["gene", "human_transcript", "fasta"])
    selected = _select_alignments(manifest, only_passing, longest_only)

    tasks = []
    for row in selected:
        gene, tx = row["gene"], row["human_transcript"]
        fasta = Path(row.get("selected_fasta") or row["fasta"])
        if str(fasta) == "NA":
            fasta = Path(row["fasta"])
        alignment = fasta.with_name(fasta.name.replace(".cds.fa", ".codon.aln.fa"))
        if not alignment.exists():
            log(f"  {gene} {tx}: no codon alignment - skipped")
            continue
        tree_src = row.get("selected_newick") or row.get("newick", "NA")
        if tree_src in ("NA", "", None):
            tree_src = row.get("newick", "NA")
        if tree_src == "NA":
            log(f"  {gene} {tx}: no pruned tree - skipped")
            continue
        work = alignment.parent / "hyphy"
        work.mkdir(parents=True, exist_ok=True)
        tree = _tree_for_alignment(alignment, Path(tree_src),
                                   work / "tree_for_hyphy.nwk")
        if tree is None:
            log(f"  {gene} {tx}: tree and alignment taxa do not match - skipped")
            continue
        for method in methods:
            tasks.append({"gene": gene, "tx": tx, "alignment": alignment,
                          "tree": tree, "work": work, "method": method})

    if not tasks:
        log("no alignments to analyse")
        return write_tsv(outdir / "hyphy_gene_level.tsv", [], GENE_FIELDS)

    log(f"running {len(tasks)} analyses over {len(selected)} alignment(s), "
        f"{jobs} at a time x {threads} core(s)")

    results: Dict[int, Optional[Path]] = {}
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        futures = {}
        for i, t in enumerate(tasks):
            futures[pool.submit(run_one, t["alignment"], t["tree"], t["work"],
                                t["method"], threads, timeout)] = i
        done = 0
        for fut in as_completed(futures):
            i = futures[fut]
            try:
                results[i] = fut.result()
            except Exception as exc:                       # noqa: BLE001
                log(f"  {tasks[i]['method']} raised for {tasks[i]['gene']} "
                    f"{tasks[i]['tx']}: {exc}")
                results[i] = None
            done += 1
            if done % 10 == 0 or done == len(tasks):
                log(f"  {done}/{len(tasks)} analyses finished")

    gene_rows: List[dict] = []
    branch_rows: List[dict] = []
    site_rows: List[dict] = []

    for i, t in enumerate(tasks):
        gene, tx, method = t["gene"], t["tx"], t["method"]
        json_path = results.get(i)
        if json_path is None:
            gene_rows.append({"gene": gene, "human_transcript": tx,
                              "method": method, "runtime_status": "FAILED",
                              "significant_at_0.05": "NA"})
            continue
        doc = _load(json_path)
        if doc is None:
            gene_rows.append({"gene": gene, "human_transcript": tx,
                              "method": method,
                              "runtime_status": "UNPARSEABLE_JSON",
                              "significant_at_0.05": "NA"})
            continue
        if method == "absrel":
            g, b = parse_absrel(doc, gene, tx)
            gene_rows.append(g)
            branch_rows.extend(b)
        elif method == "busted":
            gene_rows.append(parse_busted(doc, gene, tx))
        else:
            g, srows = parse_sites(doc, gene, tx, method)
            gene_rows.append(g)
            site_rows.extend(srows)

    gene_rows.sort(key=lambda r: (r["gene"], r["human_transcript"], r["method"]))
    write_tsv(outdir / "hyphy_gene_level.tsv", gene_rows, GENE_FIELDS)
    write_tsv(outdir / "hyphy_absrel_branches.tsv", branch_rows, BRANCH_FIELDS)
    write_tsv(outdir / "hyphy_selected_sites.tsv", site_rows, SITE_FIELDS)

    n_sig = sum(1 for r in gene_rows if r.get("significant_at_0.05") == 1)
    n_failed = sum(1 for r in gene_rows if r.get("runtime_status") == "FAILED")
    log(f"HyPhy: {len(gene_rows)} analyses, {n_sig} significant, {n_failed} failed")
    return outdir / "hyphy_gene_level.tsv"
