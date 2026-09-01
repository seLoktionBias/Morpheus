"""Step 3 - pairwise BLAST and one-to-one human/bat transcript assignment.

The assignment is run twice, because two different questions need two different
search scopes.

REGION_RESTRICTED
    Only projections inside the locus of the projected longest human isoform.
    This answers "what does this gene's own syntenic locus produce?", which is
    what transcript status must report: a paralog elsewhere in the genome cannot
    stand in for a transcript the locus itself has lost.

UNRESTRICTED
    Projections of the gene anywhere in the genome, but still delimited to the
    gene itself - a projection that belongs to a sibling paralog is excluded by
    the whole-proteome screen wherever it sits. This answers "does the animal
    make this protein at all?", which is the right scope for alignment and
    selection analysis.

The two genuinely disagree, and the disagreement is the biology. OAS1 in
Phyllostomus discolor is the motivating case: the transcript carrying the
C-terminal CaaX domain is degraded at the syntenic locus, but a copy on another
scaffold carries it at ~98% identity to human. Region-restricted, the animal
looks like it has lost the domain; unrestricted, it plainly has not. Conversely,
dropping the region restriction entirely would let the IFITM tandem array mix
paralogs together, which is why the screen still delimits by gene.

`compare.py` reports, per transcript, whether the two scopes chose the same
model.

Assignment happens in two passes, in this order.

  1. PROJECTION IDENTITY. TOGA2 names every projection
     `<human_transcript>#<gene>#<chain>`, so a projection already records which
     human transcript it was built from. That is not an estimate of
     correspondence, it is the construction, and it is the best evidence
     available. Each human transcript takes the projection built from it.

  2. SIMILARITY MATCHING for the remainder. Many human transcripts are never
     projected at all, and those are genuinely ambiguous. For them the full
     pairwise BLAST matrix is solved as a maximum-weight bipartite matching, so
     one query transcript serves at most one human transcript and vice versa.

Order matters here, but not iteration order: the Hungarian solution maximises
the total over all pairings at once and is provably independent of the order
rows and columns are presented in. What pass 1 does is stop similarity from
overriding a fact TOGA2 already established.

Without pass 1 the similarity matching recovered projection identity about three
times in four and overrode it in the rest. Without the one-to-one constraint in
pass 2, roughly twice as many assignments are made, because the same handful of
query models get copied into every isoform's file - the median group has six
human transcripts and two distinct query models to share between them.

Every row records which pass produced it in `assignment_basis`.

The whole-proteome screen still annotates every candidate, but it does not
filter here: inside the locus the region is the evidence, and the screen verdict
is reported as a column for inspection.
"""
from __future__ import annotations

import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .common import (die, fasta_iter, log, read_tsv, strip_version, to_float,
                     to_int, translate, write_fasta, write_tsv)
from .matching import max_weight_matching

BLAST_FIELDS = ("qseqid sseqid pident positive length qlen slen "
                "evalue bitscore qstart qend sstart send")

PAIR_FIELDS = ["species", "gene", "human_transcript", "candidate_id", "pool",
               "is_retro", "identity_global", "pident_local", "ppos",
               "alignment_coverage", "bitscore", "evalue",
               "structure_lcs_coverage", "structure_term_applied",
               "exon_label_jaccard", "length_similarity", "similarity_score"]

ASSIGN_FIELDS = ["species", "gene", "human_transcript", "candidate_id",
                 "projection", "toga_transcript", "toga_gene_label",
                 "toga_gene_label_matches_target", "toga_status",
                 "orthology_class", "pool", "chrom", "start", "end", "strand",
                 "upstream_gene", "downstream_gene",
                 "distance_to_home_locus_bp", "is_retro", "cds_bp",
                 "seq_length", "screen_verdict", "best_human_gene",
                 "pident", "identity_global", "ppos", "alignment_coverage",
                 "bitscore", "structure_lcs_coverage", "length_similarity",
                 "similarity_score", "assignment_basis",
                 "sequence_best_candidate", "sequence_best_identity",
                 "assigned_identity", "rank_for_human_transcript",
                 "best_alternative_candidate", "best_alternative_score",
                 "score_margin", "assignment_status", "note"]

# Search scopes. FAMILY-pool projections are never transcript candidates in
# either scope; they exist only for copy counting.
SCOPES = {
    "region_restricted": ("IN_REGION",),
    "unrestricted": ("IN_REGION", "OFF_REGION"),
}

# Screen verdicts that keep a candidate attached to the target gene. Applied in
# the unrestricted scope, where region no longer does the delimiting.
KEEP_VERDICTS = {"CONSISTENT_WITH_TARGET_GENE",
                 "CONSISTENT_BY_SYNTENY_WITHIN_FAMILY"}

# Weights: sequence similarity dominates, because the task is to pick the most
# similar CDS. Structure and length break ties between near-equal hits.
W_PIDENT, W_PPOS, W_COV, W_STRUCT, W_LEN = 0.45, 0.10, 0.25, 0.12, 0.08


# Two defensible ways to rank a candidate, kept side by side because they
# answer different questions and genuinely disagree.
#
# sequence_similarity
#     Which query sequence most resembles this human transcript? A retrocopy is
#     a reverse-transcribed mRNA with no introns, so its exon-structure term is
#     necessarily zero; scoring it there would penalise a processed copy for
#     being processed, which is circular. The term is dropped for retrocopies
#     and the remaining weights renormalised. This is the view that recovers the
#     Phyllostomus discolor OAS1 retrogene - the only model retaining the
#     C-terminal CaaX motif.
#
# structure_aware
#     Which query model has the *architecture* of the real gene? A processed
#     pseudogene can be the closest sequence and still never be transcribed,
#     having no promoter, enhancer or TSS. Here the structure term applies to
#     every candidate, so an intronless copy pays the full 0.12 for having no
#     exon structure to match, and the multi-exon model wins.
#
#     Named for what it does. It weighs *exon* structure and nothing positional:
#     no scoring term reads a coordinate, a flanking gene or a distance from the
#     home locus. The positional constraint on orthology is the other axis --
#     the region_restricted scope -- and keeping the two separate is the point.
#     (This policy was called synteny_aware before v2.1.0, which promised
#     positional reasoning the ranking never did.)
POLICIES = ("sequence_similarity", "structure_aware")


def score_components(identity_global: float, ppos: float, coverage: float,
                     structure: float, length_similarity: float,
                     is_retro: bool, policy: str = "sequence_similarity"
                     ) -> Tuple[float, dict]:
    """Combine the similarity terms into one score under the given policy."""
    terms = {"identity": (W_PIDENT, identity_global),
             "positives": (W_PPOS, ppos / 100.0),
             "coverage": (W_COV, coverage),
             "length": (W_LEN, length_similarity)}
    # Under structure_aware the structure term always applies - that is the whole
    # point of the policy. Under sequence_similarity it is dropped for the one
    # class of candidate that cannot possibly satisfy it.
    apply_structure = (policy == "structure_aware") or not is_retro
    if apply_structure:
        terms["structure"] = (W_STRUCT, structure)
    total_w = sum(w for w, _ in terms.values())
    score = sum(w * v for w, v in terms.values()) / total_w
    return score, {k: round(w * v / total_w, 4) for k, (w, v) in terms.items()}


# --------------------------------------------------------------------------
# BLAST
# --------------------------------------------------------------------------


def _require(*exes: str) -> None:
    missing = [e for e in exes if shutil.which(e) is None]
    if missing:
        die(f"executable(s) not on PATH: {', '.join(missing)}. "
            f"Activate the analysis environment first.")


def run_pairwise_blast(candidates_tsv, candidate_fasta, reference_dir, outdir,
                       threads: int = 4, evalue: str = "1e-3") -> Path:
    """blastx every bat candidate against the human proteins of its own gene."""
    _require("makeblastdb", "blastx")
    outdir = Path(outdir)
    work = outdir / "blast_work"
    work.mkdir(parents=True, exist_ok=True)

    cands = read_tsv(candidates_tsv, required=["candidate_id", "gene"])
    # Only the locus pool is a transcript candidate. Off-region and paralog-family
    # projections exist for copy counting and are not BLASTed here.
    gene_of_candidate = {r["candidate_id"]: r["gene"] for r in cands
                         if r.get("pool") in ("IN_REGION", "OFF_REGION")}
    seqs = {h.split()[0]: s for h, s in fasta_iter(candidate_fasta)}

    human_prot = {h.split()[0]: s for h, s in
                  fasta_iter(Path(reference_dir) / "human_transcript_proteins.fa")}
    ref_tx = read_tsv(Path(reference_dir) / "transcripts.tsv")
    human_by_gene: Dict[str, List[str]] = defaultdict(list)
    for r in ref_tx:
        human_by_gene[r["gene"]].append(r["transcript_id"])

    out_path = outdir / "pairwise_blastx.tsv"
    with open(out_path, "w") as combined:
        combined.write("\t".join(BLAST_FIELDS.split()) + "\n")
        for gene, human_ids in sorted(human_by_gene.items()):
            q_ids = [cid for cid, g in gene_of_candidate.items()
                     if g == gene and cid in seqs]
            if not q_ids or not human_ids:
                log(f"  {gene}: nothing to BLAST "
                    f"({len(q_ids)} candidates, {len(human_ids)} human transcripts)")
                continue

            db_fa = work / f"{gene}.human_proteins.fa"
            write_fasta(db_fa, [(t, human_prot[t]) for t in human_ids if t in human_prot])
            q_fa = work / f"{gene}.candidates.fa"
            write_fasta(q_fa, [(c, seqs[c]) for c in sorted(q_ids)])

            db = str(work / f"{gene}.db")
            subprocess.run(["makeblastdb", "-in", str(db_fa), "-dbtype", "prot",
                            "-out", db], check=True, stdout=subprocess.DEVNULL)
            hits = work / f"{gene}.blastx.tsv"
            subprocess.run(
                ["blastx", "-query", str(q_fa), "-db", db, "-out", str(hits),
                 "-evalue", evalue, "-seg", "no",
                 "-max_target_seqs", str(max(10, len(human_ids) * 2)),
                 "-num_threads", str(threads),
                 "-outfmt", "6 " + BLAST_FIELDS],
                check=True)
            n = 0
            with open(hits) as fh:
                for line in fh:
                    combined.write(line)
                    n += 1
            log(f"  {gene}: {len(q_ids)} candidates x {len(human_ids)} human "
                f"transcripts -> {n} HSP rows")
    log(f"pairwise BLAST written to {out_path}")
    return out_path


def _merge_hsps(hsps: List[dict], qlen_aa: float, slen: float) -> dict:
    """Combine every non-overlapping HSP of one pair into one measure.

    The single best HSP is a *local* view. Two models can have near-identical
    best HSPs while one of them carries far more indels and loses a terminal
    domain altogether; a local figure cannot see that, and ranks them the wrong
    way round. So HSPs are accepted greedily by bitscore, skipping any that
    overlap a query or subject range already taken, and the identity is then
    expressed over the *longer* of the two sequences - the analogue of counting
    gaps in a global alignment.
    """
    taken_q: List[Tuple[float, float]] = []
    taken_s: List[Tuple[float, float]] = []

    def overlaps(rng, taken):
        return any(min(rng[1], t[1]) - max(rng[0], t[0]) > 0 for t in taken)

    identities = aligned = positives = 0.0
    best_local = best_bits = 0.0
    evalue = "NA"
    for h in sorted(hsps, key=lambda x: -x["bitscore"]):
        if overlaps(h["qrange"], taken_q) or overlaps(h["srange"], taken_s):
            continue
        taken_q.append(h["qrange"]); taken_s.append(h["srange"])
        identities += h["pident"] / 100.0 * h["length"]
        positives += h["positive"]
        aligned += h["length"]
        if h["bitscore"] > best_bits:
            best_bits, best_local, evalue = h["bitscore"], h["pident"], h["evalue"]

    longer = max(qlen_aa, slen, 1.0)
    return {
        # identity over the longer sequence: partial or indel-riddled matches
        # are penalised the way a global alignment would penalise them
        "identity_global": identities / longer,
        "pident_local": best_local,
        "ppos": 100.0 * positives / aligned if aligned else 0.0,
        "coverage": min(aligned / max(qlen_aa, 1.0), aligned / max(slen, 1.0)),
        "aligned_aa": aligned,
        "evalue": evalue,
        "bitscore": best_bits,
    }


def load_best_hsps(blast_tsv) -> Dict[Tuple[str, str], dict]:
    """Aggregated alignment statistics per (candidate, human transcript)."""
    raw: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    lens: Dict[Tuple[str, str], Tuple[float, float]] = {}
    with open(blast_tsv) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 13:
                continue
            q, sub = f[0], f[1]
            aln = to_float(f[4])
            if aln <= 0:
                continue
            key = (q, sub)
            lens[key] = (max(1.0, to_float(f[5]) / 3.0), max(1.0, to_float(f[6])))
            qs, qe = sorted((to_float(f[9]), to_float(f[10])))
            ss, se = sorted((to_float(f[11]), to_float(f[12])))
            raw[key].append({
                "pident": to_float(f[2]), "positive": to_float(f[3]),
                "length": aln, "evalue": f[7], "bitscore": to_float(f[8]),
                "qrange": (qs, qe), "srange": (ss, se),
            })
    return {k: _merge_hsps(v, *lens[k]) for k, v in raw.items()}


# --------------------------------------------------------------------------
# structure similarity
# --------------------------------------------------------------------------


def _lcs_length(a: Sequence[str], b: Sequence[str]) -> int:
    if not a or not b:
        return 0
    prev = [0] * (len(b) + 1)
    for x in a:
        cur = [0]
        for j, y in enumerate(b, 1):
            cur.append(prev[j - 1] + 1 if x == y else max(prev[j], cur[j - 1]))
        prev = cur
    return prev[-1]


def structure_similarity(bat_labels: Sequence[str], human_labels: Sequence[str]) -> dict:
    """Ordered and set-wise agreement between two exon-label chains."""
    lcs = _lcs_length(bat_labels, human_labels)
    bat_set, human_set = set(bat_labels), set(human_labels)
    union = bat_set | human_set
    return {
        "structure_lcs_coverage": lcs / max(1, len(human_labels)),
        "exon_label_jaccard": len(bat_set & human_set) / max(1, len(union)),
    }


def _length_similarity(a: int, b: int) -> float:
    return 1.0 - abs(a - b) / max(1, a, b)


# --------------------------------------------------------------------------
# assignment
# --------------------------------------------------------------------------


def build_identity_map(ref_tx: List[dict]) -> Dict[str, str]:
    """Any human transcript id -> the representative we actually analyse.

    Transcripts with an identical CDS were collapsed to one representative, so a
    projection built from a collapsed id still belongs to its representative.
    """
    rep_of: Dict[str, str] = {}
    for r in ref_tx:
        rep = strip_version(r["transcript_id"])
        rep_of[rep] = rep
        for field in ("identical_model_transcripts",
                      "identical_cds_other_model_transcripts"):
            for other in (r.get(field) or "NONE").split(","):
                other = strip_version(other)
                if other and other != "NONE":
                    rep_of[other] = rep
    return rep_of


def assign(candidates_tsv, blast_tsv, reference_dir, outdir,
           screen_tsv=None, min_score: float = 0.30,
           min_pident: float = 40.0, scope: str = "region_restricted",
           identity_first: bool = True,
           policy: str = "sequence_similarity") -> Path:
    if policy not in POLICIES:
        die(f"unknown policy '{policy}'; expected one of {', '.join(POLICIES)}")
    """Assign query transcripts to human transcripts within one search scope."""
    if scope not in SCOPES:
        die(f"unknown scope '{scope}'; expected one of {', '.join(SCOPES)}")
    pools = SCOPES[scope]
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cands = read_tsv(candidates_tsv, required=["candidate_id", "species", "gene"])
    # Annotation only - the region restriction is what keeps paralogs out here.
    screen = {r["candidate_id"]: r for r in read_tsv(screen_tsv)} if screen_tsv else {}
    ref_tx = read_tsv(Path(reference_dir) / "transcripts.tsv")
    human_by_gene: Dict[str, List[dict]] = defaultdict(list)
    for r in ref_tx:
        human_by_gene[r["gene"]].append(r)
    rep_of = build_identity_map(ref_tx)

    hsps = load_best_hsps(blast_tsv)

    groups: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    excluded: List[dict] = []
    for r in cands:
        if to_int(r.get("seq_length")) <= 0:
            continue
        if r.get("pool") not in pools:
            excluded.append({**r, "reason": f"pool_not_searched_in_{scope}"})
            continue
        # Inside the locus the region delimits the gene. Outside it, only the
        # screen can, so a projection that belongs to a sibling is dropped.
        if r.get("pool") != "IN_REGION":
            verdict = screen.get(r["candidate_id"], {}).get("screen_verdict", "NA")
            if verdict not in KEEP_VERDICTS:
                excluded.append({**r, "screen_verdict": verdict,
                                 "reason": "off_locus_and_not_delimited_to_this_gene"})
                continue
        groups[(r["species"], r["gene"])].append(r)

    write_tsv(outdir / f"candidates_excluded_{scope}__{policy}.tsv", excluded,
              ["species", "gene", "candidate_id", "projection", "toga_gene_label",
               "pool", "chrom", "start", "end", "distance_to_home_locus_bp",
               "upstream_gene", "downstream_gene", "screen_verdict", "reason"])
    log(f"[{scope}__{policy}] {sum(len(v) for v in groups.values())} candidates "
        f"searched; {len(excluded)} excluded")

    pair_rows: List[dict] = []
    assign_rows: List[dict] = []

    for (species, gene), members in sorted(groups.items()):
        humans = sorted(human_by_gene.get(gene, []), key=lambda r: r["transcript_id"])
        if not humans:
            continue
        members = sorted(members, key=lambda r: r["candidate_id"])

        # ---- full pairwise similarity matrix ---------------------------
        matrix: List[List[float]] = []
        detail: List[List[dict]] = []
        # identity alone, ignoring structure and length, so a disagreement
        # between "best sequence match" and "what the scorer chose" is visible
        seq_only: Dict[int, List[float]] = defaultdict(list)
        for i, h in enumerate(humans):
            h_labels = [x for x in h["exon_structure"].split(",") if x]
            h_len = to_int(h["cds_sequence_length"])
            row_scores, row_detail = [], []
            for c in members:
                hsp = hsps.get((c["candidate_id"], h["transcript_id"]))
                ident = hsp["identity_global"] if hsp else 0.0
                pid_local = hsp["pident_local"] if hsp else 0.0
                ppos = hsp["ppos"] if hsp else 0.0
                cov = hsp["coverage"] if hsp else 0.0
                bits = hsp["bitscore"] if hsp else 0.0
                evalue = hsp["evalue"] if hsp else "NA"

                b_labels = [x for x in c["exon_structure"].split(",") if x]
                sm = structure_similarity(b_labels, h_labels)
                len_sim = _length_similarity(to_int(c["seq_length"]), h_len)
                is_retro = str(c.get("is_retro", "0")) == "1"

                score, _parts = score_components(
                    ident, ppos, cov, sm["structure_lcs_coverage"], len_sim,
                    is_retro, policy)
                # A pair with no alignment support, or an implausible protein
                # identity, is never allowed to win a slot.
                if bits <= 0 or pid_local < min_pident:
                    score = 0.0

                d = {"species": species, "gene": gene,
                     "human_transcript": h["transcript_id"],
                     "candidate_id": c["candidate_id"], "pool": c["pool"],
                     "is_retro": int(is_retro),
                     "identity_global": round(ident, 4),
                     "pident_local": round(pid_local, 3), "ppos": round(ppos, 3),
                     "alignment_coverage": round(cov, 4),
                     "bitscore": round(bits, 1), "evalue": evalue,
                     "structure_lcs_coverage": round(sm["structure_lcs_coverage"], 4),
                     "structure_term_applied":
                         int(policy == "structure_aware" or not is_retro),
                     "exon_label_jaccard": round(sm["exon_label_jaccard"], 4),
                     "length_similarity": round(len_sim, 4),
                     "similarity_score": round(score, 4)}
                row_scores.append(score)
                row_detail.append(d)
                # sequence-only view, for the disagreement flag below
                seq_only[i].append(ident if bits > 0 else 0.0)
                if bits > 0:
                    pair_rows.append(d)
            matrix.append(row_scores)
            detail.append(row_detail)

        # ---- pass 1: projection identity --------------------------------
        # A candidate maps to exactly one human transcript - the one it was
        # projected from - so identity is one-to-one by construction. Where two
        # chains project the same transcript, the better-scoring one wins and
        # the other returns to the pool for pass 2.
        identity_of: Dict[int, int] = {}          # candidate index -> human index
        human_index = {h["transcript_id"]: i for i, h in enumerate(humans)}
        claims: Dict[int, List[int]] = defaultdict(list)
        if identity_first:
            for j, c in enumerate(members):
                target = rep_of.get(strip_version(c.get("toga_transcript", "")))
                i = human_index.get(target) if target else None
                if i is not None:
                    claims[i].append(j)

        identity_pairs: Dict[int, int] = {}
        contested: Dict[int, List[int]] = {}
        for i, js in claims.items():
            # No similarity floor here: a degraded projection is still that
            # transcript's projection, and the status table needs to say so.
            js_sorted = sorted(js, key=lambda j: -matrix[i][j])
            identity_pairs[i] = js_sorted[0]
            identity_of[js_sorted[0]] = i
            if len(js_sorted) > 1:
                contested[i] = js_sorted[1:]

        # ---- pass 2: similarity matching for whatever is left ------------
        free_humans = [i for i in range(len(humans)) if i not in identity_pairs]
        taken = set(identity_pairs.values())
        free_cands = [j for j in range(len(members)) if j not in taken]

        chosen = dict(identity_pairs)
        if free_humans and free_cands:
            sub = [[matrix[i][j] for j in free_cands] for i in free_humans]
            for a, b in max_weight_matching(sub, min_profit=min_score):
                chosen[free_humans[a]] = free_cands[b]

        for i, h in enumerate(humans):
            ranked = sorted(range(len(members)), key=lambda j: -matrix[i][j])
            top = ranked[0] if ranked else None
            j = chosen.get(i)

            # The best candidate this human transcript did *not* get: the
            # second choice when it was assigned, the blocked first choice when
            # it was not. Same meaning either way.
            alt = next((k for k in ranked if k != j), None)
            alt_id = members[alt]["candidate_id"] if alt is not None else "NA"
            alt_score = matrix[i][alt] if alt is not None else 0.0

            if j is None:
                best_score = matrix[i][top] if top is not None else 0.0
                assign_rows.append({
                    "species": species, "gene": gene,
                    "human_transcript": h["transcript_id"],
                    "candidate_id": "NA", "projection": "NA",
                    "assignment_status": "UNASSIGNED",
                    "assignment_basis": "none",
                    "sequence_best_candidate": "NA",
                    "sequence_best_identity": "NA",
                    "assigned_identity": "NA",
                    "similarity_score": round(best_score, 4),
                    "best_alternative_candidate": alt_id,
                    "best_alternative_score": round(alt_score, 4),
                    "score_margin": "NA",
                    "note": ("no_projection_of_this_transcript_and_"
                             "no_candidate_above_min_score"
                             if best_score <= min_score
                             else "no_projection_of_this_transcript_and_best_"
                                  "candidate_taken_by_another"),
                })
                continue

            c = members[j]
            d = detail[i][j]
            by_identity = i in identity_pairs

            # Which candidate has the highest raw sequence identity, regardless
            # of structure, length or which slot it ended up in?
            seq_rank = sorted(range(len(members)), key=lambda k: -seq_only[i][k])
            seq_best = seq_rank[0] if seq_rank else None
            seq_best_id = members[seq_best]["candidate_id"] if seq_best is not None else "NA"
            seq_best_val = seq_only[i][seq_best] if seq_best is not None else 0.0

            # Notes accumulate: a match can be both a reassignment and
            # off-region, and losing either fact would be misleading.
            notes = []
            if by_identity:
                if i in contested:
                    notes.append(f"{len(contested[i]) + 1}_chains_project_this_"
                                 f"transcript_best_scoring_kept")
                if top is not None and top != j:
                    # Recorded, never acted on: TOGA2 built this projection from
                    # this transcript, and a higher BLAST score against a
                    # different model does not undo that.
                    notes.append("similarity_would_have_chosen_a_different_model")
            elif top is not None and top != j:
                notes.append("reassigned_from_greedy_best_to_satisfy_one_to_one_constraint")
            if str(c.get("toga_gene_label", "")).upper() != gene.upper():
                notes.append("toga2_filed_this_projection_under_a_different_gene")
            # Surface it whenever the assigned model is not the one with the
            # highest sequence identity - the reader can then judge for
            # themselves rather than discovering it by hand-aligning.
            if seq_best is not None and seq_best != j and seq_best_val > 0:
                margin = seq_best_val - d["identity_global"]
                if margin > 0.01:
                    notes.append(f"a_different_model_has_higher_sequence_identity"
                                 f"_by_{margin:.3f}")
            if not notes:
                notes.append("optimal_one_to_one_match")
            note = ";".join(notes)

            assign_rows.append({
                **{k: c.get(k, "NA") for k in
                   ("projection", "toga_transcript", "toga_gene_label",
                    "toga_status", "orthology_class", "pool", "chrom", "start",
                    "end", "strand", "upstream_gene", "downstream_gene",
                    "distance_to_home_locus_bp", "is_retro",
                    "cds_bp", "seq_length")},
                "species": species, "gene": gene,
                "human_transcript": h["transcript_id"],
                "candidate_id": c["candidate_id"],
                "screen_verdict": screen.get(c["candidate_id"], {}).get("screen_verdict", "NA"),
                "best_human_gene": screen.get(c["candidate_id"], {}).get("best_human_gene", "NA"),
                "toga_gene_label_matches_target":
                    int(str(c.get("toga_gene_label", "")).upper() == gene.upper()),
                "pident": d["pident_local"],
                "identity_global": d["identity_global"], "ppos": d["ppos"],
                "alignment_coverage": d["alignment_coverage"],
                "bitscore": d["bitscore"],
                "structure_lcs_coverage": d["structure_lcs_coverage"],
                "length_similarity": d["length_similarity"],
                "similarity_score": d["similarity_score"],
                "assignment_basis": ("projection_identity" if by_identity
                                     else "similarity_matching"),
                "sequence_best_candidate": seq_best_id,
                "sequence_best_identity": round(seq_best_val, 4),
                "assigned_identity": d["identity_global"],
                "rank_for_human_transcript": ranked.index(j) + 1,
                "best_alternative_candidate": alt_id,
                "best_alternative_score": round(alt_score, 4),
                "score_margin": round(matrix[i][j] - alt_score, 4),
                "assignment_status": "ASSIGNED",
                "note": note,
            })

    tag = f"{scope}__{policy}"
    write_tsv(outdir / f"pairwise_similarity_{tag}.tsv", pair_rows, PAIR_FIELDS)
    out_path = outdir / f"transcript_assignments_{tag}.tsv"
    write_tsv(out_path, assign_rows, ASSIGN_FIELDS)

    n_assigned = sum(1 for r in assign_rows if r["assignment_status"] == "ASSIGNED")
    n_ident = sum(1 for r in assign_rows
                  if r.get("assignment_basis") == "projection_identity")
    n_sim = sum(1 for r in assign_rows
                if r.get("assignment_basis") == "similarity_matching")
    n_off = sum(1 for r in assign_rows if r["assignment_status"] == "ASSIGNED"
                and r.get("pool") == "OFF_REGION")
    n_override = sum(1 for r in assign_rows
                     if "similarity_would_have_chosen" in (r.get("note") or ""))
    n_retro = sum(1 for r in assign_rows if r.get("is_retro") == "1"
                  and r["assignment_status"] == "ASSIGNED")
    log(f"[{tag}] assigned {n_assigned}/{len(assign_rows)} human-transcript "
        f"slots: {n_ident} by projection identity, {n_sim} by similarity "
        f"matching; {n_off} from outside the locus, {n_retro} retro/processed")
    log(f"[{tag}] {n_override} identity assignments where similarity would "
        f"have chosen differently (recorded, not acted on)")
    _sanity_check(assign_rows, outdir, tag)
    return out_path


def _sanity_check(rows: List[dict], outdir: Path, scope: str) -> None:
    """One bat candidate must never serve two human transcripts."""
    seen: Dict[Tuple[str, str, str], List[str]] = defaultdict(list)
    for r in rows:
        if r["assignment_status"] == "ASSIGNED":
            seen[(r["species"], r["gene"], r["candidate_id"])].append(r["human_transcript"])
    violations = [{"species": k[0], "gene": k[1], "candidate_id": k[2],
                   "human_transcripts": ",".join(v)}
                  for k, v in seen.items() if len(v) > 1]
    write_tsv(outdir / f"one_to_one_violations_{scope}.tsv", violations,
              ["species", "gene", "candidate_id", "human_transcripts"])
    if violations:
        log(f"WARNING: [{scope}] {len(violations)} one-to-one violations")
    else:
        log(f"[{scope}] one-to-one constraint holds for every species/gene group")
