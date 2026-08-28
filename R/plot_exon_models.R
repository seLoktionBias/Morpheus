#!/usr/bin/env Rscript
#
# plot_exon_models.R - UCSC-style exon models for one gene. OPTIONAL: the
# pipeline does not call it, so run it by hand for whichever gene you want.
#
# Exons keep their true genomic width; introns are compressed by one shared
# transform per gene, so a two-exon model separated by 40 kb and a fifteen-exon
# model fit the same panel and the exons stay the thing you can see. The
# transform is built from the union of every exon in the gene, so all
# transcripts share one x-axis and equivalent exons line up vertically.
#
# Every gene is drawn 5' to 3', left to right, whatever strand it is on, so the
# figures stay comparable; on the minus strand the axis therefore counts down.
#
# Only the longest/anchor model carries the full exon1..exonN labels. On the
# other rows the shared exons are identified by colour and by lining up under
# the reference, so repeating the names on every row adds clutter and nothing
# else; those rows label only their novel exons.
#
# It reads either table, detecting the columns:
#
#   human   results/01_human_reference/human_isoform_exons.tsv
#   query   results/02_bat_search/all_species_candidate_exons.tsv
#
# Human models for one gene:
#   Rscript pipeline/lib/R/plot_exon_models.R --gene MX1
#
# Query models for one species and gene:
#   Rscript pipeline/lib/R/plot_exon_models.R \
#     --plot-tsv results/02_bat_search/all_species_candidate_exons.tsv \
#     --gene OAS1 --species Phyllostomus_discolor
#
# Options:
#   --gene NAME           gene to draw (required unless --all-genes)
#   --species NAME        restrict to one species (query tables)
#   --plot-tsv FILE       exon table; defaults to the human one
#   --context-tsv FILE    for the neighbouring-gene header
#   --outdir DIR          default results/07_plots/exon_models
#   --format pdf|png|both default both
#   --intron-mode         log | sqrt | linear | none   (default log)
#   --intron-factor N     intron compression strength  (default 45)
#   --min-exon-width N    smallest drawn exon, in transformed units (default 8)
#   --no-exon-labels      drop the per-exon labels
#   --label-all-rows      repeat exonN labels on every transcript, not just the reference
#   --no-flip-minus       keep minus-strand genes in ascending genomic order
#   --no-collapse-identical   keep transcripts whose models are identical
#   --all-genes           loop over every gene in the table

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NA) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i[length(i)] == length(args)) stop("missing value after ", flag)
  args[i[length(i)] + 1]
}
has_flag <- function(flag) any(args == flag)

plot_tsv       <- get_arg("--plot-tsv", "results/01_human_reference/human_isoform_exons.tsv")
context_tsv    <- get_arg("--context-tsv", "results/01_human_reference/gene_context.tsv")
gene_filter    <- get_arg("--gene", NA)
species_filter <- get_arg("--species", NA)
outdir         <- get_arg("--outdir", "results/07_plots/exon_models")
fmt            <- get_arg("--format", "both")
intron_mode    <- get_arg("--intron-mode", "log")
intron_factor  <- as.numeric(get_arg("--intron-factor", "45"))
min_exon_width <- as.numeric(get_arg("--min-exon-width", "8"))
label_cex      <- as.numeric(get_arg("--label-cex", "0.58"))
tx_cex         <- as.numeric(get_arg("--tx-cex", "0.70"))
main_cex       <- as.numeric(get_arg("--main-cex", "1.08"))
arrow_count    <- as.integer(get_arg("--arrow-count", "7"))
width          <- as.numeric(get_arg("--width", "13"))
height         <- as.numeric(get_arg("--height", NA))

show_exon_labels   <- !has_flag("--no-exon-labels")
label_all_rows     <- has_flag("--label-all-rows")
flip_minus         <- !has_flag("--no-flip-minus")
collapse_identical <- !has_flag("--no-collapse-identical")
all_genes          <- has_flag("--all-genes")

if (is.na(gene_filter) && !all_genes)
  stop("give --gene NAME, or --all-genes to draw every gene in the table")
if (!file.exists(plot_tsv)) stop("cannot find --plot-tsv: ", plot_tsv)

read_tsv <- function(path)
  read.delim(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
             check.names = FALSE, fill = TRUE)

pick_col <- function(df, candidates, required = TRUE, what = "column") {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0) return(hit[1])
  if (required) stop("missing ", what, "; tried: ", paste(candidates, collapse = ", "))
  NA
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nchar(x) == 0, "plot", x)
}

# A candidate id carries species and gene for uniqueness, which makes a useless
# row label; the TOGA projection is what a reader wants to see.
clean_tx_label <- function(x) {
  x <- as.character(x)
  x <- sub("#retro$", " (retro)", x)
  x
}

# --------------------------------------------------------------------------
# intron-compressed coordinate transform
# --------------------------------------------------------------------------

compress_gap <- function(gap, mode = "log", factor = 45) {
  gap <- pmax(0, gap)
  switch(mode,
         none = gap,
         linear = gap / factor,
         sqrt = sqrt(gap),
         log10(gap + 1) * factor)
}

interval_union <- function(starts, ends) {
  o <- order(starts, ends)
  starts <- starts[o]; ends <- ends[o]
  us <- numeric(0); ue <- numeric(0)
  for (i in seq_along(starts)) {
    if (length(us) == 0 || starts[i] > ue[length(ue)] + 1) {
      us <- c(us, starts[i]); ue <- c(ue, ends[i])
    } else {
      ue[length(ue)] <- max(ue[length(ue)], ends[i])
    }
  }
  data.frame(start = us, end = ue)
}

# One transform per gene, built from the union of all exons, so every
# transcript is drawn on the same axis and shared exons align down the page.
make_transform <- function(starts, ends, mode = "log", factor = 45) {
  u <- interval_union(starts, ends)
  x0 <- numeric(nrow(u)); cursor <- 0
  for (i in seq_len(nrow(u))) {
    if (i == 1) {
      x0[i] <- 0
    } else {
      prev_len <- u$end[i - 1] - u$start[i - 1] + 1
      gap <- u$start[i] - u$end[i - 1] - 1
      cursor <- cursor + prev_len + compress_gap(gap, mode, factor)
      x0[i] <- cursor
    }
  }
  function(pos) {
    pos <- as.numeric(pos)
    vapply(pos, function(p) {
      idx <- which(p >= u$start & p <= u$end)
      if (length(idx) > 0) {
        i <- idx[1]
        x0[i] + (p - u$start[i])
      } else if (p < u$start[1]) {
        -(u$start[1] - p)
      } else if (p > u$end[nrow(u)]) {
        last <- nrow(u)
        x0[last] + (u$end[last] - u$start[last] + 1) +
          compress_gap(p - u$end[last], mode, factor)
      } else {
        left <- max(which(u$end < p)); right <- left + 1
        gap_start <- u$end[left] + 1
        gap_end <- u$start[right] - 1
        gap_len <- max(1, gap_end - gap_start + 1)
        x0[left] + (u$end[left] - u$start[left] + 1) +
          ((p - gap_start) / gap_len) * compress_gap(gap_len, mode, factor)
      }
    }, numeric(1))
  }
}

# --------------------------------------------------------------------------
# palette
# --------------------------------------------------------------------------

EXON_FILL <- c(longest = "#4C78A8",       # shared with the longest isoform
               human_novel = "#F58518",   # isoform-specific human exon
               bat_novel = "#54A24B",     # query-specific exon
               other = "#9E9E9E")

exon_colour <- function(type, label) {
  type <- tolower(as.character(type)); label <- as.character(label)
  ifelse(type %in% c("longest", "longest_isoform", "reference", "anchor") |
           grepl("^exon[0-9]+$", label), EXON_FILL[["longest"]],
    ifelse(grepl("^human_novel", label) | type %in% c("human_novel", "isoform_specific"),
           EXON_FILL[["human_novel"]],
      ifelse(grepl("^bat_novel", label) | type == "bat_novel",
             EXON_FILL[["bat_novel"]], EXON_FILL[["other"]])))
}

# --------------------------------------------------------------------------
# load
# --------------------------------------------------------------------------

x <- read_tsv(plot_tsv)
if (nrow(x) == 0) stop("no rows in ", plot_tsv)

gene_col   <- pick_col(x, c("gene", "target_gene", "gene_name"), what = "gene column")
tx_col     <- pick_col(x, c("projection", "transcript_id", "candidate_id"), what = "transcript column")
uid_col    <- pick_col(x, c("candidate_id", "transcript_id"), required = FALSE)
chrom_col  <- pick_col(x, c("chrom", "chr", "scaffold"), required = FALSE)
start_col  <- pick_col(x, c("start", "exon_start"), what = "start column")
end_col    <- pick_col(x, c("end", "exon_end"), what = "end column")
strand_col <- pick_col(x, c("strand"), required = FALSE)
rank_col   <- pick_col(x, c("transcript_exon_rank", "exon_rank"), required = FALSE)
label_col  <- pick_col(x, c("exon_label", "label"), what = "exon label column")
type_col   <- pick_col(x, c("exon_label_type", "exon_type"), required = FALSE)
sp_col     <- pick_col(x, c("species"), required = FALSE)
anchor_col <- pick_col(x, c("is_anchor", "is_longest"), required = FALSE)
pool_col   <- pick_col(x, c("pool"), required = FALSE)
prev_col   <- pick_col(x, c("previous_gene", "upstream_gene"), required = FALSE)
next_col   <- pick_col(x, c("next_gene", "downstream_gene"), required = FALSE)

d <- data.frame(
  gene    = x[[gene_col]],
  tx_id   = if (!is.na(uid_col)) x[[uid_col]] else x[[tx_col]],
  tx_lab  = clean_tx_label(x[[tx_col]]),
  chrom   = if (!is.na(chrom_col)) x[[chrom_col]] else "NA",
  start   = suppressWarnings(as.numeric(x[[start_col]])),
  end     = suppressWarnings(as.numeric(x[[end_col]])),
  strand  = if (!is.na(strand_col)) x[[strand_col]] else "+",
  rank    = if (!is.na(rank_col)) suppressWarnings(as.numeric(x[[rank_col]])) else NA,
  label   = x[[label_col]],
  type    = if (!is.na(type_col)) x[[type_col]] else "other",
  species = if (!is.na(sp_col)) x[[sp_col]] else "Homo_sapiens",
  anchor  = if (!is.na(anchor_col)) suppressWarnings(as.numeric(x[[anchor_col]])) else 0,
  pool    = if (!is.na(pool_col)) x[[pool_col]] else NA,
  prev_g  = if (!is.na(prev_col)) x[[prev_col]] else NA,
  next_g  = if (!is.na(next_col)) x[[next_col]] else NA,
  stringsAsFactors = FALSE
)
d <- d[!is.na(d$start) & !is.na(d$end) & !is.na(d$gene) & d$gene != "", ]
lo <- pmin(d$start, d$end); d$end <- pmax(d$start, d$end); d$start <- lo

# The query table holds every pool; only the gene's own locus belongs in a
# gene-model figure, otherwise unrelated copies pile onto one axis.
if (!all(is.na(d$pool))) d <- d[is.na(d$pool) | d$pool == "IN_REGION", ]
if (!is.na(species_filter)) {
  d <- d[d$species == species_filter, ]
  if (nrow(d) == 0) stop("no rows for --species ", species_filter)
}

ctx <- if (!is.na(context_tsv) && file.exists(context_tsv))
  tryCatch(read_tsv(context_tsv), error = function(e) NULL) else NULL

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# draw one gene
# --------------------------------------------------------------------------

plot_gene <- function(g) {
  gx <- d[d$gene == g, , drop = FALSE]
  if (nrow(gx) == 0) return(character(0))

  txs <- unique(gx$tx_id)
  key <- setNames(character(length(txs)), txs)
  lab <- setNames(character(length(txs)), txs)
  anc <- setNames(logical(length(txs)), txs)
  for (t in txs) {
    z <- gx[gx$tx_id == t, , drop = FALSE]
    z <- z[order(z$rank, z$start, z$end), ]
    key[t] <- paste(paste(z$start, z$end, z$label, sep = ":"), collapse = "|")
    lab[t] <- z$tx_lab[1]
    anc[t] <- any(z$anchor == 1)
  }

  keep <- txs
  if (collapse_identical) keep <- txs[!duplicated(key[txs])]
  gx <- gx[gx$tx_id %in% keep, , drop = FALSE]
  tx_order <- keep[order(!anc[keep], lab[keep])]

  species <- unique(gx$species)[1]
  chrom <- unique(gx$chrom)[1]
  strand <- unique(gx$strand)[1]
  if (is.na(strand) || strand == "") strand <- "+"

  trans0 <- make_transform(gx$start, gx$end, intron_mode, intron_factor)
  # Draw every gene in reading order. Negating the transform mirrors the panel,
  # so a minus-strand gene starts at exon1 on the left like every other figure;
  # the axis labels stay true genomic coordinates and simply descend.
  reverse_axis <- flip_minus && identical(strand, "-")
  trans <- if (reverse_axis) function(p) -trans0(p) else trans0

  if (reverse_axis) {
    gx$x1 <- trans(gx$end); gx$x2 <- trans(gx$start)
  } else {
    gx$x1 <- trans(gx$start); gx$x2 <- trans(gx$end)
  }

  # a 5 bp exon is invisible next to a 1.5 kb one; give every exon a floor
  narrow <- (gx$x2 - gx$x1) < min_exon_width
  mid <- (gx$x1 + gx$x2) / 2
  gx$x1[narrow] <- mid[narrow] - min_exon_width / 2
  gx$x2[narrow] <- mid[narrow] + min_exon_width / 2

  xlim <- range(c(gx$x1, gx$x2), finite = TRUE)
  pad <- diff(xlim) * 0.08; if (!is.finite(pad) || pad == 0) pad <- 100
  xlim <- c(xlim[1] - pad, xlim[2] + pad)

  n_tx <- length(tx_order)
  ylim <- c(0.2, n_tx + 2.2)

  # The reference row keeps the full exon vocabulary; the rest label only what
  # is new, since colour and vertical alignment already say the rest.
  reference_tx <- if (any(anc[tx_order])) tx_order[which(anc[tx_order])[1]] else tx_order[1]

  prev_g <- gx$prev_g[1]; next_g <- gx$next_g[1]
  if (!is.null(ctx)) {
    cg <- pick_col(ctx, c("gene", "target_gene"), required = FALSE)
    if (!is.na(cg)) {
      cc <- ctx[ctx[[cg]] == g, , drop = FALSE]
      # A query context table holds every species, so the gene alone is not a
      # key; without the species filter the header shows another animal's
      # neighbours.
      csp <- pick_col(cc, c("species"), required = FALSE)
      if (!is.na(csp)) cc <- cc[cc[[csp]] == species, , drop = FALSE]
      # and every locus of that gene, so prefer the home locus row
      chome <- pick_col(cc, c("is_home_locus"), required = FALSE)
      if (!is.na(chome) && any(cc[[chome]] == 1))
        cc <- cc[cc[[chome]] == 1, , drop = FALSE]
      cp <- pick_col(cc, c("upstream_gene", "previous_gene"), required = FALSE)
      cn <- pick_col(cc, c("downstream_gene", "next_gene"), required = FALSE)
      if (nrow(cc) > 0 && !is.na(cp)) prev_g <- cc[[cp]][1]
      if (nrow(cc) > 0 && !is.na(cn)) next_g <- cc[[cn]][1]
    }
  }
  if (is.na(prev_g) || prev_g == "") prev_g <- "NA"
  if (is.na(next_g) || next_g == "") next_g <- "NA"

  h <- if (!is.na(height)) height else max(3.6, 0.42 * n_tx + 3.1)
  stem <- file.path(outdir, safe_name(paste(species, g, sep = "__")))
  devices <- if (fmt == "both") c("pdf", "png") else fmt
  written <- character(0)

  for (dev_fmt in devices) {
    if (dev_fmt == "pdf") pdf(paste0(stem, ".pdf"), width = width, height = h, onefile = FALSE)
    else png(paste0(stem, ".png"), width = width, height = h, units = "in", res = 300)

    op <- par(no.readonly = TRUE)
    par(mar = c(5.2, 12.5, 4.6, 2.0), xpd = NA, bg = "white")

    plot(NA, NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "",
         yaxt = "n", xaxt = "n", bty = "n",
         main = paste0(g, "  |  ", species), cex.main = main_cex)

    ticks <- pretty(range(c(gx$start, gx$end)), n = 5)
    ticks <- ticks[ticks >= min(gx$start) & ticks <= max(gx$end)]
    axis(1, at = trans(ticks), labels = format(ticks, big.mark = ",", trim = TRUE),
         cex.axis = 0.72, col = "grey40", col.axis = "grey25")
    mtext(sprintf("%s  |  strand %s  |  introns compressed (%s x%g); exon widths are real%s",
                  chrom, strand, intron_mode, intron_factor,
                  if (reverse_axis) "  |  drawn 5' to 3', so coordinates descend" else ""),
          side = 3, line = 0.55, cex = 0.72, col = "grey30")

    # neighbouring genes, so the locus can be read at a glance
    text(xlim[1], ylim[2] - 0.42, prev_g, adj = c(0, 0.5), cex = 0.82, col = "grey30")
    text(mean(xlim), ylim[2] - 0.42, g, adj = c(0.5, 0.5), cex = 0.95, font = 2)
    text(xlim[2], ylim[2] - 0.42, next_g, adj = c(1, 0.5), cex = 0.82, col = "grey30")

    ymap <- setNames(rev(seq_len(n_tx)), tx_order)
    rect_h <- 0.23

    for (t in tx_order) {
      z <- gx[gx$tx_id == t, , drop = FALSE]
      z <- z[order(z$rank, z$start, z$end), ]
      y <- ymap[t]

      x_lo <- min(z$x1); x_hi <- max(z$x2)
      segments(x_lo, y, x_hi, y, lty = 2, lwd = 1, col = "grey45")

      if (arrow_count > 0 && x_hi > x_lo) {
        xs <- seq(x_lo, x_hi, length.out = arrow_count + 2)
        xs <- xs[-c(1, length(xs))]
        dx <- (x_hi - x_lo) * 0.015; if (!is.finite(dx) || dx <= 0) dx <- 5
        # the panel is already in reading order, so transcription runs rightwards
        # on both strands
        for (a in xs)
          arrows(a - dx, y, a + dx, y, length = 0.045, lwd = 0.8, col = "grey45")
      }

      rect(z$x1, y - rect_h, z$x2, y + rect_h,
           col = exon_colour(z$type, z$label), border = "black", lwd = 0.45)

      # horizontal labels, staggered so neighbours cannot collide
      if (show_exon_labels) {
        novel <- grepl("^(human|bat)_novel", z$label)
        keep <- if (label_all_rows || identical(t, reference_tx)) rep(TRUE, nrow(z)) else novel
        if (any(keep)) {
          zz <- z[keep, , drop = FALSE]
          offs <- rep(c(0.16, 0.34), length.out = nrow(zz))
          text((zz$x1 + zz$x2) / 2, y + rect_h + offs, labels = zz$label,
               cex = label_cex, adj = c(0.5, 0), col = "grey15")
        }
      }

      text(xlim[1], y, labels = lab[t], adj = c(1.04, 0.5), cex = tx_cex,
           font = if (anc[t]) 2 else 1)
    }

    # A fractional inset scales with the panel, so on a tall gene the legend
    # walks off the page. Offset a fixed distance in inches instead.
    legend_inset <- 0.66 / par("pin")[2]
    legend("bottom", horiz = TRUE, bty = "n", inset = c(0, -legend_inset), cex = 0.74,
           fill = EXON_FILL[c("longest", "human_novel", "bat_novel")],
           border = "black",
           legend = c("shared with longest isoform", "isoform-specific exon",
                      "query-specific exon"))

    par(op); dev.off()
    written <- c(written, paste0(stem, ".", dev_fmt))
  }
  for (f in written) cat("wrote", f, "\n")
  written
}

genes <- if (all_genes) sort(unique(d$gene)) else gene_filter
missing <- setdiff(genes, unique(d$gene))
if (length(missing))
  stop("no exons for gene(s) ", paste(missing, collapse = ", "),
       ". Available: ", paste(sort(unique(d$gene)), collapse = ", "))

invisible(lapply(genes, plot_gene))
