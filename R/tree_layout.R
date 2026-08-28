# Shared helpers for the phylogeny-plus-matrix figures.
#
# Deliberately base R + ape + ggplot2 only: no tidyverse dependency, so the
# plots run in any environment that can already read a tree.

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(grid)
})

get_arg <- function(args, flag, default = NA_character_) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i[length(i)] == length(args)) stop("missing value after ", flag)
  args[i[length(i)] + 1]
}

norm_species <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

read_table_tsv <- function(path) {
  if (!file.exists(path)) stop("missing input table: ", path)
  utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                    quote = "", na.strings = c("", "NA"), check.names = FALSE)
}

# Lay a cladogram out in ggplot coordinates and return the pieces needed to
# draw it beside a matrix: the branch segments and the tip positions.
tree_layout <- function(tree, compress = 0.15) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot.phylo(tree, plot = FALSE, type = "cladogram",
             use.edge.length = FALSE, direction = "rightwards")
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)

  edges <- as.data.frame(tree$edge)
  names(edges) <- c("parent", "child")
  edges$x_parent <- pp$xx[edges$parent] * compress
  edges$x_child <- pp$xx[edges$child] * compress
  edges$y_parent <- pp$yy[edges$parent]
  edges$y_child <- pp$yy[edges$child]

  segments <- rbind(
    data.frame(x = edges$x_parent, xend = edges$x_child,
               y = edges$y_child, yend = edges$y_child),
    data.frame(x = edges$x_parent, xend = edges$x_parent,
               y = edges$y_parent, yend = edges$y_child)
  )

  n_tips <- length(tree$tip.label)
  tips <- data.frame(
    species = tree$tip.label,
    species_key = norm_species(tree$tip.label),
    x = pp$xx[seq_len(n_tips)] * compress,
    y = pp$yy[seq_len(n_tips)],
    stringsAsFactors = FALSE
  )
  list(segments = segments, tips = tips)
}

# Column positions for a matrix whose columns are grouped by gene, with a gap
# between groups.
column_layout <- function(gene, column, gap = 1.2, offset = 0) {
  keys <- unique(data.frame(gene = gene, column = column, stringsAsFactors = FALSE))
  keys <- keys[order(keys$gene, keys$column), , drop = FALSE]
  genes <- unique(keys$gene)

  cols <- do.call(rbind, lapply(seq_along(genes), function(i) {
    sub <- keys[keys$gene == genes[i], , drop = FALSE]
    sub$index_within_gene <- seq_len(nrow(sub))
    sub
  }))

  n_per_gene <- table(factor(cols$gene, levels = genes))
  starts <- cumsum(c(0, head(as.numeric(n_per_gene), -1))) + seq_along(genes) * gap
  names(starts) <- genes

  cols$x <- starts[cols$gene] + cols$index_within_gene - 1 + offset
  groups <- data.frame(
    gene = genes,
    start = as.numeric(starts) + offset,
    end = as.numeric(starts) + as.numeric(n_per_gene) - 1 + offset,
    stringsAsFactors = FALSE
  )
  groups$mid <- (groups$start + groups$end) / 2
  list(columns = cols, groups = groups)
}

# Space to leave between the tip labels and the matrix, in x-units.
#
# A fixed gap breaks as soon as the number of matrix columns changes: the same
# 4.4 units that clear the labels on a 21-column figure put them underneath the
# tiles on a 135-column one. So solve for the gap that makes the longest label
# fit, given the page width and how many columns share it.
label_gap <- function(tip_labels, n_columns, tree_span, width_inches,
                      font_size_mm, minimum = 2) {
  label_pt <- max(nchar(tip_labels)) * font_size_mm * 2.845 * 0.55
  page_pt <- width_inches * 72
  if (label_pt >= page_pt) return(minimum)
  gap <- label_pt * (tree_span + n_columns) / (page_pt - label_pt)
  max(gap, minimum)
}

save_plot <- function(p, outdir, stem, width, height) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  # theme_void() leaves the canvas transparent, which renders as black text on
  # black in most viewers. Paint it white explicitly.
  p <- p + theme(plot.background = element_rect(fill = "white", colour = NA),
                 panel.background = element_rect(fill = "white", colour = NA))
  # Outline every legend key. Pale fills - the near-white middle of a diverging
  # ramp, or a grey "absent" category - are invisible against a white legend
  # background and read as a missing entry rather than a colour. The key
  # background is drawn with a border as well as the glyph, because theme_void()
  # blanks legend.key and override.aes alone does not survive a layer that fixes
  # its own colour.
  p <- p + guides(fill = guide_legend(
                    override.aes = list(colour = "grey25", linewidth = 0.35)),
                  shape = "none") +
    theme(legend.key = element_rect(fill = NA, colour = "grey25", linewidth = 0.35),
          legend.key.spacing.x = unit(0.25, "cm"))
  pdf_path <- file.path(outdir, paste0(stem, ".pdf"))
  png_path <- file.path(outdir, paste0(stem, ".png"))
  ggsave(pdf_path, p, width = width, height = height, limitsize = FALSE)
  ggsave(png_path, p, width = width, height = height, dpi = 300, limitsize = FALSE)
  cat("wrote", pdf_path, "\n")
  cat("wrote", png_path, "\n")
}
