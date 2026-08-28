#!/usr/bin/env Rscript
# Gene copy number across the phylogeny.
#
# One column per gene, one row per species ordered by the species tree. Cell
# colour is the copy number; the count is printed whenever it differs from one,
# so gains and losses read at a glance.
#
# The default metric is total_copies, because that is what copy number variation
# means: every distinct genomic locus that really belongs to the gene. A copy
# whose ORF is disrupted, uncertain, or retro/processed is still a copy, so
# rather than dropping it from the count the cell is marked with a dot. Use
# --metric functional_copies for the stricter count.
#
#   Rscript plot_copy_number.R --matrix copy_number_matrix.tsv \
#           --tree bat1k_tree.nwk --outdir results/07_plots

suppressPackageStartupMessages(library(ape))
args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "tree_layout.R"))

matrix_file <- get_arg(args, "--matrix")
tree_file   <- get_arg(args, "--tree")
outdir      <- get_arg(args, "--outdir")
width       <- as.numeric(get_arg(args, "--width", "16"))
height      <- as.numeric(get_arg(args, "--height", "24"))
compress    <- as.numeric(get_arg(args, "--tree-compress", "0.15"))
stem        <- get_arg(args, "--stem", "gene_copy_number_phylogeny")
metric      <- get_arg(args, "--metric", "total_copies")
METRICS <- c("total_copies", "copies_excluding_retro", "functional_copies")
if (!metric %in% METRICS)
  stop("--metric must be one of: ", paste(METRICS, collapse = ", "))

if (is.na(matrix_file) || is.na(tree_file) || is.na(outdir))
  stop("usage: --matrix FILE --tree FILE --outdir DIR")

# One colour per copy number, no lumped top category.
#
# Up to 12 levels the reference palette is used step for step. Beyond that,
# interpolating the reference directly is a poor choice: its middle steps are
# almost white, so a long scale spends a third of its range invisible against
# the page. The ramp is instead rebuilt from blue to dark orange through the
# reference's own saturated anchors, dropping the near-white steps, which keeps
# every level distinguishable however far the counts run.
count_colours <- function(levels) {
  reference <- if (requireNamespace("paletteer", quietly = TRUE)) {
    as.character(paletteer::paletteer_d("colorBlindness::Blue2DarkOrange12Steps"))
  } else {
    c("#1E8E99", "#51C3CC", "#99F9FF", "#B2FCFF", "#CCFEFF", "#E5FFFF",
      "#FFE5CC", "#FFCA99", "#FFAD65", "#FF8E32", "#CC5800", "#993F00")
  }
  n <- length(levels)
  if (n <= length(reference)) {
    return(stats::setNames(reference[seq_len(n)], levels))
  }
  # luminance of each reference colour; drop the washed-out middle before
  # interpolating so the extended ramp stays legible end to end
  rgb_m <- grDevices::col2rgb(reference)
  luma <- (0.2126 * rgb_m[1, ] + 0.7152 * rgb_m[2, ] + 0.0722 * rgb_m[3, ]) / 255
  anchors <- reference[luma < 0.92]
  if (length(anchors) < 2) anchors <- c("#1E8E99", "#993F00")
  stats::setNames(grDevices::colorRampPalette(anchors)(n), levels)
}

dat <- read_table_tsv(matrix_file)
for (col in c("species", "gene", metric))
  if (!col %in% names(dat)) stop("copy number matrix lacks column: ", col)

dat$species_key <- norm_species(dat$species)
dat$count <- as.integer(dat[[metric]])
dat$functional_copies <- as.integer(dat$functional_copies)
dat$total_copies <- as.integer(dat$total_copies)
if (!"copies_excluding_retro" %in% names(dat))
  dat$copies_excluding_retro <- dat$total_copies
# a copy that is retro/processed, disrupted, or uncertain is still a copy, but
# it is worth flagging on the figure
dat$n_not_functional <- pmax(0L, dat$total_copies - dat$functional_copies)

tree <- read.tree(tree_file)
# Human is the reference, not a surveyed genome: drop it from a bat copy-number
# figure unless it was actually counted.
if (!"Homo_sapiens" %in% dat$species_key && "Homo_sapiens" %in% tree$tip.label)
  tree <- drop.tip(tree, "Homo_sapiens")

layout <- tree_layout(tree, compress = compress)
tips <- layout$tips

tree_max_x <- max(tips$x)
species_label_x <- tree_max_x + 0.15

genes <- sort(unique(dat$gene))
n_columns <- max(column_layout(genes, genes, gap = 0.35)$columns$x) + 1
box_start <- tree_max_x + label_gap(tips$species, n_columns, tree_max_x,
                                    width, font_size_mm = 2.2)
cl <- column_layout(genes, genes, gap = 0.35, offset = box_start)
columns <- cl$columns

grid_df <- expand.grid(species_key = tips$species_key, gene = genes,
                       stringsAsFactors = FALSE)
plot_df <- merge(grid_df,
                 dat[, c("species_key", "gene", "count", "total_copies",
                         "functional_copies", "n_not_functional")],
                 by = c("species_key", "gene"), all.x = TRUE)
for (col in c("count", "total_copies", "functional_copies", "n_not_functional"))
  plot_df[[col]][is.na(plot_df[[col]])] <- 0L
COUNT_LEVELS <- as.character(seq(0, max(plot_df$count, na.rm = TRUE)))
COUNT_COLOURS <- count_colours(COUNT_LEVELS)
plot_df$bin <- factor(as.character(plot_df$count), levels = COUNT_LEVELS)
plot_df <- merge(plot_df, columns[, c("column", "x")],
                 by.x = "gene", by.y = "column")
plot_df <- merge(plot_df, tips[, c("species_key", "y")], by = "species_key")

labels_df <- plot_df[plot_df$count != 1, ]
flag_df <- plot_df[plot_df$n_not_functional > 0, ]

y_min <- min(tips$y); y_max <- max(tips$y)
metric_label <- switch(metric,
  total_copies = "Gene copies",
  copies_excluding_retro = "Gene copies (excluding retro/processed)",
  functional_copies = "Functional gene copies")

p <- ggplot() +
  geom_segment(data = layout$segments,
               aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.9) +
  geom_text(data = tips, aes(x = species_label_x, y = y, label = species),
            hjust = 0, size = 2.2) +
  geom_tile(data = plot_df, aes(x = x, y = y, fill = bin),
            width = 0.88, height = 0.8, colour = "white", linewidth = 0.12,
            show.legend = TRUE, key_glyph = ggplot2::draw_key_rect) +
  geom_text(data = labels_df, aes(x = x, y = y, label = count),
            size = 2.1, colour = "grey10") +
  geom_point(data = flag_df, aes(x = x + 0.34, y = y + 0.29),
             shape = 16, size = 0.55, colour = "grey20") +
  geom_text(data = columns, aes(x = x, y = y_max + 1.2, label = column),
            angle = 90, hjust = 0, size = 3, fontface = "bold") +
  scale_fill_manual(values = COUNT_COLOURS, drop = FALSE,
                    limits = COUNT_LEVELS) +
  coord_cartesian(ylim = c(y_min - 3, y_max + 10), clip = "off") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(fill = metric_label,
       caption = paste0("dot = at least one copy at that locus is retro/processed, ",
                        "ORF-disrupted, or called uncertain by TOGA2")) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.1, "cm"),
        plot.caption = element_text(size = 11, hjust = 0.5),
        plot.margin = margin(60, 20, 40, 20))

save_plot(p, outdir, stem, width, height)

out_long <- plot_df[order(plot_df$gene, plot_df$species_key),
                    c("gene", "species_key", "count", "total_copies",
                      "functional_copies", "n_not_functional")]
utils::write.table(out_long, file.path(outdir, paste0(stem, "_matrix.tsv")),
                   sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", file.path(outdir, paste0(stem, "_matrix.tsv")), "\n")
