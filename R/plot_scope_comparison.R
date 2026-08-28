#!/usr/bin/env Rscript
# Does the search scope change which transcript you recover?
#
# One column per human transcript, grouped by gene; one row per species. The
# colour says whether the region-restricted search (the gene's own syntenic
# locus) and the unrestricted search (the gene anywhere it occurs, paralogs still
# excluded) chose the same model.
#
# Cells are outlined where the unrestricted search recovers a complete ORF that
# the syntenic locus does not - the OAS1-in-Phyllostomus situation, where the
# animal makes the protein but not from the expected place.
#
#   Rscript plot_scope_comparison.R --comparison scope_comparison.tsv \
#           --tree bat1k_tree.nwk --outdir results/07_plots

suppressPackageStartupMessages(library(ape))
args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "tree_layout.R"))

comp_file   <- get_arg(args, "--comparison")
tree_file   <- get_arg(args, "--tree")
outdir      <- get_arg(args, "--outdir")
width       <- as.numeric(get_arg(args, "--width", "22"))
height      <- as.numeric(get_arg(args, "--height", "24"))
compress    <- as.numeric(get_arg(args, "--tree-compress", "0.15"))
stem        <- get_arg(args, "--stem", "transcript_scope_comparison")
min_species <- as.numeric(get_arg(args, "--min-species", "50"))
only_changed <- "--only-changed" %in% args

if (is.na(comp_file) || is.na(tree_file) || is.na(outdir))
  stop("usage: --comparison FILE --tree FILE --outdir DIR")

LEVELS <- c("SAME_MODEL", "DIFFERENT_LOCUS", "ONLY_UNRESTRICTED",
            "ONLY_REGION_RESTRICTED", "NEITHER")
COLOURS <- c(SAME_MODEL = "#bdbdbd",
             DIFFERENT_LOCUS = "#d94801",
             ONLY_UNRESTRICTED = "#6a51a3",
             ONLY_REGION_RESTRICTED = "#238b45",
             NEITHER = "grey93")
LABELS <- c(SAME_MODEL = "same model in both",
            DIFFERENT_LOCUS = "different model outside the locus",
            ONLY_UNRESTRICTED = "found only outside the locus",
            ONLY_REGION_RESTRICTED = "found only in the locus",
            NEITHER = "no model in either")

dat <- read_table_tsv(comp_file)
for (col in c("gene", "human_transcript", "species", "outcome"))
  if (!col %in% names(dat)) stop("comparison table lacks column: ", col)

dat$species_key <- norm_species(dat$species)
dat$transcript <- sub("\\..*$", "", dat$human_transcript)
dat$outcome[!dat$outcome %in% LEVELS] <- "NEITHER"
dat <- dat[dat$species_key != "Homo_sapiens", , drop = FALSE]

# Same column filter as the status figures, so the three plots line up.
found <- dat[dat$outcome != "NEITHER", ]
per_tx <- table(paste(found$gene, found$transcript, sep = "\r"))
keep <- names(per_tx)[per_tx >= min_species]
if (only_changed) {
  changed <- dat[dat$outcome %in% c("DIFFERENT_LOCUS", "ONLY_UNRESTRICTED",
                                    "ONLY_REGION_RESTRICTED"), ]
  keep <- intersect(keep, unique(paste(changed$gene, changed$transcript, sep = "\r")))
}
dat <- dat[paste(dat$gene, dat$transcript, sep = "\r") %in% keep, , drop = FALSE]
if (nrow(dat) == 0) stop("nothing left to plot at --min-species ", min_species)
cat(sprintf("kept %d transcript columns\n", length(keep)))

# where the unrestricted search rescues a complete ORF the locus has lost
dat$rescued <- (dat$outcome %in% c("DIFFERENT_LOCUS", "ONLY_UNRESTRICTED") &
                dat$unrestricted_cds_status == "complete" &
                dat$region_cds_status != "complete")

tree <- read.tree(tree_file)
if ("Homo_sapiens" %in% tree$tip.label) tree <- drop.tip(tree, "Homo_sapiens")
layout <- tree_layout(tree, compress = compress)
tips <- layout$tips

tree_max_x <- max(tips$x)
species_label_x <- tree_max_x + 0.15
n_columns <- max(column_layout(dat$gene, dat$transcript, gap = 1.2)$columns$x) + 1
box_start <- tree_max_x + label_gap(tips$species, n_columns, tree_max_x,
                                    width, font_size_mm = 2.2)
cl <- column_layout(dat$gene, dat$transcript, gap = 1.2, offset = box_start)
columns <- cl$columns; groups <- cl$groups

grid_df <- merge(data.frame(species_key = tips$species_key, stringsAsFactors = FALSE),
                 columns[, c("gene", "column", "x")], by = NULL)
names(grid_df)[names(grid_df) == "column"] <- "transcript"
plot_df <- merge(grid_df, dat[, c("gene", "transcript", "species_key", "outcome", "rescued")],
                 by = c("gene", "transcript", "species_key"), all.x = TRUE)
plot_df$outcome[is.na(plot_df$outcome)] <- "NEITHER"
plot_df$rescued[is.na(plot_df$rescued)] <- FALSE
plot_df$outcome <- factor(plot_df$outcome, levels = LEVELS)
plot_df <- merge(plot_df, tips[, c("species_key", "y")], by = "species_key")

y_min <- min(tips$y); y_max <- max(tips$y)
n_changed <- sum(plot_df$outcome %in% c("DIFFERENT_LOCUS", "ONLY_UNRESTRICTED",
                                        "ONLY_REGION_RESTRICTED"))
n_rescued <- sum(plot_df$rescued)

p <- ggplot() +
  geom_segment(data = layout$segments,
               aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.9) +
  geom_text(data = tips, aes(x = species_label_x, y = y, label = species),
            hjust = 0, size = 2.2) +
  geom_tile(data = plot_df, aes(x = x, y = y, fill = outcome),
            width = 0.85, height = 0.75, colour = "white", linewidth = 0.10,
            key_glyph = ggplot2::draw_key_rect) +
  geom_tile(data = plot_df[plot_df$rescued, ], aes(x = x, y = y),
            width = 0.85, height = 0.75, fill = NA, colour = "black",
            linewidth = 0.45) +
  geom_rect(data = groups,
            aes(xmin = start - 0.55, xmax = end + 0.55,
                ymin = y_min - 0.5, ymax = y_max + 0.5),
            fill = NA, colour = "black", linewidth = 0.45) +
  geom_text(data = columns, aes(x = x, y = y_max + 1.4, label = column),
            angle = 90, hjust = 0, size = 2.5) +
  geom_text(data = groups, aes(x = mid, y = y_min - 3.2, label = gene),
            angle = 90, hjust = 1, fontface = "bold", size = 4) +
  scale_fill_manual(values = COLOURS, labels = LABELS, drop = FALSE,
                    limits = LEVELS) +
  coord_cartesian(ylim = c(y_min - 6, y_max + 8), clip = "off") +
  labs(fill = "Region-restricted vs unrestricted search",
       caption = sprintf(paste("black outline: a complete ORF recovered only outside the gene's own locus",
                               "(%d cells).  %d of %d cells differ between scopes."),
                         n_rescued, n_changed, nrow(plot_df))) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 13),
        legend.key.size = unit(1.1, "cm"),
        plot.caption = element_text(size = 12, hjust = 0.5),
        plot.margin = margin(80, 20, 60, 20))

save_plot(p, outdir, stem, width, height)
utils::write.table(plot_df[order(plot_df$gene, plot_df$transcript, plot_df$species_key),
                           c("gene", "transcript", "species_key", "outcome", "rescued")],
                   file.path(outdir, paste0(stem, "_matrix.tsv")),
                   sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", file.path(outdir, paste0(stem, "_matrix.tsv")), "\n")
