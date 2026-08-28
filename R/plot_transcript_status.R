#!/usr/bin/env Rscript
# Transcript status across the phylogeny.
#
# One column per human transcript, grouped by gene; one row per species, ordered
# by the species tree. Cell colour is the status of the bat CDS assigned to that
# human transcript.
#
# Transcripts recovered in only a handful of species are mostly grey and crowd
# out the ones carrying signal, so a column is kept only when it was recovered
# in at least --min-species species. Set --min-species 0 to draw everything.
#
#   Rscript plot_transcript_status.R --status transcript_status.tsv \
#           --tree bat1k_tree.nwk --outdir results/07_plots
#   ... --min-species 0        # keep every transcript column

suppressPackageStartupMessages(library(ape))
args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "tree_layout.R"))

status_file <- get_arg(args, "--status")
tree_file   <- get_arg(args, "--tree")
outdir      <- get_arg(args, "--outdir")
width       <- as.numeric(get_arg(args, "--width", "22"))
height      <- as.numeric(get_arg(args, "--height", "24"))
compress    <- as.numeric(get_arg(args, "--tree-compress", "0.15"))
stem        <- get_arg(args, "--stem", "transcript_status_phylogeny")
min_species <- as.numeric(get_arg(args, "--min-species", "50"))

if (is.na(status_file) || is.na(tree_file) || is.na(outdir))
  stop("usage: --status FILE --tree FILE --outdir DIR")

STATUS_LEVELS <- c("complete", "partial", "fragmented", "pseudogenized", "not_found")
STATUS_COLOURS <- c(complete = "#1b9e77", partial = "#d95f02",
                    fragmented = "#7570b3", pseudogenized = "#e7298a",
                    not_found = "grey85")
# When a species somehow has two rows for one transcript, keep the worst.
SEVERITY <- c(not_found = 0, complete = 1, fragmented = 2, partial = 3,
              pseudogenized = 4)

dat <- read_table_tsv(status_file)
for (col in c("gene", "human_transcript", "species", "status"))
  if (!col %in% names(dat)) stop("status table lacks column: ", col)

dat$species_key <- norm_species(dat$species)
dat$transcript <- sub("\\..*$", "", dat$human_transcript)
dat$status[!dat$status %in% STATUS_LEVELS] <- "not_found"

key <- paste(dat$gene, dat$transcript, dat$species_key, sep = "\r")
dat <- dat[order(key, -SEVERITY[dat$status]), ]
dat <- dat[!duplicated(paste(dat$gene, dat$transcript, dat$species_key, sep = "\r")), ]

# Drop sparse transcript columns. "Present" means a sequence was recovered at
# all, whatever its status; not_found rows are absences by construction.
present <- dat[dat$status != "not_found" & dat$species_key != "Homo_sapiens", ]
per_tx <- table(paste(present$gene, present$transcript, sep = "\r"))
keep_key <- names(per_tx)[per_tx >= min_species]
dropped <- length(unique(paste(dat$gene, dat$transcript, sep = "\r"))) - length(keep_key)
dat <- dat[paste(dat$gene, dat$transcript, sep = "\r") %in% keep_key, , drop = FALSE]
if (nrow(dat) == 0)
  stop("no transcript reaches --min-species ", min_species, "; lower it")
cat(sprintf("kept %d transcript columns present in >= %d species (dropped %d)\n",
            length(keep_key), min_species, dropped))

tree <- read.tree(tree_file)
layout <- tree_layout(tree, compress = compress)
tips <- layout$tips

tree_max_x <- max(tips$x)
species_label_x <- tree_max_x + 0.15

# lay the columns out once at offset 0 just to learn how many x-units they span
n_columns <- max(column_layout(dat$gene, dat$transcript, gap = 1.2)$columns$x) + 1
box_start <- tree_max_x + label_gap(tips$species, n_columns, tree_max_x,
                                    width, font_size_mm = 2.2)

cl <- column_layout(dat$gene, dat$transcript, gap = 1.2, offset = box_start)
columns <- cl$columns
groups <- cl$groups

# Full species x transcript grid, so a species with no model shows as not_found.
grid_df <- merge(
  data.frame(species_key = tips$species_key, stringsAsFactors = FALSE),
  columns[, c("gene", "column", "x")], by = NULL
)
names(grid_df)[names(grid_df) == "column"] <- "transcript"

observed <- dat[, c("gene", "transcript", "species_key", "status")]
plot_df <- merge(grid_df, observed,
                 by = c("gene", "transcript", "species_key"), all.x = TRUE)
plot_df$status[is.na(plot_df$status)] <- "not_found"
plot_df$status <- factor(plot_df$status, levels = STATUS_LEVELS)
plot_df <- merge(plot_df, tips[, c("species_key", "y")], by = "species_key")

y_min <- min(tips$y); y_max <- max(tips$y)

p <- ggplot() +
  geom_segment(data = layout$segments,
               aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.9) +
  geom_text(data = tips, aes(x = species_label_x, y = y, label = species),
            hjust = 0, size = 2.2) +
  geom_tile(data = plot_df, aes(x = x, y = y, fill = status),
            width = 0.85, height = 0.75, colour = "white", linewidth = 0.10,
            key_glyph = ggplot2::draw_key_rect) +
  geom_rect(data = groups,
            aes(xmin = start - 0.55, xmax = end + 0.55,
                ymin = y_min - 0.5, ymax = y_max + 0.5),
            fill = NA, colour = "black", linewidth = 0.45) +
  geom_text(data = columns, aes(x = x, y = y_max + 1.4, label = column),
            angle = 90, hjust = 0, size = 2.5) +
  geom_text(data = groups, aes(x = mid, y = y_min - 3.2, label = gene),
            angle = 90, hjust = 1, fontface = "bold", size = 4) +
  scale_fill_manual(values = STATUS_COLOURS, drop = FALSE,
                    limits = STATUS_LEVELS) +
  coord_cartesian(ylim = c(y_min - 6, y_max + 8), clip = "off") +
  labs(fill = "Transcript status") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1.5, "cm"),
        legend.spacing.x = unit(0.8, "cm"),
        plot.margin = margin(80, 20, 80, 20))

save_plot(p, outdir, stem, width, height)

out_long <- plot_df[order(plot_df$gene, plot_df$transcript, plot_df$species_key),
                    c("gene", "transcript", "species_key", "status")]
utils::write.table(out_long, file.path(outdir, paste0(stem, "_matrix.tsv")),
                   sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", file.path(outdir, paste0(stem, "_matrix.tsv")), "\n")
