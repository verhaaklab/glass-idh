# GLASS-I: CN Segments Astrocytomas Recurrence
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# Load libraries
library(circlize)
library(tidyverse)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(grid)
library(ComplexHeatmap)

# ── Load cytoband data and define chromosome arms ──
cyto <- read.cytoband(species = "hg19")$df
colnames(cyto) <- c("chrom", "start", "end", "band", "gieStain")
cyto$arm <- ifelse(grepl("^p", cyto$band), "p", "q")

cyto <- cyto %>%
  group_by(chr = chrom, arm) %>%
  summarise(start = min(start), end = max(end), .groups = "drop") %>%
  filter(chr %in% paste0("chr", 1:22))  # autosomes only

cyto$chr_arm <- paste0(cyto$chr, cyto$arm)

# ── Make 1Mb bins ──
genome_bins <- tileGenome(
  seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg38)[paste0("chr", 1:22)],
  tilewidth = 1e6,
  cut.last.tile.in.chrom = TRUE
)

# ── Segment processing function ──
process_segments <- function(file) {
  seg <- read.csv(file)
  seg$chrom <- ifelse(grepl("^chr", seg$chrom), seg$chrom, paste0("chr", seg$chrom))
  seg <- seg %>%
    filter(!is.na(log2_copy_ratio)) %>%
    mutate(log2_copy_ratio = pmax(pmin(log2_copy_ratio, 2), -2))
  
  gr_seg <- GRanges(seqnames = seg$chrom,
                    ranges = IRanges(start = seg$start, end = seg$end),
                    score = seg$log2_copy_ratio)
  
  ov <- findOverlaps(genome_bins, gr_seg)
  df <- data.frame(bin_id = queryHits(ov),
                   score = mcols(gr_seg)$score[subjectHits(ov)])
  
  df %>%
    group_by(bin_id) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")
}

# ── Load and process recurrence dataset ──
mean_scores <- process_segments("~/gistic_circos_input_noncodel_recurrence.csv")                         ## created with prepare_gistic_input.R

# ── Create bin-level dataframe ──
bin_coords <- as.data.frame(genome_bins)
recurrence_seg <- bin_coords[mean_scores$bin_id, c("seqnames", "start", "end")]
recurrence_seg$score <- mean_scores$mean_score
colnames(recurrence_seg)[1] <- "chrom"

# ── Assign arm color to each bin ──
gr_bins <- GRanges(seqnames = recurrence_seg$chrom,
                   ranges = IRanges(start = recurrence_seg$start, end = recurrence_seg$end))

gr_arms <- GRanges(seqnames = cyto$chr,
                   ranges = IRanges(start = cyto$start, end = cyto$end),
                   arm = cyto$chr_arm)

# Find overlaps and get arm labels
hits <- findOverlaps(gr_bins, gr_arms)
recurrence_seg$arm <- NA
recurrence_seg$arm[queryHits(hits)] <- mcols(gr_arms)$arm[subjectHits(hits)]

# Compute mean recurrence score per arm for coloring
df_arm <- recurrence_seg %>%
  filter(!is.na(arm)) %>%
  group_by(arm) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

# Join color info
recurrence_seg <- left_join(recurrence_seg, df_arm, by = "arm") %>%
  filter(!is.na(mean_score))

# ── Define color function ──
col_fun <- colorRamp2(
  c(-0.4, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.4),
  c("#2066a8", "#8ec1da", "#cde1ec", "white", "white", "white", "#f6d6c2", "#d47264", "#a00000")
)

# ── Plot ──
circos.clear()
circos.par("track.margin" = c(0.01, 0.01), cell.padding = c(0, 0, 0, 0))
circos.initializeWithIdeogram(species = "hg19", chromosome.index = paste0("chr", c(1:22)))

# Add blank track for layering polygons
circos.track(ylim = c(-0.8, 0.8), track.height = 0.4, panel.fun = function(x, y) {})

# Plot 1Mb bins as polygons, colored by arm-level mean recurrence score
for (i in seq_len(nrow(recurrence_seg))) {
  row <- recurrence_seg[i, ]
  chrom <- row$chrom
  start <- as.numeric(row$start)
  end   <- as.numeric(row$end)
  score <- as.numeric(row$score)
  arm_score <- as.numeric(row$mean_score)
  
  if (!is.na(start) && !is.na(end) && end > start) {
    circos.polygon(
      x = c(start, end, end, start),
      y = c(0, 0, score, score),
      sector.index = chrom,
      col = col_fun(arm_score),
      border = NA
    )
  }
}

# ── Add center label ──
grid.text("log2(CN ratio)\nRecurrence\nShape = Bin, Color = Arm",
          x = 0.5, y = 0.5,
          gp = gpar(col = "#800074", fontsize = 14, fontface = "bold"))

# ── Add legend ──
legend_log2 <- Legend(
  col_fun = col_fun,
  title = "log2(CN ratio)\n(chromosomal arm level)",
  at = c(-0.4, -0.1, 0, 0.1, 0.4),
  direction = "horizontal",
  title_position = "topcenter",
  legend_width = unit(4, "cm"),
  grid_height = unit(4, "mm")
)
draw(legend_log2, x = unit(0.15, "npc"), y = unit(0.05, "npc"), just = c("center", "bottom"))