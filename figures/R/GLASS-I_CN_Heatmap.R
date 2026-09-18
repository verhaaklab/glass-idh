# GLASS-I: Copy Number Heatmap
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# ── Load required libraries ──
library(tidyverse)
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)
library(BSgenome.Hsapiens.UCSC.hg38)
library(colorspace)
library(grid)

# ── Load cytoband data ──
cyto <- read.cytoband(species = "hg19")$df
colnames(cyto) <- c("chrom", "start", "end", "band", "gieStain")
cyto$arm <- ifelse(grepl("^p", cyto$band), "p", "q")

cyto <- cyto %>%
  group_by(chr = chrom, arm) %>%
  summarise(start = min(start), end = max(end), .groups = "drop") %>%
  filter(chr %in% paste0("chr", 1:22))

cyto$chr_arm <- paste0(cyto$chr, cyto$arm)

gr_arms <- GRanges(seqnames = cyto$chr,
                   ranges = IRanges(start = cyto$start, end = cyto$end),
                   arm = cyto$chr_arm)

# ── Chromosome arm order ──
chroms <- paste0("chr", 1:22)
arms <- c("p", "q")
arm_order <- as.vector(t(outer(chroms, arms, paste0)))

# ── Function to process one CNV file ──
process_file <- function(file_path, group_label) {
  seg <- read.csv(file_path)
  seg$chrom <- ifelse(grepl("^chr", seg$chrom), seg$chrom, paste0("chr", seg$chrom))
  seg <- seg %>% filter(!is.na(log2_copy_ratio))
  
  gr_seg <- GRanges(seqnames = seg$chrom,
                    ranges = IRanges(start = seg$start, end = seg$end),
                    score = seg$log2_copy_ratio,
                    sample = seg$aliquot_barcode)
  
  hits <- findOverlaps(gr_seg, gr_arms)
  df <- data.frame(
    sample = mcols(gr_seg)$sample[queryHits(hits)],
    arm = mcols(gr_arms)$arm[subjectHits(hits)],
    score = mcols(gr_seg)$score[queryHits(hits)]
  )
  
  df %>%
    group_by(sample, arm) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    mutate(group = group_label)
}

# ── Load all 4 input files ──
df1 <- process_file("~/gistic_circos_input_codel_primary.csv", "Oligo. Primary")                            ## created with prepare_gistic_input.R
df2 <- process_file("~/gistic_circos_input_codel_recurrence.csv", "Oligo. Recurrence")                      ## created with prepare_gistic_input.R
df3 <- process_file("~/gistic_circos_input_noncodel_primary.csv", "Astro. Primary")                         ## created with prepare_gistic_input.R
df4 <- process_file("~/gistic_circos_input_noncodel_recurrence.csv", "Astro. Recurrence")                   ## created with prepare_gistic_input.R

# ── Combine data ──
all_df <- bind_rows(df1, df2, df3, df4)

# ── Row group labels for heatmap ──
group_df <- all_df %>% select(sample, group) %>% distinct()

# ── Create wide matrix ──
wide_df <- all_df %>%
  select(-group) %>%
  pivot_wider(names_from = arm, values_from = mean_score) %>%
  column_to_rownames("sample")

# ── Ensure all arms are present ──
for (arm in arm_order) {
  if (!(arm %in% colnames(wide_df))) {
    wide_df[[arm]] <- NA
  }
}
wide_df <- wide_df[, arm_order]

# ── Remove arms with all NA values ──
wide_df <- wide_df[, colSums(!is.na(wide_df)) > 0]

# ── Transpose ──
transposed_df <- t(as.matrix(wide_df))

# ── Update arm order and color ──
arm_order_filtered <- rownames(transposed_df)
base_colors <- qualitative_hcl(length(chroms), palette = "Dark 3")
names(base_colors) <- chroms

arm_colors <- unlist(lapply(chroms, function(chr) {
  c(lighten(base_colors[chr], 0.4), darken(base_colors[chr], 0.2))
}))
names(arm_colors) <- as.vector(t(outer(chroms, arms, paste0)))
arm_colors_filtered <- arm_colors[arm_order_filtered]

# ── Row annotation (for chrom arms, now rows) ──
right_anno <- rowAnnotation(
  " " = anno_simple(x = arm_order_filtered, col = arm_colors_filtered)
)

# ── Heatmap color function ──
col_fun <- colorRamp2(
  c(-0.6, -0.2, 0, 0.2, 0.6),
  c("#2066a8", "#8ec1da", "white", "#f6d6c2", "#a00000")
)

# ── Match group labels ──
row_groups <- group_df$group[match(colnames(transposed_df), group_df$sample)]

# ── Draw heatmap ──
Heatmap(transposed_df,
        name = "log2 CN ratio",
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 6),
        column_split = factor(row_groups, levels = c(
          "Oligo. Primary", "Oligo. Recurrence",
          "Astro. Primary", "Astro. Recurrence"
        )),
        column_gap = unit(4, "mm"),
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        border = TRUE,
        right_annotation = right_anno,
        heatmap_legend_param = list(
          title = "log2(CN ratio)",
          at = c(-0.4, -0.1, 0, 0.1, 0.4)
        ))