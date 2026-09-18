# ============================================================
# FGA calculation and annotation
#
# Inputs already loaded / created:
#   segment data from sequenza calls
#   gold set with IDH-mut subtypes included
#
# Final output:
#   fga_calc_annotated.txt
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ------------------------------------------------------------
# 1. Parameters
# ------------------------------------------------------------

log2_threshold <- 0.2

# Approximate callable genome sizes
callable_bp_wxs <- 50e6
callable_bp_wgs <- 3e9


# ------------------------------------------------------------
# 2. Combine segmentation data
# ------------------------------------------------------------

seg_data <- bind_rows(
  codel_p,
  codel_r,
  noncodel_p,
  noncodel_r
)


# ------------------------------------------------------------
# 3. Determine sequencing platform
#
# Example aliquot barcode:
# fields separated by "-"
# sixth field expected to be WGS or WXS
# ------------------------------------------------------------

seg_data <- seg_data %>%
  mutate(
    platform = sapply(
      strsplit(as.character(aliquot_barcode), "-"),
      function(x) {
        if (length(x) >= 6) x[6] else NA_character_
      }
    )
  )


# ------------------------------------------------------------
# 4. Calculate FGA for each sample
#
# Altered segment:
#   abs(log2_copy_ratio) > 0.2
#
# FGA:
#   total altered bp / callable bp
# ------------------------------------------------------------

cnv_summary <- seg_data %>%
  mutate(
    seg_size = end - start,
    altered = abs(log2_copy_ratio) > log2_threshold
  ) %>%
  group_by(
    case_barcode,
    aliquot_barcode,
    platform
  ) %>%
  summarise(
    n_segments = sum(altered, na.rm = TRUE),

    altered_bp = sum(
      seg_size[altered],
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  mutate(
    total_bp = case_when(
      platform == "WGS" ~ callable_bp_wgs,
      platform == "WXS" ~ callable_bp_wxs,
      TRUE ~ NA_real_
    ),

    fga = altered_bp / total_bp
  )


# ------------------------------------------------------------
# 5. Annotate using gold cohort information
#
# gold_long_subtypes:
#   case_barcode
#   tumor_barcode
#   seq
#   idh_codel_subtype
# ------------------------------------------------------------

cnv_summary_gold <- gold_long_subtypes %>%
  inner_join(
    cnv_summary,
    by = c(
      "case_barcode" = "case_barcode",
      "tumor_barcode" = "aliquot_barcode",
      "seq" = "platform"
    )
  )


# ------------------------------------------------------------
# 6. Normalize FGA separately for WGS and WXS
#
# This reproduces:
#
# group_by(seq) %>%
# mutate(
#   fga_norm = (fga - min(fga)) /
#              (max(fga) - min(fga))
# )
# ------------------------------------------------------------

cnv_summary_gold <- cnv_summary_gold %>%
  group_by(seq) %>%
  mutate(
    fga_norm = if_else(
      max(fga, na.rm = TRUE) > min(fga, na.rm = TRUE),

      (fga - min(fga, na.rm = TRUE)) /
        (max(fga, na.rm = TRUE) - min(fga, na.rm = TRUE)),

      NA_real_
    )
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 7. Create primary-sample FGA lookup
#
# Your history referenced:
#   cnv_summary_filter_primary
#
# but did not show how it was created.
#
# Here we define primary samples using sample_type if that
# column exists.
# ------------------------------------------------------------

if ("sample_type" %in% colnames(cnv_summary_gold)) {

  cnv_summary_filter_primary <- cnv_summary_gold %>%
    filter(
      tolower(sample_type) %in%
        c("primary", "initial")
    )

} else {

  # If gold_long_subtypes contains one row per tumor already,
  # use the existing table as the lookup.
  cnv_summary_filter_primary <- cnv_summary_gold
}


# ------------------------------------------------------------
# 8. Primary FGA lookup table
# ------------------------------------------------------------

primary_fga <- cnv_summary_filter_primary %>%
  select(
    case_barcode,
    aliquot_barcode = tumor_barcode,
    fga_norm
  ) %>%
  distinct() %>%
  rename(
    fga_norm_primary = fga_norm
  )


# ------------------------------------------------------------
# 9. Add primary-sample normalized FGA
# ------------------------------------------------------------

cnv_summary_gold <- cnv_summary_gold %>%
  left_join(
    primary_fga,
    by = c(
      "case_barcode",
      "tumor_barcode" = "aliquot_barcode"
    )
  )


# ------------------------------------------------------------
# 10. Restrict to IDH-mutant tumors
#
# Original:
#   filter(!idh_codel_subtype == "IDHwt")
# ------------------------------------------------------------

cnv_summary_gold_mut <- cnv_summary_gold %>%
  filter(
    !is.na(idh_codel_subtype),
    idh_codel_subtype != "IDHwt"
  )


# ------------------------------------------------------------
# 11. Re-normalize FGA within IDH-mutant samples only
#
# This reproduces fga_norm.z from the R history.
# ------------------------------------------------------------

cnv_summary_gold_mut <- cnv_summary_gold_mut %>%
  group_by(seq) %>%
  mutate(
    fga_norm_idhmut = if_else(
      max(fga, na.rm = TRUE) > min(fga, na.rm = TRUE),

      (fga - min(fga, na.rm = TRUE)) /
        (max(fga, na.rm = TRUE) - min(fga, na.rm = TRUE)),

      NA_real_
    )
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 12. Calculate delta
#
# Equivalent to:
#   fga_norm.y - fga_norm.z
# ------------------------------------------------------------

cnv_summary_gold_mut <- cnv_summary_gold_mut %>%
  mutate(
    delta = fga_norm_primary - fga_norm_idhmut
  )


# ------------------------------------------------------------
# 13. Final FGA
#
# Original:
#   FGA = coalesce(fga_norm.y, fga_norm.z)
#
# Prefer primary normalized FGA when available.
# Otherwise use IDH-mutant normalized FGA.
# ------------------------------------------------------------

cnv_summary_gold_mut <- cnv_summary_gold_mut %>%
  mutate(
    FGA = coalesce(
      fga_norm_primary,
      fga_norm_idhmut
    )
  )


# ------------------------------------------------------------
# 14. Create final annotated FGA table
# ------------------------------------------------------------

fga_calc_annotated <- cnv_summary_gold_mut %>%
  select(
    case_barcode,
    aliquot_barcode = tumor_barcode,
    seq,
    idh_codel_subtype,
    n_segments,
    altered_bp,
    total_bp,
    fga,
    fga_norm,
    fga_norm_primary,
    fga_norm_idhmut,
    delta,
    FGA
  ) %>%
  distinct()


# ------------------------------------------------------------
# 15. Write final output
# ------------------------------------------------------------

write.table(
  fga_calc_annotated,
  file = "fga_calc_annotated.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


# ------------------------------------------------------------
# 16. Also reproduce the smaller table from the original
# history, if desired
# ------------------------------------------------------------

fga_seqz <- fga_calc_annotated %>%
  select(
    case_barcode,
    aliquot_barcode,
    idh_codel_subtype,
    FGA
  )

write.csv(
  fga_seqz,
  file = "fga_seqz_gold_IDHmut.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# Done
# ------------------------------------------------------------

message("FGA calculation complete.")
message("Output: fga_calc_annotated.txt")
message("Output: fga_seqz_gold_IDHmut.csv")