# ============================================================
# WGD calling from Sequenza segment table
# Date: 2026.01.21
# Author: Tamrin Chowdhury
#
#
# INPUT columns expected:
#  pair_barcode, chrom, pos, major_cn, minor_cn, cellularity, ploidy
# (other columns can exist; they'll be ignored)
# Generate WGD_sequenza_all.csv directly from glass5 database
# ============================================================

suppressPackageStartupMessages({
  library(DBI)
  library(dplyr)
  library(stringr)
})

# ------------------------- SETTINGS -------------------------
DSN <- "glass5"
OUT_FILE <- "WGD_sequenza_all.csv"

WGD_FRAC_CUTOFF <- 0.50
PURITY_CUTOFF <- 0.40
PLOIDY_CUTOFF <- 3
# ------------------------------------------------------------


is_autosome <- function(chr) {
  chr <- str_remove(as.character(chr), regex("^chr", ignore_case = TRUE))
  x <- suppressWarnings(as.integer(chr))
  !is.na(x) & x >= 1 & x <= 22
}


parse_pos <- function(pos) {
  x <- str_remove_all(as.character(pos), "[\\[\\]\\(\\)\\s]")
  p <- str_split_fixed(x, ",", 2)

  tibble(
    start = suppressWarnings(as.numeric(p[, 1])),
    end   = suppressWarnings(as.numeric(p[, 2]))
  )
}


barcode_field <- function(x, n) {
  p <- str_split_fixed(as.character(x), "-", n)
  out <- p[, n]
  out[out == ""] <- NA_character_
  out
}


main <- function() {

  # ---- Connect once and retrieve only the tables needed ----
  con <- dbConnect(odbc::odbc(), dsn = DSN)
  on.exit(dbDisconnect(con), add = TRUE)

  seqz_seg <- dbGetQuery(
    con,
    "SELECT * FROM variants.seqz_seg"
  )

  seqz_params <- dbGetQuery(
    con,
    "SELECT * FROM variants.seqz_params"
  )

  subtypes <- dbGetQuery(
    con,
    "SELECT case_barcode, idh_codel_subtype FROM clinical.subtypes"
  )

  aliquots <- dbGetQuery(
    con,
    paste(
      "SELECT aliquot_barcode, sample_barcode, aliquot_analysis_type",
      "FROM biospecimen.aliquots"
    )
  )


  # ---- Sequenza: calculate length-weighted frac_major_ge2 ----
  segs <- merge(
    seqz_seg,
    seqz_params,
    by = "pair_barcode",
    all.x = TRUE
  )

  required <- c(
    "pair_barcode", "chrom", "pos",
    "major_cn", "minor_cn",
    "cellularity", "ploidy"
  )

  missing <- setdiff(required, names(segs))
  if (length(missing) > 0) {
    stop("Missing required column(s): ", paste(missing, collapse = ", "))
  }

  dat <- segs %>%
    transmute(
      pair_barcode = as.character(pair_barcode),
      chrom = str_remove(
        as.character(chrom),
        regex("^chr", ignore_case = TRUE)
      ),
      pos = as.character(pos),
      major_cn = as.numeric(major_cn),
      minor_cn = as.numeric(minor_cn),
      purity = as.numeric(cellularity),
      ploidy = as.numeric(ploidy)
    ) %>%
    bind_cols(parse_pos(.$pos)) %>%
    mutate(seg_len = end - start) %>%
    filter(
      is_autosome(chrom),
      is.finite(start),
      is.finite(end),
      is.finite(seg_len),
      seg_len > 0,
      is.finite(major_cn),
      is.finite(minor_cn)
    )

  if (nrow(dat) == 0) {
    stop("No valid autosomal Sequenza segments remained after filtering.")
  }

  summ <- dat %>%
    group_by(pair_barcode) %>%
    summarise(
      total_len = sum(seg_len, na.rm = TRUE),
      frac_major_ge2 =
        sum(if_else(major_cn >= 2, seg_len, 0), na.rm = TRUE) /
        total_len,
      purity = median(purity, na.rm = TRUE),
      ploidy = median(ploidy, na.rm = TRUE),
      .groups = "drop"
    )

  f1 <- barcode_field(summ$pair_barcode, 1)
  f2 <- barcode_field(summ$pair_barcode, 2)
  f3 <- barcode_field(summ$pair_barcode, 3)
  f4 <- barcode_field(summ$pair_barcode, 4)
  f5 <- barcode_field(summ$pair_barcode, 5)
  f8 <- barcode_field(summ$pair_barcode, 8)

  summ <- summ %>%
    mutate(
      case_barcode = str_c(f1, f2, f3, sep = "-"),
      sample_barcode = str_c(f1, f2, f3, f4, sep = "-"),
      seq = f8,
      aliquot = if_else(!is.na(f5), paste0(f5, "D"), NA_character_)
    ) %>%
    left_join(
      subtypes,
      by = "case_barcode"
    )


  # ---- Map pair_barcode-derived IDs to aliquot_barcode ----
  aliquot_map <- aliquots %>%
    transmute(
      aliquot_barcode = as.character(aliquot_barcode),
      sample_barcode = as.character(sample_barcode),
      seq = as.character(aliquot_analysis_type),
      aliquot = barcode_field(aliquot_barcode, 5)
    )

  summ_all <- summ %>%
    left_join(
      aliquot_map,
      by = c("sample_barcode", "seq", "aliquot")
    )


  # Final table and final WGD rule
  wgd_sequenza_all <- summ_all %>%
    select(
      case_barcode,
      aliquot_barcode,
      sample_barcode,
      pair_barcode,
      frac_major_ge2,
      purity,
      ploidy,
      idh_codel_subtype
    ) %>%
    mutate(
      WGD_status_pre = if_else(
        frac_major_ge2 >= WGD_FRAC_CUTOFF,
        "WGD",
        "No_WGD"
      ),
      confidence = if_else(
        purity >= PURITY_CUTOFF & ploidy > PLOIDY_CUTOFF,
        "High",
        "Low"
      ),
      WGD_status_final = case_when(
        WGD_status_pre == "WGD" & confidence == "High" ~ "WGD",
        WGD_status_pre == "No_WGD" ~ "No_WGD",
        TRUE ~ "WGD_uncertain"
      )
    )


  # ---- Write final file ----
  write.csv(
    wgd_sequenza_all,
    file = OUT_FILE,
    row.names = FALSE
  )

  message("Wrote: ", normalizePath(OUT_FILE, mustWork = FALSE))
  message("Rows: ", nrow(wgd_sequenza_all))
  print(table(wgd_sequenza_all$WGD_status_final, useNA = "ifany"))

  invisible(wgd_sequenza_all)
}


wgd_sequenza_all <- main()
