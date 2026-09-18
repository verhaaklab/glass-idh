# =============================================================================
# Create IDH-mutant GLASS gold-set purity and coverage table
#
# This script reads only these three source files:
#   - gold (From Glass5 analysis.gold_set)
#   - seqz_ascat_gold_mut (Table combining sequenza purity from glass5 analysis.seqz_params with ascat purity)
#   - glass5_coverage_all (Table combining WGS coverage values from glass5 database and WXS coverage values calculated with samtools depth)
#
# Supported input formats:
#   .csv, .tsv, .txt, and .rds
#
# Output:
#   ~/revision/idh_mut_gold_purity_coverage.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------

work_dir <- "~/revision"

# The script will try the exact path first, followed by:
# .csv, .tsv, .txt, and .rds
gold_file <- file.path(work_dir, "gold")
purity_file <- file.path(work_dir, "seqz_ascat_gold_mut")
coverage_file <- file.path(work_dir, "glass5_coverage_all")

output_file <- file.path(
  work_dir,
  "idh_mut_gold_purity_coverage.csv"
)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

find_input_file <- function(path_without_required_extension) {
  candidates <- unique(c(
    path_without_required_extension,
    paste0(path_without_required_extension, ".csv"),
    paste0(path_without_required_extension, ".tsv"),
    paste0(path_without_required_extension, ".txt"),
    paste0(path_without_required_extension, ".rds")
  ))

  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0) {
    stop(
      "Input file not found. Tried:\n",
      paste(candidates, collapse = "\n"),
      call. = FALSE
    )
  }

  if (length(existing) > 1) {
    warning(
      "Multiple matching files were found. Using:\n",
      existing[[1]]
    )
  }

  existing[[1]]
}

read_input_file <- function(path) {
  extension <- tolower(tools::file_ext(path))

  if (extension == "rds") {
    return(readRDS(path))
  }

  if (extension == "csv") {
    return(read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  if (extension %in% c("tsv", "txt")) {
    return(read.delim(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  # For extension-free files, try comma-separated first and then tab-separated.
  csv_attempt <- try(
    read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    silent = TRUE
  )

  if (
    !inherits(csv_attempt, "try-error") &&
    ncol(csv_attempt) > 1
  ) {
    return(csv_attempt)
  }

  read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

require_columns <- function(data, required, object_name) {
  missing_columns <- setdiff(required, names(data))

  if (length(missing_columns) > 0) {
    stop(
      object_name,
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))

  if (all(is.na(x))) {
    return(NA_real_)
  }

  max(x, na.rm = TRUE)
}

# -----------------------------------------------------------------------------
# Read source files
# -----------------------------------------------------------------------------

dir.create(
  work_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

gold_path <- find_input_file(gold_file)
purity_path <- find_input_file(purity_file)
coverage_path <- find_input_file(coverage_file)

message("Reading gold file: ", gold_path)
gold <- read_input_file(gold_path)

message("Reading purity file: ", purity_path)
seqz_ascat_gold_mut <- read_input_file(purity_path)

message("Reading coverage file: ", coverage_path)
glass5_coverage_all <- read_input_file(coverage_path)

# -----------------------------------------------------------------------------
# Validate source tables
# -----------------------------------------------------------------------------

if (!is.data.frame(gold)) {
  stop("The gold input must contain a data frame.", call. = FALSE)
}

if (ncol(gold) < 3) {
  stop(
    "The gold table must contain at least three columns. ",
    "Column 3 is expected to contain the initial aliquot barcode.",
    call. = FALSE
  )
}

require_columns(
  seqz_ascat_gold_mut,
  c(
    "case_barcode",
    "aliquot_barcode",
    "purity",
    "ascat_purity"
  ),
  "seqz_ascat_gold_mut"
)

require_columns(
  glass5_coverage_all,
  c(
    "case_barcode",
    "sample_name",
    "mean_coverage",
    "idh_codel_subtype",
    "gold"
  ),
  "glass5_coverage_all"
)

# -----------------------------------------------------------------------------
# Identify initial gold-set aliquots
#
# This follows the original history, where column 3 of gold represented the
# initial aliquot barcode.
# -----------------------------------------------------------------------------

initial_barcodes <- unique(
  trimws(as.character(gold[[3]]))
)

initial_barcodes <- initial_barcodes[
  !is.na(initial_barcodes) &
    initial_barcodes != ""
]

message(
  "Initial aliquot barcodes identified: ",
  length(initial_barcodes)
)

# -----------------------------------------------------------------------------
# Create case-level purity table
# -----------------------------------------------------------------------------

idh_mut_gold_purity <- seqz_ascat_gold_mut %>%
  mutate(
    sample_type = if_else(
      trimws(as.character(aliquot_barcode)) %in% initial_barcodes,
      "initial",
      "recurrent"
    ),
    purity = suppressWarnings(
      as.numeric(as.character(purity))
    ),
    ascat_purity = suppressWarnings(
      as.numeric(as.character(ascat_purity))
    )
  ) %>%
  select(
    case_barcode,
    sample_type,
    purity,
    ascat_purity
  ) %>%
  pivot_wider(
    names_from = sample_type,
    values_from = c(
      purity,
      ascat_purity
    ),
    names_sep = "_",
    values_fn = safe_max
  )

expected_purity_columns <- c(
  "purity_initial",
  "purity_recurrent",
  "ascat_purity_initial",
  "ascat_purity_recurrent"
)

for (column_name in expected_purity_columns) {
  if (!column_name %in% names(idh_mut_gold_purity)) {
    idh_mut_gold_purity[[column_name]] <- NA_real_
  }
}

idh_mut_gold_purity <- idh_mut_gold_purity %>%
  transmute(
    case_barcode,
    seqz_initial = purity_initial,
    seqz_recurrent = purity_recurrent,
    ascat_initial = ascat_purity_initial,
    ascat_recurrent = ascat_purity_recurrent
  )

# -----------------------------------------------------------------------------
# Create case-level coverage table
# -----------------------------------------------------------------------------

idh_mut_gold_coverage <- glass5_coverage_all %>%
  filter(
    !is.na(idh_codel_subtype),
    trimws(as.character(idh_codel_subtype)) != "IDHwt",
    tolower(trimws(as.character(gold))) == "yes"
  ) %>%
  mutate(
    sample_type = if_else(
      trimws(as.character(sample_name)) %in% initial_barcodes,
      "initial",
      "recurrent"
    ),
    mean_coverage = suppressWarnings(
      as.numeric(as.character(mean_coverage))
    )
  ) %>%
  select(
    case_barcode,
    sample_type,
    mean_coverage
  ) %>%
  pivot_wider(
    names_from = sample_type,
    values_from = mean_coverage,
    names_prefix = "cov_mean_",
    values_fn = safe_max
  )

expected_coverage_columns <- c(
  "cov_mean_initial",
  "cov_mean_recurrent"
)

for (column_name in expected_coverage_columns) {
  if (!column_name %in% names(idh_mut_gold_coverage)) {
    idh_mut_gold_coverage[[column_name]] <- NA_real_
  }
}

idh_mut_gold_coverage <- idh_mut_gold_coverage %>%
  select(
    case_barcode,
    cov_mean_initial,
    cov_mean_recurrent
  )

# -----------------------------------------------------------------------------
# Merge purity and coverage
#
# inner_join reproduces the behavior of merge() from the original history:
# only cases present in both tables are retained.
# -----------------------------------------------------------------------------

idh_mut_gold_purity_coverage <- inner_join(
  idh_mut_gold_purity,
  idh_mut_gold_coverage,
  by = "case_barcode"
) %>%
  arrange(case_barcode)

# -----------------------------------------------------------------------------
# Write final CSV
# -----------------------------------------------------------------------------

write.csv(
  idh_mut_gold_purity_coverage,
  file = output_file,
  row.names = FALSE,
  na = ""
)

message(
  "Rows written: ",
  nrow(idh_mut_gold_purity_coverage)
)

message(
  "Output file: ",
  normalizePath(
    output_file,
    winslash = "/",
    mustWork = FALSE
  )
)
