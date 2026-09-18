# =============================================================================
# Build IDH-mutant GLASS WGS TelSeq dataset
#
# Output:
#   idh_mut_gold_wgs_tel_length_norm.csv
#
# Workflow:
#   1. Query the GLASS gold set, existing TelSeq results, and molecular subtypes.
#   2. Retain IDH-mutant, 1p/19q-codeleted or non-codeleted WGS aliquots.
#   3. Recalculate TelSeq values from the raw missing-sample TelSeq output.
#   4. Fill only missing TelSeq fields.
#   5. Export one final CSV file.
# =============================================================================

suppressPackageStartupMessages({
  library(DBI)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------

dsn_name <- "glass5"

work_dir <- "~/glass.oligo/draft/revision"

raw_missing_telseq_file <- file.path(
  work_dir,
  "telseq_glass5_missing-gold_20251023.tsv"
)

output_file <- file.path(
  work_dir,
  "idh_mut_gold_wgs_tel_length_norm.csv"
)

# TelSeq constants used in the original analysis
TELSEQ_K <- 7L
GENOME_SIZE <- 332720800
READ_LENGTH_CONSTANT <- 46000

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

extract_barcode_field <- function(x, field_number) {
  vapply(
    strsplit(as.character(x), "-", fixed = TRUE),
    function(parts) {
      if (length(parts) >= field_number) {
        parts[[field_number]]
      } else {
        NA_character_
      }
    },
    character(1)
  )
}

safe_weighted_mean <- function(x, w) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  w <- suppressWarnings(as.numeric(as.character(w)))

  keep <- is.finite(x) & is.finite(w) & w > 0

  if (!any(keep)) {
    return(NA_real_)
  }

  stats::weighted.mean(x[keep], w[keep])
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

# -----------------------------------------------------------------------------
# Validate paths
# -----------------------------------------------------------------------------

dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(raw_missing_telseq_file)) {
  stop(
    "Raw missing-sample TelSeq file was not found:\n",
    raw_missing_telseq_file,
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# Query GLASS database
# -----------------------------------------------------------------------------

message("Connecting to ODBC DSN: ", dsn_name)

con <- DBI::dbConnect(
  odbc::odbc(),
  dsn = dsn_name
)

on.exit(
  try(DBI::dbDisconnect(con), silent = TRUE),
  add = TRUE
)

gold <- DBI::dbGetQuery(
  con,
  "SELECT * FROM analysis.gold_set"
)

telseq_existing <- DBI::dbGetQuery(
  con,
  "SELECT * FROM analysis.telseq"
)

subtypes <- DBI::dbGetQuery(
  con,
  "SELECT * FROM clinical.subtypes"
)

# -----------------------------------------------------------------------------
# Convert the paired gold-set aliquot columns to long format
# -----------------------------------------------------------------------------

if (ncol(gold) < 4) {
  stop(
    "analysis.gold_set must contain at least four columns. ",
    "The original workflow expects columns 3 and 4 to be the two aliquots.",
    call. = FALSE
  )
}

gold_a <- gold[, c(1, 2, 3), drop = FALSE]
gold_b <- gold[, c(1, 2, 4), drop = FALSE]

names(gold_a)[3] <- "aliquot_barcode"
names(gold_b)[3] <- "aliquot_barcode"
names(gold_b)[1:2] <- names(gold_a)[1:2]

gold_long <- bind_rows(gold_a, gold_b) %>%
  filter(
    !is.na(aliquot_barcode),
    aliquot_barcode != ""
  ) %>%
  distinct()

require_columns(
  gold_long,
  c("case_barcode", "aliquot_barcode"),
  "Long-format gold set"
)

require_columns(
  telseq_existing,
  "aliquot_barcode",
  "analysis.telseq"
)

require_columns(
  subtypes,
  c("case_barcode", "idh_codel_subtype"),
  "clinical.subtypes"
)

# Prevent accidental row expansion during joins
if (anyDuplicated(telseq_existing$aliquot_barcode)) {
  stop(
    "analysis.telseq contains duplicate aliquot_barcode values. ",
    "Resolve the duplicates before running this script.",
    call. = FALSE
  )
}

subtypes_one_row <- subtypes %>%
  distinct(case_barcode, .keep_all = TRUE)

# -----------------------------------------------------------------------------
# Build IDH-mutant WGS dataset
# -----------------------------------------------------------------------------

telseq_idh_wgs <- gold_long %>%
  left_join(
    telseq_existing,
    by = "aliquot_barcode"
  ) %>%
  left_join(
    subtypes_one_row,
    by = "case_barcode"
  ) %>%
  mutate(
    sequencing_type = extract_barcode_field(aliquot_barcode, 6)
  ) %>%
  filter(
    idh_codel_subtype %in% c(
      "IDHmut-codel",
      "IDHmut-noncodel"
    ),
    sequencing_type == "WGS"
  )

message(
  "IDH-mutant WGS aliquots identified: ",
  nrow(telseq_idh_wgs)
)

# -----------------------------------------------------------------------------
# Calculate TelSeq values from the raw missing-sample output
# -----------------------------------------------------------------------------

raw_telseq <- read.delim(
  raw_missing_telseq_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

require_columns(
  raw_telseq,
  c(
    "Sample",
    "Total",
    "Mapped",
    "Duplicates",
    "GC4",
    "GC5"
  ),
  "Raw missing-sample TelSeq file"
)

tel_columns <- grep(
  "^TEL[0-9]+$",
  names(raw_telseq),
  value = TRUE
)

if (length(tel_columns) == 0) {
  stop(
    "No TEL columns were found. Expected columns such as TEL0, TEL1, ...",
    call. = FALSE
  )
}

tel_indices <- as.integer(sub("^TEL", "", tel_columns))

raw_telseq$K <- TELSEQ_K
raw_telseq$G <- GENOME_SIZE
raw_telseq$c <- READ_LENGTH_CONSTANT

numeric_columns <- unique(
  c(
    "Total",
    "Mapped",
    "Duplicates",
    "GC4",
    "GC5",
    tel_columns,
    "K",
    "G",
    "c"
  )
)

raw_telseq[numeric_columns] <- lapply(
  raw_telseq[numeric_columns],
  function(x) suppressWarnings(as.numeric(as.character(x)))
)

tel_matrix <- as.matrix(
  raw_telseq[, tel_columns, drop = FALSE]
)

tel_count <- vapply(
  seq_len(nrow(raw_telseq)),
  function(i) {
    qualifying_columns <- which(
      tel_indices >= raw_telseq$K[[i]]
    )

    if (length(qualifying_columns) == 0) {
      return(NA_real_)
    }

    sum(
      tel_matrix[i, qualifying_columns, drop = TRUE],
      na.rm = TRUE
    )
  },
  numeric(1)
)

gc_count <- raw_telseq$GC4 + raw_telseq$GC5

tel_length_norm <- ifelse(
  is.finite(gc_count) & gc_count > 0,
  (tel_count / gc_count) *
    (raw_telseq$G / raw_telseq$c),
  NA_real_
)

raw_telseq_enhanced <- raw_telseq %>%
  mutate(
    tel = tel_count,
    gc = gc_count,
    length = tel_length_norm
  )

calculated_telseq <- raw_telseq_enhanced %>%
  group_by(Sample) %>%
  summarise(
    total_reads = round(sum(Total, na.rm = TRUE)),
    mapped_reads = round(sum(Mapped, na.rm = TRUE)),
    duplicate_reads = round(sum(Duplicates, na.rm = TRUE)),
    tel = safe_weighted_mean(tel, Total),
    gc = safe_weighted_mean(gc, Total),
    length = safe_weighted_mean(length, Total),
    K = first(K),
    G = first(G),
    c = first(c),
    .groups = "drop"
  ) %>%
  rename(
    aliquot_barcode = Sample
  ) %>%
  select(
    aliquot_barcode,
    total_reads,
    mapped_reads,
    duplicate_reads,
    tel,
    K,
    G,
    c,
    gc,
    length
  )

if (anyDuplicated(calculated_telseq$aliquot_barcode)) {
  stop(
    "Calculated TelSeq table contains duplicate aliquot barcodes.",
    call. = FALSE
  )
}

message(
  "Recalculated TelSeq aliquots available: ",
  nrow(calculated_telseq)
)

# -----------------------------------------------------------------------------
# Fill missing TelSeq values without overwriting existing database values
# -----------------------------------------------------------------------------

columns_to_fill <- c(
  "total_reads",
  "mapped_reads",
  "duplicate_reads",
  "tel",
  "K",
  "G",
  "c",
  "gc",
  "length"
)

telseq_idh_wgs_filled <- telseq_idh_wgs %>%
  left_join(
    calculated_telseq,
    by = "aliquot_barcode",
    suffix = c("", ".calculated")
  )

for (column_name in columns_to_fill) {
  calculated_column <- paste0(
    column_name,
    ".calculated"
  )

  if (!calculated_column %in% names(telseq_idh_wgs_filled)) {
    next
  }

  calculated_values <- suppressWarnings(
    as.numeric(
      as.character(
        telseq_idh_wgs_filled[[calculated_column]]
      )
    )
  )

  if (column_name %in% names(telseq_idh_wgs_filled)) {
    existing_values <- suppressWarnings(
      as.numeric(
        as.character(
          telseq_idh_wgs_filled[[column_name]]
        )
      )
    )

    telseq_idh_wgs_filled[[column_name]] <- dplyr::coalesce(
      existing_values,
      calculated_values
    )
  } else {
    telseq_idh_wgs_filled[[column_name]] <- calculated_values
  }
}

telseq_idh_wgs_filled <- telseq_idh_wgs_filled %>%
  select(
    -ends_with(".calculated")
  ) %>%
  mutate(
    tel_length_norm = length
  ) %>%
  arrange(
    case_barcode,
    aliquot_barcode
  )

# -----------------------------------------------------------------------------
# Export final CSV
# -----------------------------------------------------------------------------

write.csv(
  telseq_idh_wgs_filled,
  file = output_file,
  row.names = FALSE,
  na = ""
)

message("Final rows written: ", nrow(telseq_idh_wgs_filled))
message("Output file: ", normalizePath(output_file, winslash = "/", mustWork = FALSE))
