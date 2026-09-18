#!/usr/bin/env bash
#SBATCH --job-name=mut_sig_pipeline
#SBATCH --partition=day
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# Coherent mutational signature pipeline for:
#   1) RepeatMasker/ENCODE repeat-region filtering from GLASS5
#   2) SigProfiler input generation and WGS/WXS split
#   3) SigProfilerMatrixGenerator + SigProfilerExtractor SBS/DBS/ID extraction
#   4) Palimpsest reference deconvolution and per-mutation signature origins
#
# This script writes small helper R/Python scripts into WORKDIR/_pipeline_scripts
# so the full workflow can be launched from one command while keeping each stage
# resume-aware via marker files.

set -Eeuo pipefail
IFS=$'\n\t'

SCRIPT_NAME=$(basename "$0")

usage() {
  cat <<'EOF'
Usage:
  mutational_signature_pipeline.sh --workdir WORKDIR --sql final_query.sql [options]

Required for full filter/input run:
  --workdir PATH                  Main pipeline working directory.
  --sql PATH                      SQL query used for repeat-region filtering.

Required only if running palimpsest:
  --reference-rdata PATH          RData containing Palimpsest reference signatures.

Common options:
  --steps LIST                    Comma list: filter,input,matrix,extract,palimpsest,all [default: all]
  --data-types LIST               Comma list: wgs,wxs [default: wgs]
  --sig-types LIST                Comma list: SBS,DBS,ID [default: SBS,DBS,ID]
  --dsn NAME                      ODBC DSN for GLASS5 database [default: glass5]
  --subtype VALUE                 clinical.subtypes idh_codel_subtype filter [default: IDHmut-codel]
  --project-name VALUE            Project label written to SigProfiler input [default: GLASS5]
  --genome VALUE                  Genome label for SigProfiler input [default: GRCh37]
  --matrix-project VALUE          MatrixGenerator project/cohort prefix [default: codel]
  --min-signatures INT            SigProfilerExtractor minimum signatures [default: 1]
  --max-signatures INT            SigProfilerExtractor maximum signatures [default: 20]
  --cushion INT                   MatrixGenerator cushion [default: 100]
  --r-env NAME                    Conda environment with R/Palimpsest [default: palimpsest]
  --sigprofiler-env NAME          Conda environment with SigProfiler [default: sigprofiler]
  --miniconda-module NAME         Module to load before conda activation [default: miniconda]
  --no-module-load                Do not run module load before conda activate.
  --no-conda                      Do not activate conda environments.
  --force                         Rerun completed stages and overwrite markers.
  --dry-run                       Print planned actions without running them.

Palimpsest-specific options:
  --sbs-ref-object NAME           SBS reference object in --reference-rdata [default: SBS_cosmic3.2_sigs]
  --dbs-ref-object NAME           DBS reference object in --reference-rdata [default: dbs_ref_codel]
  --id-ref-object NAME            ID reference object in --reference-rdata [default: ID_cosmic]
  --id83-reference-file PATH      Optional SigProfiler COSMIC_ID83_Signatures.txt; overrides --id-ref-object for ID.
  --ref-fasta PATH                Optional hg19/GRCh37 FASTA for Palimpsest ID annotation.
  --sample-info-object NAME       Optional sample metadata object in RData for shared/private fractions [default: codel_info]

Useful examples:

  # Run filtering + input creation only on SLURM
  sbatch mutational_signature_pipeline.sh \
    --workdir /path/to/mut_sig_run \
    --sql /path/to/final_query.sql \
    --steps filter,input \
    --r-env palimpsest

  # Run WGS SigProfiler extraction only after input files are ready
  sbatch mutational_signature_pipeline.sh \
    --workdir /path/to/mut_sig_run \
    --steps matrix,extract \
    --data-types wgs \
    --sig-types SBS,DBS,ID \
    --sigprofiler-env sigprofiler

  # Run Palimpsest deconvolution only
  sbatch mutational_signature_pipeline.sh \
    --workdir /path/to/mut_sig_run \
    --steps palimpsest \
    --data-types wgs,wxs \
    --sig-types SBS,DBS,ID \
    --reference-rdata /path/to/signature_references.RData \
    --r-env palimpsest
EOF
}

WORKDIR=""
SQL=""
STEPS="all"
DATA_TYPES="wgs"
SIG_TYPES="SBS,DBS,ID"
DSN="glass5"
SUBTYPE="IDHmut-codel"
PROJECT_NAME="GLASS5"
GENOME="GRCh37"
MATRIX_PROJECT="codel"
MIN_SIGNATURES="1"
MAX_SIGNATURES="20"
CUSHION="100"
R_ENV="palimpsest"
SIGPROFILER_ENV="sigprofiler"
MINICONDA_MODULE="miniconda"
LOAD_MODULE="1"
USE_CONDA="1"
FORCE="0"
DRY_RUN="0"
REFERENCE_RDATA=""
SBS_REF_OBJECT="SBS_cosmic3.2_sigs"
DBS_REF_OBJECT="dbs_ref_codel"
ID_REF_OBJECT="ID_cosmic"
ID83_REFERENCE_FILE=""
REF_FASTA=""
SAMPLE_INFO_OBJECT="codel_info"
PLOT_MATRIX="true"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workdir) WORKDIR="$2"; shift 2 ;;
    --sql) SQL="$2"; shift 2 ;;
    --steps) STEPS="$2"; shift 2 ;;
    --data-types) DATA_TYPES="$2"; shift 2 ;;
    --sig-types) SIG_TYPES="$2"; shift 2 ;;
    --dsn) DSN="$2"; shift 2 ;;
    --subtype) SUBTYPE="$2"; shift 2 ;;
    --project-name) PROJECT_NAME="$2"; shift 2 ;;
    --genome) GENOME="$2"; shift 2 ;;
    --matrix-project) MATRIX_PROJECT="$2"; shift 2 ;;
    --min-signatures) MIN_SIGNATURES="$2"; shift 2 ;;
    --max-signatures) MAX_SIGNATURES="$2"; shift 2 ;;
    --cushion) CUSHION="$2"; shift 2 ;;
    --r-env) R_ENV="$2"; shift 2 ;;
    --sigprofiler-env) SIGPROFILER_ENV="$2"; shift 2 ;;
    --miniconda-module) MINICONDA_MODULE="$2"; shift 2 ;;
    --no-module-load) LOAD_MODULE="0"; shift ;;
    --no-conda) USE_CONDA="0"; shift ;;
    --force) FORCE="1"; shift ;;
    --dry-run) DRY_RUN="1"; shift ;;
    --reference-rdata) REFERENCE_RDATA="$2"; shift 2 ;;
    --sbs-ref-object) SBS_REF_OBJECT="$2"; shift 2 ;;
    --dbs-ref-object) DBS_REF_OBJECT="$2"; shift 2 ;;
    --id-ref-object) ID_REF_OBJECT="$2"; shift 2 ;;
    --id83-reference-file) ID83_REFERENCE_FILE="$2"; shift 2 ;;
    --ref-fasta) REF_FASTA="$2"; shift 2 ;;
    --sample-info-object) SAMPLE_INFO_OBJECT="$2"; shift 2 ;;
    --no-matrix-plots) PLOT_MATRIX="false"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "$WORKDIR" ]]; then
  echo "ERROR: --workdir is required." >&2
  usage
  exit 2
fi

WORKDIR=$(readlink -m "$WORKDIR")
SCRIPT_DIR="$WORKDIR/_pipeline_scripts"
LOG_DIR="$WORKDIR/logs"
FILTER_DIR="$WORKDIR/filter_step"
INPUT_DIR="$WORKDIR/sigprofiler_input"
SP_OUTDIR="$WORKDIR/sigprofiler_extraction"
PAL_OUTDIR="$WORKDIR/palimpsest_deconvolution"
MARKER_DIR="$WORKDIR/.pipeline_markers"
FILTER_RDATA="$FILTER_DIR/filtered_repeat_region_variants.RData"
FILTER_STATS="$FILTER_DIR/filter_repeat_region_stats.tsv"
INPUT_SUMMARY="$INPUT_DIR/sigprofiler_input_summary.tsv"

mkdir -p "$WORKDIR" "$SCRIPT_DIR" "$LOG_DIR" "$FILTER_DIR" "$INPUT_DIR" "$SP_OUTDIR" "$PAL_OUTDIR" "$MARKER_DIR"

if [[ -n "$SQL" ]]; then
  SQL=$(readlink -m "$SQL")
fi
if [[ -n "$REFERENCE_RDATA" ]]; then
  REFERENCE_RDATA=$(readlink -m "$REFERENCE_RDATA")
fi
if [[ -n "$ID83_REFERENCE_FILE" ]]; then
  ID83_REFERENCE_FILE=$(readlink -m "$ID83_REFERENCE_FILE")
fi
if [[ -n "$REF_FASTA" ]]; then
  REF_FASTA=$(readlink -m "$REF_FASTA")
fi

log() {
  local ts
  ts=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[$ts] $*"
}

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

run_cmd() {
  log "+ $*"
  if [[ "$DRY_RUN" == "1" ]]; then
    return 0
  fi
  "$@"
}

marker_path() {
  local name="$1"
  echo "$MARKER_DIR/${name}.done"
}

is_done() {
  local name="$1"
  [[ "$FORCE" != "1" && -s "$(marker_path "$name")" ]]
}

mark_done() {
  local name="$1"
  if [[ "$DRY_RUN" == "1" ]]; then
    return 0
  fi
  printf '%s\n' "completed_at=$(date -Iseconds)" > "$(marker_path "$name")"
}

contains_item() {
  local needle="$1"
  local haystack=",${2},"
  [[ "$haystack" == *",${needle},"* ]]
}

has_step() {
  local step="$1"
  [[ "$STEPS" == "all" || "$STEPS" == "ALL" ]] && return 0
  contains_item "$step" "$STEPS"
}

activate_env() {
  local env_name="$1"
  if [[ "$USE_CONDA" != "1" ]]; then
    return 0
  fi
  if [[ "$LOAD_MODULE" == "1" ]]; then
    # Some clusters use "module" as a shell function; ignore if absent.
    if command -v module >/dev/null 2>&1 || type module >/dev/null 2>&1; then
      module load "$MINICONDA_MODULE" || true
    fi
  fi
  if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
    conda activate "$env_name"
  else
    fail "conda not found. Use --no-conda if the environment is already active."
  fi
}

write_helpers() {
  local FILTER_R="$SCRIPT_DIR/01_filter_repeat_regions.R"
  local INPUT_R="$SCRIPT_DIR/02_make_sigprofiler_input.R"
  local SP_PY="$SCRIPT_DIR/03_run_sigprofiler.py"
  local PAL_R="$SCRIPT_DIR/04_run_palimpsest_deconvolution.R"

  cat > "$FILTER_R" <<'RSCRIPT'
suppressPackageStartupMessages({
  library(DBI)
  library(odbc)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) stop("Missing value for ", flag)
  args[[hit + 1]]
}

workdir <- normalizePath(get_arg("--workdir"), mustWork = FALSE)
sql_path <- normalizePath(get_arg("--sql"), mustWork = TRUE)
dsn <- get_arg("--dsn", "glass5")
out_rdata <- normalizePath(get_arg("--out-rdata"), mustWork = FALSE)
stats_path <- normalizePath(get_arg("--stats", file.path(workdir, "filter_repeat_region_stats.tsv")), mustWork = FALSE)

dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_rdata), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(stats_path), showWarnings = FALSE, recursive = TRUE)
setwd(workdir)

message("Connecting to ODBC DSN: ", dsn)
con <- DBI::dbConnect(odbc::odbc(), dsn = dsn)
on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

sql <- paste(readLines(sql_path, warn = FALSE), collapse = "\n")
message("Running repeat-region SQL query: ", sql_path)
df <- DBI::dbGetQuery(con, sql)
names(df) <- tolower(names(df))

required <- c("rm_type", "en_type")
missing <- setdiff(required, names(df))
if (length(missing) > 0) stop("SQL output missing required columns: ", paste(missing, collapse = ", "))

df <- df %>% mutate(filter_call = ifelse(is.na(.data$rm_type) & is.na(.data$en_type), "pass", "fail"))
df_filtered <- df %>% filter(.data$filter_call == "pass")

stats <- df %>% count(.data$filter_call, name = "n_variants")
write.table(stats, file = stats_path, sep = "\t", row.names = FALSE, quote = FALSE)
save(df, df_filtered, file = out_rdata)

message("Saved filtered RData: ", out_rdata)
message("Saved stats TSV: ", stats_path)
RSCRIPT

  cat > "$INPUT_R" <<'RSCRIPT'
suppressPackageStartupMessages({
  library(DBI)
  library(odbc)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) stop("Missing value for ", flag)
  args[[hit + 1]]
}

rdata <- normalizePath(get_arg("--filter-rdata"), mustWork = TRUE)
outdir <- normalizePath(get_arg("--outdir"), mustWork = FALSE)
dsn <- get_arg("--dsn", "glass5")
subtype_value <- get_arg("--subtype", "IDHmut-codel")
project_name <- get_arg("--project-name", "GLASS5")
genome <- get_arg("--genome", "GRCh37")
summary_path <- normalizePath(get_arg("--summary", file.path(outdir, "sigprofiler_input_summary.tsv")), mustWork = FALSE)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "all"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "wgs"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "wxs"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "unknown"), showWarnings = FALSE, recursive = TRUE)

load(rdata)
if (!exists("df_filtered")) stop("Object df_filtered not found in ", rdata)
names(df_filtered) <- tolower(names(df_filtered))

message("Connecting to ODBC DSN: ", dsn)
con <- DBI::dbConnect(odbc::odbc(), dsn = dsn)
on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)
subtype <- DBI::dbGetQuery(con, "SELECT * FROM clinical.subtypes")
names(subtype) <- tolower(names(subtype))

if (!"case_barcode" %in% names(df_filtered)) stop("df_filtered missing case_barcode")
if (!"case_barcode" %in% names(subtype)) stop("clinical.subtypes missing case_barcode")
if (!"idh_codel_subtype" %in% names(subtype)) stop("clinical.subtypes missing idh_codel_subtype")

df_merged <- merge(df_filtered, subtype, by = "case_barcode", all.x = TRUE)
df_sub <- df_merged[df_merged$idh_codel_subtype == subtype_value, , drop = FALSE]
message("Rows after subtype filter [", subtype_value, "]: ", nrow(df_sub))

required <- c("tumor_barcode_a", "tumor_barcode_b", "variant_type", "chrom", "pos", "ref", "alt", "mutect2_call_a", "mutect2_call_b")
missing <- setdiff(required, names(df_sub))
if (length(missing) > 0) stop("Filtered data missing required columns: ", paste(missing, collapse = ", "))

make_arm <- function(dat, sample_col, call_col) {
  keep <- !is.na(dat[[sample_col]]) & !is.na(dat[[call_col]]) & dat[[call_col]] == 1
  out <- dat[keep, c(sample_col, "variant_type", "chrom", "pos", "ref", "alt"), drop = FALSE]
  names(out) <- c("Sample", "mut_type", "chrom", "pos", "ref", "alt")
  out
}

df_long <- bind_rows(
  make_arm(df_sub, "tumor_barcode_a", "mutect2_call_a"),
  make_arm(df_sub, "tumor_barcode_b", "mutect2_call_b")
) %>% distinct()

if (nrow(df_long) == 0) stop("No variants remained after mutect2_call_a/b and subtype filtering.")

extract_pos <- function(x, which_num = 1L) {
  vals <- regmatches(as.character(x), gregexpr("[0-9]+", as.character(x)))
  vapply(vals, function(z) {
    if (length(z) >= which_num) as.integer(z[[which_num]]) else NA_integer_
  }, integer(1))
}

df_long$pos_start <- extract_pos(df_long$pos, 1L)
df_long$pos_end <- extract_pos(df_long$pos, 2L)
df_long$pos_end[is.na(df_long$pos_end)] <- df_long$pos_start[is.na(df_long$pos_end)]
if (any(is.na(df_long$pos_start))) stop("Could not parse numeric position from some pos values.")

df_out <- df_long %>%
  transmute(
    Project = project_name,
    Sample = .data$Sample,
    ID = ".",
    Genome = genome,
    mut_type = .data$mut_type,
    chrom = .data$chrom,
    pos_start = .data$pos_start,
    pos_end = .data$pos_end,
    ref = .data$ref,
    alt = .data$alt,
    Type = "SOMATIC"
  ) %>% distinct()

split_df <- split(df_out, df_out$Sample)
summary <- data.frame(Sample = character(), data_type = character(), n_variants = integer(), file = character())

for (sample_name in names(split_df)) {
  dat <- split_df[[sample_name]]
  safe_sample <- gsub("[^A-Za-z0-9_.-]", "_", sample_name)
  all_file <- file.path(outdir, "all", paste0(safe_sample, ".txt"))
  write.table(dat, file = all_file, sep = "\t", row.names = FALSE, quote = FALSE)

  dtype <- if (grepl("WGS", sample_name, ignore.case = TRUE)) {
    "wgs"
  } else if (grepl("WXS|WES|EXOME", sample_name, ignore.case = TRUE)) {
    "wxs"
  } else {
    "unknown"
  }
  typed_file <- file.path(outdir, dtype, paste0(safe_sample, ".txt"))
  write.table(dat, file = typed_file, sep = "\t", row.names = FALSE, quote = FALSE)
  summary <- rbind(summary, data.frame(Sample = sample_name, data_type = dtype, n_variants = nrow(dat), file = typed_file))
  message("Wrote ", dtype, " sample: ", sample_name, " -> ", typed_file)
}

write.table(summary, file = summary_path, sep = "\t", row.names = FALSE, quote = FALSE)
message("Saved input summary: ", summary_path)
RSCRIPT

  cat > "$SP_PY" <<'PYTHON'
#!/usr/bin/env python3
import argparse
import glob
import os
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Run SigProfilerMatrixGenerator and SigProfilerExtractor on one input directory.")
    p.add_argument("--input-dir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--project", default="codel")
    p.add_argument("--genome", default="GRCh37")
    p.add_argument("--exome", action="store_true")
    p.add_argument("--sig-types", default="SBS,DBS,ID")
    p.add_argument("--min-signatures", type=int, default=1)
    p.add_argument("--max-signatures", type=int, default=20)
    p.add_argument("--cushion", type=int, default=100)
    p.add_argument("--plot", default="true", choices=["true", "false"])
    p.add_argument("--skip-matrix", action="store_true")
    p.add_argument("--skip-extract", action="store_true")
    p.add_argument("--force", action="store_true")
    return p.parse_args()


def matrix_candidates(input_dir: Path, project: str, sig_type: str):
    suffix = {"SBS": "SBS96", "DBS": "DBS78", "ID": "ID83"}[sig_type]
    subdir = {"SBS": "SBS", "DBS": "DBS", "ID": "ID"}[sig_type]
    return [
        input_dir / "output" / subdir / f"{project}.{suffix}.all",
        input_dir / "output" / subdir / f"{input_dir.name}.{suffix}.all",
    ]


def find_matrix(input_dir: Path, project: str, sig_type: str):
    for candidate in matrix_candidates(input_dir, project, sig_type):
        if candidate.exists() and candidate.stat().st_size > 0:
            return candidate
    suffix = {"SBS": "SBS96", "DBS": "DBS78", "ID": "ID83"}[sig_type]
    matches = sorted(glob.glob(str(input_dir / "output" / "**" / f"*.{suffix}.all"), recursive=True))
    if matches:
        return Path(matches[0])
    return None


def main():
    args = parse_args()
    input_dir = Path(args.input_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    txt_files = list(input_dir.glob("*.txt"))
    if not txt_files:
        raise SystemExit(f"No .txt SigProfiler input files found in {input_dir}")

    sig_types = [x.strip().upper() for x in args.sig_types.split(",") if x.strip()]
    valid = {"SBS", "DBS", "ID"}
    bad = sorted(set(sig_types) - valid)
    if bad:
        raise SystemExit(f"Unsupported signature type(s): {bad}")

    if not args.skip_matrix:
        print(f"[SigProfiler] Creating matrices for {input_dir}", flush=True)
        from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
        matGen.SigProfilerMatrixGeneratorFunc(
            args.project,
            args.genome,
            str(input_dir),
            plot=(args.plot == "true"),
            exome=args.exome,
            bed_file=None,
            chrom_based=False,
            tsb_stat=False,
            seqInfo=False,
            cushion=args.cushion,
        )

    if args.skip_extract:
        print("[SigProfiler] Skipping extraction by request.", flush=True)
        return

    from SigProfilerExtractor import sigpro as sig

    for sig_type in sig_types:
        matrix = find_matrix(input_dir, args.project, sig_type)
        if matrix is None:
            raise SystemExit(f"Could not find {sig_type} matrix under {input_dir}/output")

        sig_outdir = outdir / sig_type
        marker = sig_outdir / f".{sig_type}.extract.done"
        if marker.exists() and not args.force:
            print(f"[SigProfiler] Skipping completed extraction: {sig_type} ({marker})", flush=True)
            continue

        sig_outdir.mkdir(parents=True, exist_ok=True)
        print(f"[SigProfiler] Extracting {sig_type} signatures from matrix: {matrix}", flush=True)
        sig.sigProfilerExtractor(
            "matrix",
            str(sig_outdir),
            str(matrix),
            reference_genome=args.genome,
            exome=args.exome,
            minimum_signatures=args.min_signatures,
            maximum_signatures=args.max_signatures,
        )
        marker.write_text("done\n")
        print(f"[SigProfiler] Completed {sig_type}: {sig_outdir}", flush=True)


if __name__ == "__main__":
    main()
PYTHON
  chmod +x "$SP_PY"

  cat > "$PAL_R" <<'RSCRIPT'
suppressPackageStartupMessages({
  library(Palimpsest)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) stop("Missing value for ", flag)
  args[[hit + 1]]
}

input_dir <- normalizePath(get_arg("--input-dir"), mustWork = TRUE)
outdir <- normalizePath(get_arg("--outdir"), mustWork = FALSE)
reference_rdata <- normalizePath(get_arg("--reference-rdata"), mustWork = TRUE)
sig_types <- toupper(strsplit(get_arg("--sig-types", "SBS,DBS,ID"), ",")[[1]])
sbs_ref_object <- get_arg("--sbs-ref-object", "SBS_cosmic3.2_sigs")
dbs_ref_object <- get_arg("--dbs-ref-object", "dbs_ref_codel")
id_ref_object <- get_arg("--id-ref-object", "ID_cosmic")
id83_reference_file <- get_arg("--id83-reference-file", "")
ref_fasta <- get_arg("--ref-fasta", "")
sample_info_object <- get_arg("--sample-info-object", "codel_info")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("Loading reference RData: ", reference_rdata)
load(reference_rdata)

# Load the hg19 package explicitly if available. Palimpsest histories used BSgenome.Hsapiens.UCSC.hg19.
if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
  suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
  ref_genome <- BSgenome.Hsapiens.UCSC.hg19
} else {
  stop("Package BSgenome.Hsapiens.UCSC.hg19 is required for this GRCh37/hg19 Palimpsest workflow.")
}

read_sigprofiler_inputs <- function(input_dir) {
  files <- list.files(input_dir, pattern = "\\.txt$", full.names = TRUE)
  if (length(files) == 0) stop("No .txt input files found in ", input_dir)
  dfs <- lapply(files, function(f) {
    x <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
    x$.source_file <- basename(f)
    x
  })
  dat <- bind_rows(dfs)
  names(dat) <- tolower(names(dat))
  required <- c("sample", "mut_type", "chrom", "pos_start", "ref", "alt")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) stop("Input files missing required columns: ", paste(missing, collapse = ", "))
  dat %>% transmute(
    Sample = .data$sample,
    Type = toupper(.data$mut_type),
    CHROM = as.character(.data$chrom),
    POS = as.integer(.data$pos_start),
    REF = as.character(.data$ref),
    ALT = as.character(.data$alt)
  ) %>% distinct()
}

normalize_vcf <- function(vcf) {
  vcf$Type[vcf$Type == "SNP"] <- "SNV"
  vcf$Type[vcf$Type == "DNP"] <- "DBS"
  vcf$Type[vcf$Type == "DNV"] <- "DBS"
  vcf$Type[vcf$Type == "INS"] <- "INS"
  vcf$Type[vcf$Type == "DEL"] <- "DEL"
  vcf$CHROM <- sub("^chr", "", vcf$CHROM, ignore.case = TRUE)
  vcf$CHROM[vcf$CHROM == "23"] <- "X"
  vcf$CHROM[vcf$CHROM == "24"] <- "Y"
  vcf
}

read_id83_reference <- function(path) {
  raw <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"MutationType" %in% names(raw)) stop("ID83 reference file must contain MutationType column.")
  parts <- strsplit(raw$MutationType, ":")
  if (any(lengths(parts) < 4)) stop("Unexpected MutationType format in ID83 reference file.")
  p <- do.call(rbind, lapply(parts, function(x) x[1:4]))
  p <- as.data.frame(p, stringsAsFactors = FALSE)
  names(p) <- c("len", "event", "context", "repeat_len")
  p$context[p$context == "R"] <- "repeats"
  p$context[p$context == "M"] <- "MH"
  p$len[p$len == "5"] <- "5+"
  p$repeat_len[p$repeat_len == "5"] <- "5+"
  p$event[p$event %in% c("Del", "DEL")] <- "DEL"
  p$event[p$event %in% c("Ins", "INS")] <- "INS"
  category <- paste(p$event, p$context, p$len, p$repeat_len, sep = "_")
  sig_cols <- setdiff(names(raw), "MutationType")
  mat <- as.matrix(raw[, sig_cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- category
  t(mat)
}

get_reference <- function(sig_type) {
  if (sig_type == "SBS") {
    if (!exists(sbs_ref_object)) stop("SBS reference object not found: ", sbs_ref_object)
    return(get(sbs_ref_object))
  }
  if (sig_type == "DBS") {
    if (!exists(dbs_ref_object)) stop("DBS reference object not found: ", dbs_ref_object)
    return(get(dbs_ref_object))
  }
  if (sig_type == "ID") {
    if (nzchar(id83_reference_file)) return(read_id83_reference(id83_reference_file))
    if (!exists(id_ref_object)) stop("ID reference object not found: ", id_ref_object)
    return(get(id_ref_object))
  }
  stop("Unsupported signature type: ", sig_type)
}

subset_for_sig_type <- function(vcf, sig_type) {
  if (sig_type == "SBS") return(vcf %>% filter(.data$Type %in% c("SNV", "SNP")))
  if (sig_type == "DBS") return(vcf %>% filter(.data$Type %in% c("DBS", "DNP", "DNV")))
  if (sig_type == "ID") return(vcf %>% filter(grepl("INS|DEL|INDEL", .data$Type)))
  stop("Unsupported signature type: ", sig_type)
}

annotate_for_sig_type <- function(vcf, sig_type) {
  if (nrow(vcf) == 0) return(vcf)
  if (sig_type == "SBS") {
    return(annotate_VCF(vcf = vcf, ref_genome = ref_genome, add_ID_cats = FALSE, add_DBS_cats = FALSE))
  }
  if (sig_type == "DBS") {
    return(annotate_VCF(vcf = vcf, ref_genome = ref_genome, add_strand_and_SBS_cats = FALSE, add_ID_cats = FALSE, add_DBS_cats = TRUE))
  }
  if (sig_type == "ID") {
    if (nzchar(ref_fasta)) {
      return(annotate_VCF(vcf = vcf, ref_genome = ref_genome, ref_fasta = ref_fasta, add_strand_and_SBS_cats = FALSE, add_DBS_cats = FALSE, add_ID_cats = TRUE))
    }
    return(annotate_VCF(vcf = vcf, ref_genome = ref_genome, add_strand_and_SBS_cats = FALSE, add_DBS_cats = FALSE, add_ID_cats = TRUE))
  }
  stop("Unsupported signature type: ", sig_type)
}

write_any_table <- function(obj, prefix) {
  saveRDS(obj, paste0(prefix, ".rds"))
  if (is.data.frame(obj) || is.matrix(obj)) {
    write.table(as.data.frame(obj), paste0(prefix, ".tsv"), sep = "\t", quote = FALSE, row.names = TRUE)
  } else if (is.list(obj)) {
    for (nm in names(obj)) {
      if (is.data.frame(obj[[nm]]) || is.matrix(obj[[nm]])) {
        write.table(as.data.frame(obj[[nm]]), paste0(prefix, ".", nm, ".tsv"), sep = "\t", quote = FALSE, row.names = TRUE)
      }
    }
  }
}

make_fraction_vcf <- function(annotated_vcf) {
  if (!exists(sample_info_object)) return(NULL)
  info <- get(sample_info_object)
  names(info) <- tolower(names(info))
  if (!all(c("sample", "case_barcode") %in% names(info))) {
    message("Sample info object exists but lacks sample/case_barcode columns; skipping shared/private fraction mode.")
    return(NULL)
  }
  keep_cols <- unique(c("sample", "case_barcode", "surgery_number", "sample_type"))
  keep_cols <- keep_cols[keep_cols %in% names(info)]
  info <- info[, keep_cols, drop = FALSE]
  names(info)[names(info) == "sample"] <- "Sample"

  merged <- merge(annotated_vcf, info, by = "Sample", all.x = TRUE)
  if (!"case_barcode" %in% names(merged)) return(NULL)
  tmp <- merged %>%
    group_by(.data$case_barcode, .data$Type, .data$CHROM, .data$POS, .data$REF, .data$ALT) %>%
    mutate(variant_fraction = ifelse(n_distinct(.data$Sample) > 1, "shared", "private")) %>%
    ungroup()
  tmp$Sample_original <- tmp$Sample
  tmp$Sample <- ifelse(tmp$variant_fraction == "shared" & !is.na(tmp$case_barcode), paste0(tmp$case_barcode, "-shared"), tmp$Sample)
  tmp
}

run_one <- function(vcf_all, sig_type) {
  sig_type <- toupper(sig_type)
  sig_out <- file.path(outdir, sig_type)
  dir.create(sig_out, showWarnings = FALSE, recursive = TRUE)
  marker <- file.path(sig_out, paste0(".", sig_type, ".palimpsest.done"))
  if (file.exists(marker)) {
    message("Skipping already completed Palimpsest step: ", sig_type)
    return(invisible(NULL))
  }

  vcf <- subset_for_sig_type(vcf_all, sig_type)
  if (nrow(vcf) == 0) {
    message("No ", sig_type, " variants in input; skipping.")
    writeLines("no variants", marker)
    return(invisible(NULL))
  }

  message("Annotating ", sig_type, " variants: ", nrow(vcf))
  ann <- annotate_for_sig_type(vcf, sig_type)
  write.table(ann, file.path(sig_out, paste0(sig_type, "_annotated_vcf.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  ref <- get_reference(sig_type)
  input_sample <- palimpsest_input(vcf = ann, Type = sig_type)
  saveRDS(input_sample, file.path(sig_out, paste0(sig_type, "_palimpsest_input_sample.rds")))
  exp_sample <- deconvolution_fit(input_matrices = input_sample, input_signatures = ref, resdir = sig_out, input_vcf = ann, doplot = FALSE)
  write_any_table(exp_sample, file.path(sig_out, paste0(sig_type, "_signature_exposure_sample")))

  origins_sample <- tryCatch(
    signature_origins(input = ann, Type = sig_type, input_signatures = ref, signature_contribution = exp_sample),
    error = function(e) {
      message("signature_origins failed for sample-level ", sig_type, ": ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(origins_sample)) {
    write.table(origins_sample, file.path(sig_out, paste0(sig_type, "_mutation_origins_sample.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  }

  frac <- make_fraction_vcf(ann)
  if (!is.null(frac)) {
    input_frac <- palimpsest_input(vcf = frac, Type = sig_type)
    saveRDS(input_frac, file.path(sig_out, paste0(sig_type, "_palimpsest_input_fraction.rds")))
    exp_frac <- deconvolution_fit(input_matrices = input_frac, input_signatures = ref, resdir = sig_out, input_vcf = frac, doplot = FALSE)
    write_any_table(exp_frac, file.path(sig_out, paste0(sig_type, "_signature_exposure_fraction")))
    origins_frac <- tryCatch(
      signature_origins(input = frac, Type = sig_type, input_signatures = ref, signature_contribution = exp_frac),
      error = function(e) {
        message("signature_origins failed for fraction-level ", sig_type, ": ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(origins_frac)) {
      write.table(origins_frac, file.path(sig_out, paste0(sig_type, "_mutation_origins_fraction.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }

  save.image(file.path(sig_out, paste0(sig_type, "_palimpsest_workspace.RData")))
  writeLines(paste0("completed_at=", Sys.time()), marker)
  message("Completed Palimpsest ", sig_type, ": ", sig_out)
}

vcf_all <- read_sigprofiler_inputs(input_dir) %>% normalize_vcf()
write.table(vcf_all, file.path(outdir, "palimpsest_input_variants_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (sig_type in sig_types) run_one(vcf_all, sig_type)
RSCRIPT
}

write_helpers

if has_step filter || has_step input; then
  if [[ -z "$SQL" ]]; then
    fail "--sql is required for filter/input steps."
  fi
  [[ -s "$SQL" ]] || fail "SQL file not found or empty: $SQL"
fi

if has_step palimpsest; then
  [[ -n "$REFERENCE_RDATA" ]] || fail "--reference-rdata is required for palimpsest step."
  [[ -s "$REFERENCE_RDATA" ]] || fail "Reference RData not found or empty: $REFERENCE_RDATA"
fi

log "Pipeline workdir: $WORKDIR"
log "Requested steps: $STEPS"
log "Data types: $DATA_TYPES"
log "Signature types: $SIG_TYPES"

# Step 1: repeat region filtering
if has_step filter; then
  if is_done filter && [[ -s "$FILTER_RDATA" ]]; then
    log "Skipping filter step; marker and output exist. Use --force to rerun."
  else
    log "Running Step 1: repeat-region filtering"
    activate_env "$R_ENV"
    run_cmd Rscript "$SCRIPT_DIR/01_filter_repeat_regions.R" \
      --workdir "$FILTER_DIR" \
      --sql "$SQL" \
      --dsn "$DSN" \
      --out-rdata "$FILTER_RDATA" \
      --stats "$FILTER_STATS" \
      > "$LOG_DIR/01_filter_repeat_regions.log" 2>&1
    [[ "$DRY_RUN" == "1" || -s "$FILTER_RDATA" ]] || fail "Filter RData was not created: $FILTER_RDATA"
    mark_done filter
  fi
fi

# Step 2: SigProfiler input files
if has_step input; then
  if is_done input && [[ -s "$INPUT_SUMMARY" ]]; then
    log "Skipping SigProfiler input step; marker and summary exist. Use --force to rerun."
  else
    [[ -s "$FILTER_RDATA" ]] || fail "Missing filter RData for input step: $FILTER_RDATA. Run --steps filter,input or provide existing output."
    log "Running Step 2: SigProfiler input generation + WGS/WXS split"
    activate_env "$R_ENV"
    run_cmd Rscript "$SCRIPT_DIR/02_make_sigprofiler_input.R" \
      --filter-rdata "$FILTER_RDATA" \
      --outdir "$INPUT_DIR" \
      --dsn "$DSN" \
      --subtype "$SUBTYPE" \
      --project-name "$PROJECT_NAME" \
      --genome "$GENOME" \
      --summary "$INPUT_SUMMARY" \
      > "$LOG_DIR/02_make_sigprofiler_input.log" 2>&1
    [[ "$DRY_RUN" == "1" || -s "$INPUT_SUMMARY" ]] || fail "Input summary was not created: $INPUT_SUMMARY"
    mark_done input
  fi
fi

# Step 3: matrix generation and/or de novo extraction
if has_step matrix || has_step extract; then
  activate_env "$SIGPROFILER_ENV"
  IFS=',' read -r -a DTYPES <<< "$DATA_TYPES"
  for dtype_raw in "${DTYPES[@]}"; do
    dtype=$(echo "$dtype_raw" | tr '[:upper:]' '[:lower:]' | xargs)
    [[ "$dtype" == "wgs" || "$dtype" == "wxs" ]] || fail "Unsupported data type: $dtype"
    dtype_input="$INPUT_DIR/$dtype"
    [[ -d "$dtype_input" ]] || fail "Missing input directory: $dtype_input"
    if ! compgen -G "$dtype_input/*.txt" >/dev/null; then
      log "No .txt files found for $dtype in $dtype_input; skipping matrix/extract for this data type."
      continue
    fi

    marker_name="sigprofiler_${dtype}"
    if has_step matrix && has_step extract; then
      if is_done "$marker_name"; then
        log "Skipping SigProfiler matrix+extract for $dtype; marker exists. Use --force to rerun."
        continue
      fi
    fi

    exome_flag=()
    if [[ "$dtype" == "wxs" ]]; then
      exome_flag=(--exome)
    fi

    skip_matrix_flag=()
    skip_extract_flag=()
    if ! has_step matrix; then skip_matrix_flag=(--skip-matrix); fi
    if ! has_step extract; then skip_extract_flag=(--skip-extract); fi
    force_flag=()
    if [[ "$FORCE" == "1" ]]; then force_flag=(--force); fi

    log "Running Step 3 for $dtype: matrix=$(has_step matrix && echo yes || echo no), extract=$(has_step extract && echo yes || echo no)"
    run_cmd python "$SCRIPT_DIR/03_run_sigprofiler.py" \
      --input-dir "$dtype_input" \
      --outdir "$SP_OUTDIR/$dtype" \
      --project "$MATRIX_PROJECT" \
      --genome "$GENOME" \
      --sig-types "$SIG_TYPES" \
      --min-signatures "$MIN_SIGNATURES" \
      --max-signatures "$MAX_SIGNATURES" \
      --cushion "$CUSHION" \
      --plot "$PLOT_MATRIX" \
      "${exome_flag[@]}" \
      "${skip_matrix_flag[@]}" \
      "${skip_extract_flag[@]}" \
      "${force_flag[@]}" \
      > "$LOG_DIR/03_sigprofiler_${dtype}.log" 2>&1
    mark_done "$marker_name"
  done
fi

# Step 4: Palimpsest reference deconvolution
if has_step palimpsest; then
  activate_env "$R_ENV"
  IFS=',' read -r -a DTYPES <<< "$DATA_TYPES"
  for dtype_raw in "${DTYPES[@]}"; do
    dtype=$(echo "$dtype_raw" | tr '[:upper:]' '[:lower:]' | xargs)
    dtype_input="$INPUT_DIR/$dtype"
    [[ -d "$dtype_input" ]] || fail "Missing input directory: $dtype_input"
    if ! compgen -G "$dtype_input/*.txt" >/dev/null; then
      log "No .txt files found for $dtype in $dtype_input; skipping Palimpsest for this data type."
      continue
    fi
    marker_name="palimpsest_${dtype}"
    if is_done "$marker_name"; then
      log "Skipping Palimpsest deconvolution for $dtype; marker exists. Use --force to rerun."
      continue
    fi

    id83_arg=()
    if [[ -n "$ID83_REFERENCE_FILE" ]]; then id83_arg=(--id83-reference-file "$ID83_REFERENCE_FILE"); fi
    ref_fasta_arg=()
    if [[ -n "$REF_FASTA" ]]; then ref_fasta_arg=(--ref-fasta "$REF_FASTA"); fi

    log "Running Step 4: Palimpsest deconvolution for $dtype"
    run_cmd Rscript "$SCRIPT_DIR/04_run_palimpsest_deconvolution.R" \
      --input-dir "$dtype_input" \
      --outdir "$PAL_OUTDIR/$dtype" \
      --reference-rdata "$REFERENCE_RDATA" \
      --sig-types "$SIG_TYPES" \
      --sbs-ref-object "$SBS_REF_OBJECT" \
      --dbs-ref-object "$DBS_REF_OBJECT" \
      --id-ref-object "$ID_REF_OBJECT" \
      --sample-info-object "$SAMPLE_INFO_OBJECT" \
      "${id83_arg[@]}" \
      "${ref_fasta_arg[@]}" \
      > "$LOG_DIR/04_palimpsest_${dtype}.log" 2>&1
    mark_done "$marker_name"
  done
fi

log "Pipeline complete."
log "Outputs:"
log "  Filter RData:      $FILTER_RDATA"
log "  SP input:          $INPUT_DIR"
log "  SigProfiler out:   $SP_OUTDIR"
log "  Palimpsest out:    $PAL_OUTDIR"
log "  Logs:              $LOG_DIR"
