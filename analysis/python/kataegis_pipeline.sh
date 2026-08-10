#!/usr/bin/env bash
#SBATCH --job-name=kataegis_pipeline
#SBATCH --partition=day
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# Coherent kataegis / clustered mutation pipeline for GLASS5-style data:
#   1) RepeatMasker/ENCODE repeat-region filtering from GLASS5 using final_query.sql
#   2) Create SigProfiler-style per-sample mutation input files
#   3) Run SigProfilerSimulator + SigProfilerClusters clustered mutation analysis
#   4) Generate kataegis event count and detail tables
#
# The script writes helper R/Python scripts into WORKDIR/_pipeline_scripts and uses
# marker files in WORKDIR/.pipeline_markers so steps are resume-aware.

set -Eeuo pipefail
IFS=$'\n\t'

SCRIPT_NAME=$(basename "$0")

usage() {
  cat <<'EOF'
Usage:
  kataegis_pipeline.sh --workdir WORKDIR --sql final_query.sql [options]

Required for filter/input:
  --workdir PATH                  Main pipeline working directory.
  --sql PATH                      SQL query used for repeat-region filtering.

Common options:
  --steps LIST                    Comma list: filter,input,cluster,count,all [default: all]
  --data-types LIST               Comma list: wgs,wxs,unknown [default: wgs]
  --dsn NAME                      ODBC DSN for GLASS5 database [default: glass5]
  --subtype VALUE                 clinical.subtypes idh_codel_subtype filter [default: IDHmut-codel]
  --project-name VALUE            Project/cohort label [default: GLASS5]
  --genome VALUE                  Reference genome label [default: GRCh37]
  --r-env NAME                    Conda environment with R/tidyverse/DBI/odbc [default: palimpsest]
  --cluster-env NAME              Conda environment with SigProfilerSimulator/Clusters [default: sigprofiler]
  --miniconda-module NAME         Module to load before conda activation [default: miniconda]
  --no-module-load                Do not run module load before conda activate.
  --no-conda                      Do not activate conda environments.
  --force                         Rerun completed stages and overwrite marker files.
  --dry-run                       Print planned actions without running them.

Clustered mutation options:
  --sim-contexts LIST             Comma list for SigProfilerSimulator contexts [default: 96,ID]
  --simulations INT               Number of simulations [default: 100]
  --cluster-context VALUE         SigProfilerClusters context [default: 96]
  --cluster-subcontexts LIST      Comma list passed as sub-contexts [default: 6144]
  --window-size INT               SigProfilerClusters windowSize [default: 1000000]
  --max-cpu INT                   Max CPU for SigProfilerClusters [default: SLURM_CPUS_PER_TASK or 8]

Output layout:
  WORKDIR/filter_step/filtered_repeat_region_variants.RData
  WORKDIR/filter_step/filter_repeat_region_stats.tsv
  WORKDIR/sp_input/WGS/input/*.txt
  WORKDIR/sp_input/WGS/output/vcf_files_corrected/GLASS5_clustered/GLASS5_clusters_of_clusters_imd_edit.txt
  WORKDIR/r_downstream/WGS/kataegis_event_count_all_wgs.csv
  WORKDIR/r_downstream/WGS/kataegis_details_all_wgs.csv

Examples:
  # Full WGS run
  sbatch kataegis_pipeline.sh \
    --workdir /path/to/kataegis_run \
    --sql /path/to/final_query.sql \
    --steps all \
    --data-types wgs \
    --r-env palimpsest \
    --cluster-env sigprofiler

  # Run only Step 4 after SigProfilerClusters finished
  sbatch kataegis_pipeline.sh \
    --workdir /path/to/kataegis_run \
    --steps count \
    --data-types wgs
EOF
}

WORKDIR=""
SQL=""
STEPS="all"
DATA_TYPES="wgs"
DSN="glass5"
SUBTYPE="IDHmut-codel"
PROJECT_NAME="GLASS5"
GENOME="GRCh37"
R_ENV="palimpsest"
CLUSTER_ENV="sigprofiler"
MINICONDA_MODULE="miniconda"
LOAD_MODULE="1"
USE_CONDA="1"
FORCE="0"
DRY_RUN="0"
SIM_CONTEXTS="96,ID"
SIMULATIONS="100"
CLUSTER_CONTEXT="96"
CLUSTER_SUBCONTEXTS="6144"
WINDOW_SIZE="1000000"
MAX_CPU="${SLURM_CPUS_PER_TASK:-8}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workdir) WORKDIR="$2"; shift 2 ;;
    --sql) SQL="$2"; shift 2 ;;
    --steps) STEPS="$2"; shift 2 ;;
    --data-types) DATA_TYPES="$2"; shift 2 ;;
    --dsn) DSN="$2"; shift 2 ;;
    --subtype) SUBTYPE="$2"; shift 2 ;;
    --project-name) PROJECT_NAME="$2"; shift 2 ;;
    --genome) GENOME="$2"; shift 2 ;;
    --r-env) R_ENV="$2"; shift 2 ;;
    --cluster-env) CLUSTER_ENV="$2"; shift 2 ;;
    --miniconda-module) MINICONDA_MODULE="$2"; shift 2 ;;
    --no-module-load) LOAD_MODULE="0"; shift ;;
    --no-conda) USE_CONDA="0"; shift ;;
    --force) FORCE="1"; shift ;;
    --dry-run) DRY_RUN="1"; shift ;;
    --sim-contexts) SIM_CONTEXTS="$2"; shift 2 ;;
    --simulations) SIMULATIONS="$2"; shift 2 ;;
    --cluster-context) CLUSTER_CONTEXT="$2"; shift 2 ;;
    --cluster-subcontexts) CLUSTER_SUBCONTEXTS="$2"; shift 2 ;;
    --window-size) WINDOW_SIZE="$2"; shift 2 ;;
    --max-cpu) MAX_CPU="$2"; shift 2 ;;
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
INPUT_ROOT="$WORKDIR/sp_input"
DOWNSTREAM_ROOT="$WORKDIR/r_downstream"
MARKER_DIR="$WORKDIR/.pipeline_markers"
FILTER_RDATA="$FILTER_DIR/filtered_repeat_region_variants.RData"
FILTER_STATS="$FILTER_DIR/filter_repeat_region_stats.tsv"
INPUT_SUMMARY="$INPUT_ROOT/kataegis_input_summary.tsv"

mkdir -p "$WORKDIR" "$SCRIPT_DIR" "$LOG_DIR" "$FILTER_DIR" "$INPUT_ROOT" "$DOWNSTREAM_ROOT" "$MARKER_DIR"

if [[ -n "$SQL" ]]; then
  SQL=$(readlink -m "$SQL")
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

uppercase_dtype() {
  local dtype="$1"
  case "${dtype,,}" in
    wgs) echo "WGS" ;;
    wxs|wes|exome) echo "WXS" ;;
    unknown) echo "UNKNOWN" ;;
    *) fail "Unsupported data type: $dtype. Use wgs,wxs,unknown." ;;
  esac
}

lowercase_dtype() {
  local dtype="$1"
  case "${dtype,,}" in
    wgs) echo "wgs" ;;
    wxs|wes|exome) echo "wxs" ;;
    unknown) echo "unknown" ;;
    *) fail "Unsupported data type: $dtype. Use wgs,wxs,unknown." ;;
  esac
}

activate_env() {
  local env_name="$1"
  if [[ "$USE_CONDA" != "1" ]]; then
    return 0
  fi
  if [[ "$LOAD_MODULE" == "1" ]]; then
    if command -v module >/dev/null 2>&1 || type module >/dev/null 2>&1; then
      module load "$MINICONDA_MODULE" || true
    fi
  fi
  if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
    conda activate "$env_name"
  else
    fail "conda not found. Use --no-conda if the required environment is already active."
  fi
}

write_helpers() {
  local FILTER_R="$SCRIPT_DIR/01_filter_repeat_regions.R"
  local INPUT_R="$SCRIPT_DIR/02_make_kataegis_input.R"
  local CLUSTER_PY="$SCRIPT_DIR/03_run_sigprofiler_clusters.py"
  local COUNT_R="$SCRIPT_DIR/04_make_kataegis_count_table.R"

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
outroot <- normalizePath(get_arg("--outroot"), mustWork = FALSE)
dsn <- get_arg("--dsn", "glass5")
subtype_value <- get_arg("--subtype", "IDHmut-codel")
project_name <- get_arg("--project-name", "GLASS5")
genome <- get_arg("--genome", "GRCh37")
summary_path <- normalizePath(get_arg("--summary", file.path(outroot, "kataegis_input_summary.tsv")), mustWork = FALSE)

for (x in c("ALL", "WGS", "WXS", "UNKNOWN")) {
  dir.create(file.path(outroot, x, "input"), showWarnings = FALSE, recursive = TRUE)
}
dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)

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

required <- c(
  "tumor_barcode_a", "tumor_barcode_b", "variant_type", "chrom", "pos", "ref", "alt",
  "mutect2_call_a", "mutect2_call_b"
)
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

  dtype <- if (grepl("WGS", sample_name, ignore.case = TRUE)) {
    "WGS"
  } else if (grepl("WXS|WES|EXOME", sample_name, ignore.case = TRUE)) {
    "WXS"
  } else {
    "UNKNOWN"
  }

  all_file <- file.path(outroot, "ALL", "input", paste0(safe_sample, ".txt"))
  typed_file <- file.path(outroot, dtype, "input", paste0(safe_sample, ".txt"))

  write.table(dat, file = all_file, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  write.table(dat, file = typed_file, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

  summary <- rbind(summary, data.frame(Sample = sample_name, data_type = dtype, n_variants = nrow(dat), file = typed_file))
  message("Wrote ", dtype, " sample: ", sample_name, " -> ", typed_file)
}

write.table(summary, file = summary_path, sep = "\t", row.names = FALSE, quote = FALSE)
message("Saved input summary: ", summary_path)
RSCRIPT

  cat > "$CLUSTER_PY" <<'PYTHON'
#!/usr/bin/env python3
import argparse
import glob
import os
import sys
from pathlib import Path


def split_csv(value):
    return [x.strip() for x in value.split(",") if x.strip()]


def parse_args():
    p = argparse.ArgumentParser(description="Run SigProfilerSimulator and SigProfilerClusters for kataegis/clustered mutation analysis.")
    p.add_argument("--project", required=True)
    p.add_argument("--input-root", required=True, help="Directory containing input/*.txt and receiving output/.")
    p.add_argument("--genome", default="GRCh37")
    p.add_argument("--sim-contexts", default="96,ID")
    p.add_argument("--simulations", type=int, default=100)
    p.add_argument("--cluster-context", default="96")
    p.add_argument("--cluster-subcontexts", default="6144")
    p.add_argument("--window-size", type=int, default=1000000)
    p.add_argument("--max-cpu", type=int, default=8)
    return p.parse_args()


def import_tools():
    try:
        from SigProfilerSimulator import SigProfilerSimulator as sigSim
    except Exception as exc:
        raise RuntimeError(
            "Could not import SigProfilerSimulator. Make sure the active environment contains SigProfilerSimulator."
        ) from exc

    try:
        from SigProfilerClusters import SigProfilerClusters as hp
    except Exception as exc:
        raise RuntimeError(
            "Could not import SigProfilerClusters. Make sure the active environment contains SigProfilerClusters."
        ) from exc

    return sigSim, hp


def main():
    args = parse_args()
    input_root = Path(args.input_root).resolve()
    input_subdir = input_root / "input"

    if not input_subdir.is_dir():
        raise FileNotFoundError(f"Expected input directory not found: {input_subdir}")

    input_files = sorted(input_subdir.glob("*.txt"))
    if not input_files:
        raise FileNotFoundError(f"No per-sample .txt files found in: {input_subdir}")

    sigSim, hp = import_tools()

    sim_contexts = split_csv(args.sim_contexts)
    cluster_subcontexts = split_csv(args.cluster_subcontexts)

    print(f"Input root: {input_root}", flush=True)
    print(f"Samples: {len(input_files)}", flush=True)
    print(f"Running SigProfilerSimulator contexts={sim_contexts}, simulations={args.simulations}", flush=True)

    sigSim.SigProfilerSimulator(
        args.project,
        str(input_root) + os.sep,
        args.genome,
        contexts=sim_contexts,
        chrom_based=True,
        simulations=args.simulations,
    )

    print(
        f"Running SigProfilerClusters context={args.cluster_context}, "
        f"subcontexts={cluster_subcontexts}, max_cpu={args.max_cpu}",
        flush=True,
    )

    hp.analysis(
        args.project,
        args.genome,
        args.cluster_context,
        cluster_subcontexts,
        str(input_root) + os.sep,
        analysis="all",
        sortSims=True,
        subClassify=True,
        correction=True,
        calculateIMD=True,
        includedVAFs=True,
        plotIMDfigure=True,
        plotRainfall=True,
        windowSize=args.window_size,
        max_cpu=args.max_cpu,
        probability=True,
    )

    expected = input_root / "output" / "vcf_files_corrected" / f"{args.project}_clustered" / f"{args.project}_clusters_of_clusters_imd_edit.txt"
    print(f"Expected edited clusters file: {expected}", flush=True)
    if not expected.exists():
        matches = sorted(glob.glob(str(input_root / "output" / "**" / "*clusters_of_clusters*imd*.txt"), recursive=True))
        if matches:
            print("Found clustered output candidates:", flush=True)
            for m in matches:
                print(f"  {m}", flush=True)
        else:
            print("WARNING: No clusters_of_clusters IMD output was found yet. Check SigProfilerClusters logs.", flush=True)


if __name__ == "__main__":
    main()
PYTHON
  chmod +x "$CLUSTER_PY"

  cat > "$COUNT_R" <<'RSCRIPT'
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) stop("Missing value for ", flag)
  args[[hit + 1]]
}

project <- get_arg("--project", "GLASS5")
input_root <- normalizePath(get_arg("--input-root"), mustWork = TRUE)
outdir <- normalizePath(get_arg("--outdir"), mustWork = FALSE)
data_type <- get_arg("--data-type", "WGS")
clusters_file_arg <- get_arg("--clusters-file", "")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

find_clusters_file <- function(input_root, project) {
  candidates <- c(
    file.path(input_root, "output", "vcf_files_corrected", paste0(project, "_clustered"), paste0(project, "_clusters_of_clusters_imd_edit.txt")),
    file.path(input_root, "output", "vcf_files_corrected", paste0(project, "_clustered"), paste0(project, "_clusters_of_clusters_imd.txt")),
    file.path(input_root, "output", "vcf_files", paste0(project, "_clustered"), paste0(project, "_clusters_of_clusters_imd_edit.txt")),
    file.path(input_root, "output", "vcf_files", paste0(project, "_clustered"), paste0(project, "_clusters_of_clusters_imd.txt"))
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0) return(existing[[1]])

  matches <- list.files(
    file.path(input_root, "output"),
    pattern = "clusters_of_clusters.*imd.*\\.txt$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(matches) > 0) return(matches[[1]])
  NA_character_
}

clusters_file <- if (nzchar(clusters_file_arg)) normalizePath(clusters_file_arg, mustWork = TRUE) else find_clusters_file(input_root, project)
if (is.na(clusters_file) || !file.exists(clusters_file)) {
  stop("Could not find clusters_of_clusters IMD file under: ", input_root)
}

message("Reading clustered mutation file: ", clusters_file)
clusters_df <- read.delim(clusters_file, check.names = FALSE, stringsAsFactors = FALSE)

required <- c("samples", "clust_group")
missing <- setdiff(required, names(clusters_df))
if (length(missing) > 0) stop("Cluster file missing required columns: ", paste(missing, collapse = ", "))

cluster_df <- clusters_df[!is.na(clusters_df$clust_group) & clusters_df$clust_group != "", , drop = FALSE]

input_dir <- file.path(input_root, "input")
all_sample_files <- list.files(input_dir, pattern = "\\.txt$", full.names = FALSE)
all_samples <- sub("\\.txt$", "", all_sample_files)

if (length(all_samples) == 0) stop("No input sample files found in: ", input_dir)

kataegis_count <- cluster_df %>%
  distinct(.data$samples, .data$clust_group) %>%
  count(.data$samples, name = "kataegis_event_count")

kataegis_count_all <- data.frame(samples = all_samples) %>%
  left_join(kataegis_count, by = "samples") %>%
  mutate(kataegis_event_count = replace_na(.data$kataegis_event_count, 0L)) %>%
  rename(aliquot_barcode = .data$samples) %>%
  arrange(.data$aliquot_barcode)

# Add source columns to the details table so it can be safely combined later if needed.
cluster_df$data_type <- data_type
cluster_df$source_clusters_file <- clusters_file

lower_dtype <- tolower(data_type)
count_csv <- file.path(outdir, paste0("kataegis_event_count_all_", lower_dtype, ".csv"))
count_tsv <- file.path(outdir, paste0("kataegis_event_count_all_", lower_dtype, ".tsv"))
details_csv <- file.path(outdir, paste0("kataegis_details_all_", lower_dtype, ".csv"))
details_tsv <- file.path(outdir, paste0("kataegis_details_all_", lower_dtype, ".tsv"))
workspace <- file.path(outdir, paste0("kataegis_downstream_", lower_dtype, ".RData"))

write.csv(kataegis_count_all, file = count_csv, row.names = FALSE)
write.table(kataegis_count_all, file = count_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(cluster_df, file = details_csv, row.names = FALSE)
write.table(cluster_df, file = details_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
save(clusters_df, cluster_df, kataegis_count_all, file = workspace)

message("Saved count CSV: ", count_csv)
message("Saved count TSV: ", count_tsv)
message("Saved details CSV: ", details_csv)
message("Saved details TSV: ", details_tsv)
message("Saved RData: ", workspace)
RSCRIPT
}

run_filter() {
  local marker="filter"
  if is_done "$marker"; then
    log "Skipping Step 1 filter; marker exists: $(marker_path "$marker")"
    return 0
  fi
  [[ -n "$SQL" ]] || fail "--sql is required for Step 1 filter."
  [[ -s "$SQL" ]] || fail "SQL file not found or empty: $SQL"

  log "Step 1: repeat-region filtering"
  log "Output RData: $FILTER_RDATA"
  if [[ "$DRY_RUN" == "1" ]]; then
    log "DRY RUN: Rscript $SCRIPT_DIR/01_filter_repeat_regions.R ..."
    return 0
  fi

  activate_env "$R_ENV"
  Rscript "$SCRIPT_DIR/01_filter_repeat_regions.R" \
    --workdir "$FILTER_DIR" \
    --sql "$SQL" \
    --dsn "$DSN" \
    --out-rdata "$FILTER_RDATA" \
    --stats "$FILTER_STATS" \
    2>&1 | tee "$LOG_DIR/01_filter_repeat_regions.log"

  [[ -s "$FILTER_RDATA" ]] || fail "Step 1 did not create expected RData: $FILTER_RDATA"
  mark_done "$marker"
}

run_input() {
  local marker="input"
  if is_done "$marker"; then
    log "Skipping Step 2 input; marker exists: $(marker_path "$marker")"
    return 0
  fi
  [[ -s "$FILTER_RDATA" ]] || fail "Missing Step 1 RData: $FILTER_RDATA. Run --steps filter first."

  log "Step 2: create kataegis/SigProfiler-style input files"
  log "Input root: $INPUT_ROOT"
  if [[ "$DRY_RUN" == "1" ]]; then
    log "DRY RUN: Rscript $SCRIPT_DIR/02_make_kataegis_input.R ..."
    return 0
  fi

  activate_env "$R_ENV"
  Rscript "$SCRIPT_DIR/02_make_kataegis_input.R" \
    --filter-rdata "$FILTER_RDATA" \
    --outroot "$INPUT_ROOT" \
    --dsn "$DSN" \
    --subtype "$SUBTYPE" \
    --project-name "$PROJECT_NAME" \
    --genome "$GENOME" \
    --summary "$INPUT_SUMMARY" \
    2>&1 | tee "$LOG_DIR/02_make_kataegis_input.log"

  [[ -s "$INPUT_SUMMARY" ]] || fail "Step 2 did not create expected summary: $INPUT_SUMMARY"
  mark_done "$marker"
}

run_cluster_for_type() {
  local dtype_lower="$1"
  local dtype_upper
  dtype_upper=$(uppercase_dtype "$dtype_lower")
  local input_root="$INPUT_ROOT/$dtype_upper"
  local marker="cluster_${dtype_lower}"

  if is_done "$marker"; then
    log "Skipping Step 3 clustered mutation for $dtype_upper; marker exists: $(marker_path "$marker")"
    return 0
  fi

  [[ -d "$input_root/input" ]] || fail "Missing input directory: $input_root/input"
  local n_files
  n_files=$(find "$input_root/input" -maxdepth 1 -type f -name '*.txt' | wc -l | awk '{print $1}')
  if [[ "$n_files" == "0" ]]; then
    fail "No .txt input files found for $dtype_upper in $input_root/input"
  fi

  log "Step 3: SigProfilerSimulator + SigProfilerClusters for $dtype_upper ($n_files samples)"
  if [[ "$DRY_RUN" == "1" ]]; then
    log "DRY RUN: python $SCRIPT_DIR/03_run_sigprofiler_clusters.py --input-root $input_root ..."
    return 0
  fi

  activate_env "$CLUSTER_ENV"
  python "$SCRIPT_DIR/03_run_sigprofiler_clusters.py" \
    --project "$PROJECT_NAME" \
    --input-root "$input_root" \
    --genome "$GENOME" \
    --sim-contexts "$SIM_CONTEXTS" \
    --simulations "$SIMULATIONS" \
    --cluster-context "$CLUSTER_CONTEXT" \
    --cluster-subcontexts "$CLUSTER_SUBCONTEXTS" \
    --window-size "$WINDOW_SIZE" \
    --max-cpu "$MAX_CPU" \
    2>&1 | tee "$LOG_DIR/03_sigprofiler_clusters_${dtype_lower}.log"

  mark_done "$marker"
}

run_count_for_type() {
  local dtype_lower="$1"
  local dtype_upper
  dtype_upper=$(uppercase_dtype "$dtype_lower")
  local input_root="$INPUT_ROOT/$dtype_upper"
  local outdir="$DOWNSTREAM_ROOT/$dtype_upper"
  local marker="count_${dtype_lower}"

  if is_done "$marker"; then
    log "Skipping Step 4 count table for $dtype_upper; marker exists: $(marker_path "$marker")"
    return 0
  fi

  [[ -d "$input_root/input" ]] || fail "Missing input directory: $input_root/input"
  mkdir -p "$outdir"

  log "Step 4: kataegis event count table for $dtype_upper"
  log "Output directory: $outdir"
  if [[ "$DRY_RUN" == "1" ]]; then
    log "DRY RUN: Rscript $SCRIPT_DIR/04_make_kataegis_count_table.R --input-root $input_root ..."
    return 0
  fi

  activate_env "$R_ENV"
  Rscript "$SCRIPT_DIR/04_make_kataegis_count_table.R" \
    --project "$PROJECT_NAME" \
    --input-root "$input_root" \
    --outdir "$outdir" \
    --data-type "$dtype_upper" \
    2>&1 | tee "$LOG_DIR/04_kataegis_count_${dtype_lower}.log"

  local expected="$outdir/kataegis_event_count_all_${dtype_lower}.csv"
  [[ -s "$expected" ]] || fail "Step 4 did not create expected file: $expected"
  mark_done "$marker"
}

write_helpers

if has_step filter || has_step input; then
  [[ -n "$SQL" ]] || fail "--sql is required for filter/input steps."
  [[ -s "$SQL" ]] || fail "SQL file not found or empty: $SQL"
fi

log "Kataegis pipeline workdir: $WORKDIR"
log "Requested steps: $STEPS"
log "Data types: $DATA_TYPES"
log "Project: $PROJECT_NAME"
log "Genome: $GENOME"
log "Input root: $INPUT_ROOT"
log "Downstream root: $DOWNSTREAM_ROOT"

if has_step filter; then
  run_filter
fi

if has_step input; then
  run_input
fi

IFS=',' read -r -a DATA_TYPE_ARRAY <<< "$DATA_TYPES"
for raw_dtype in "${DATA_TYPE_ARRAY[@]}"; do
  dtype_lower=$(lowercase_dtype "$(echo "$raw_dtype" | xargs)")
  if has_step cluster; then
    run_cluster_for_type "$dtype_lower"
  fi
  if has_step count; then
    run_count_for_type "$dtype_lower"
  fi
done

log "Kataegis pipeline finished."
log "Main count outputs are under: $DOWNSTREAM_ROOT/<DATA_TYPE>/"
