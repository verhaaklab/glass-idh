#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list = list(
  make_option(c("--seg_file"),  type="character",
              help="Case‐level phased‐segments TSV (sample_id = pair_barcode)"),
  make_option(c("--pair"),      type="character",
              help="The pair_barcode to extract for this aliquot"),
  make_option(c("--aliquot"),   type="character",
              help="The aliquot_barcode for which we’re preparing the file"),
  make_option(c("--out_txt"),   type="character",
              help="Output ASCAT‐style TSV for this aliquot")
)

opt = parse_args(OptionParser(option_list=option_list))

# 1) Read the case‐level segments (sample_id = pair_barcode)
segs = fread(opt$seg_file,
             header=TRUE, sep="\t",
             colClasses=list(character="sample_id", character="chrom"))

# 2) Keep only rows where sample_id == the requested pair_barcode
segs = segs[sample_id == opt$pair]

# 3) Compute nMajor = max(cn_a, cn_b), nMinor = min(cn_a, cn_b)
segs[, `:=`(
  nMajor = pmax(cn_a, cn_b),
  nMinor = pmin(cn_a, cn_b)
)]

# 4) Build ASCAT‐style table: sample, chromosome, start, end, nMajor, nMinor
out_dt = segs[, .(
  sample     = opt$aliquot,
  chromosome = chrom,
  startpos   = start,
  endpos     = end,
  nMajor     = nMajor,
  nMinor     = nMinor
)]

# 5) Write to --out_txt (creating directory if necessary)
dir.create(dirname(opt$out_txt), recursive=TRUE, showWarnings=FALSE)
fwrite(out_dt,
       opt$out_txt,
       sep="\t",
       quote=FALSE,
       col.names=TRUE)

cat(sprintf("✓ Wrote ASCAT‐style segments to: %s\n", opt$out_txt))