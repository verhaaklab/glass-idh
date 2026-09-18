#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ASCAT)
    library(dplyr)
    library(tibble)
})

# Expected command-line arguments:
#  1) mode           => "single" or "multi"
#  2) tumor_baf      => single-file or merged file
#  3) tumor_logr     => single-file or merged file
#  4) normal_baf     => single-file or merged file
#  5) normal_logr    => single-file or merged file
#  6) gender         => "XX" or "XY" for single, or comma-separated for multi
#  7) gcfile         => path to GC content reference
#  8) rtfile         => path to replication timing reference
#  9) pdfdir         => output directory for plots
# 10) output_qc      => path to QC metrics output
# 11) output_seg     => path to segments output

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 11) {
    stop("Insufficient arguments. Exiting.\n")
}

mode         <- args[1]
tumor_baf    <- args[2]
tumor_logr   <- args[3]
normal_baf   <- args[4]
normal_logr  <- args[5]
gender_arg   <- args[6]   
gcfile       <- args[7]
rtfile       <- args[8]
pdfdir       <- args[9]
output_qc    <- args[10]
output_seg   <- args[11]

cat("Running ascat_analysis.R with mode =", mode, "\n")

# Parse and map gender_arg into ASCAT codes
raw_genders <- strsplit(gender_arg, ",")[[1]]
map_gender  <- function(g) {
  g0 <- tolower(g)
  if (g0 %in% c("female","f","xx")) return("XX")
  if (g0 %in% c("male",  "m","xy")) return("XY")
  warning(sprintf("Unrecognized gender '%s'; defaulting to 'XX'", g))
  return("XX")
}
mapped_genders <- vapply(raw_genders, map_gender, FUN.VALUE = character(1))

# For multi-mode, pass all; for single-mode, just the first
if (mode == "multi") {
  genders <- mapped_genders
  cat("Using per-sample genders:", paste(genders, collapse=","), "\n")
} else {
  genders <- mapped_genders[1]
  cat("Using single gender:", genders, "\n")
}

# 1) Load data
ascat.bc <- ascat.loadData(
    Tumor_BAF_file    = tumor_baf,
    Tumor_LogR_file   = tumor_logr,
    Germline_BAF_file = normal_baf,
    Germline_LogR_file= normal_logr,
    genomeVersion     = "hg19", 
    gender           = genders
)

# 2) Correct for GC content / replication timing
ascat.gcrt <- ascat.correctLogR(
    ascat.bc,
    GCcontentfile    = gcfile,
    replictimingfile = rtfile
)

# 3) Segmentation step differs for single vs multi
if (mode == "multi") {
    # multi-sample segmentation
    cat("Running ascat.asmultipcf(...)\n")
    ascat.seg <- ascat.asmultipcf(ascat.gcrt, out.dir = pdfdir, seed=42)
} else {
    # single-sample segmentation
    cat("Running ascat.aspcf(...)\n")
    ascat.seg <- ascat.aspcf(ascat.gcrt, out.dir = pdfdir, seed=42)
}

# 4) Run ASCAT
ascat.output <- ascat.runAscat(
    ascat.seg,
    gamma   = 1,
    pdfPlot = TRUE,
    img.dir = pdfdir
)

# 5) Extract QC & segments
QC <- ascat.metrics(ascat.seg, ascat.output) %>%
      tibble::rownames_to_column("aliquot_barcode")
segments <- ascat.output[["segments"]]

# 6) Write outputs
write.table(QC,      file=output_qc, sep="\t", row.names=FALSE, quote=FALSE)
write.table(segments,file=output_seg,sep="\t", row.names=FALSE, quote=FALSE)

cat("ASCAT analysis completed. Mode =", mode, "\n")


# ### Previous::
# #!/usr/bin/env Rscript

# suppressPackageStartupMessages({
#     library(optparse)
#     library(ASCAT)
# })

# # option_list <- list(
# #     make_option(c("--tumour_bam"), type = "character", help = "Path to the tumour BAM file"),
# #     make_option(c("--normal_bam"), type = "character", help = "Path to the normal BAM file"),
# #     make_option(c("--output_prefix"), type = "character", help = "Prefix for output files"),
# #     make_option(c("--allelecounter_exe"), type = "character", help = "Path to the alleleCounter executable"),
# #     make_option(c("--alleles_prefix"), type = "character", help = "Prefix for alleles file"),
# #     make_option(c("--loci_prefix"), type = "character", help = "Prefix for loci file"),
# #     make_option(c("--gender"), type = "character", help = "Gender of the sample"),
# #     make_option(c("--genomeVersion"), type = "character", default = "hg19", help = "Genome version (hg19/hg38)"),
# #     make_option(c("--nthreads"), type = "integer", default = 1, help = "Number of threads"),
# #     make_option(c("--GCcontentfile"), type = "character", help = "File containing GC content. [Required]"),
# #     make_option(c("--replictimingfile"), type = "character", help = "File containing replication timing information. [Required]")

# # )

# # opt_parser <- OptionParser(option_list=option_list)
# # opt <- parse_args(opt_parser)

# # # Setting the correct working directory (BETA VERSION, potentially adjust!)
# # current_wd <- getwd()
# # output_path <- file.path(current_wd, opt$output_prefix)
# # wd <- dirname(output_path)
# # dir.create(wd, recursive=T, showWarnings=F)
# # setwd(wd)

# # Prepare allele frequencies files
# # ascat.prepareHTS(
# #   tumourseqfile = opt$tumour_bam,
# #   normalseqfile = opt$normal_bam,
# #   tumourname = gsub(".bam$", "", basename(opt$tumour_bam)),
# #   normalname = gsub(".bam$", "", basename(opt$normal_bam)),
# #   allelecounter_exe = opt$allelecounter_exe,
# #   alleles.prefix = opt$alleles_prefix,
# #   loci.prefix = opt$loci_prefix,
# #   gender = opt$gender,
# #   genomeVersion = opt$genomeVersion,
# #   nthreads = opt$nthreads,
# #   tumourLogR_file = paste0(opt$output_prefix, "_tumour_LogR.txt"),
# #   tumourBAF_file = paste0(opt$output_prefix, "_tumour_BAF.txt"),
# #   normalLogR_file = paste0(opt$output_prefix, "_normal_LogR.txt"),
# #   normalBAF_file = paste0(opt$output_prefix, "_normal_BAF.txt")
# # )

# # # Load the data
# # ascat.bc = ascat.loadData(
# #   Tumor_LogR_file = paste0(opt$output_prefix, "_tumour_LogR.txt"),
# #   Tumor_BAF_file = paste0(opt$output_prefix, "_tumour_BAF.txt"),
# #   Germline_LogR_file = paste0(opt$output_prefix, "_normal_LogR.txt"),
# #   Germline_BAF_file = paste0(opt$output_prefix, "_normal_BAF.txt"),
# #   genomeVersion = opt$genomeVersion,
# #   gender = opt$gender
# # )

# # # Correct for GC content
# # ascat.gcrt = ascat.correctLogR(ascat.bc, GCcontentfile = opt$GCcontentfile, replictimingfile = opt$replictimingfile)

# # # ASCAT Segmentation and Purity/Ploidy estimation
# # ascat.seg = ascat.aspcf(ascat.gcrt, seed = 42)

# # # Running ASCAT
# # ascat.output = ascat.runAscat(ascat.seg, gamma = 1, pdfPlot = TRUE, img.dir = opt$output_prefix, img.prefix = "ascat_output")

# # # Extracting QC metrics and CNVs
# # QC = ascat.metrics(ascat.seg, ascat.output)
# # write.table(QC, file=paste0(opt$output_prefix, "_QC_metrics.txt"), sep="\t", quote=FALSE, row.names=FALSE)

# # segments = ascat.output[["segments"]]
# # write.table(segments, file=paste0(opt$output_prefix, "_segments.txt"), sep="\t", quote=FALSE, row.names=FALSE)

# # cat("ASCAT analysis completed successfully\n")

# # Parse command line arguments
# args <- commandArgs(trailingOnly = TRUE)
# baft <- args[1]    # BAF tumor file
# bafn <- args[2]    # BAF normal file
# logrt <- args[3]   # LogR tumor file
# logrn <- args[4]   # LogR normal file
# gender <- args[5]  # Gender parameter
# outdir <- args[6]  # Output directory
# gversion <- args[7] # Genome version (hg19)
# gcfile <- args[8]  # GC content file
# rtfile <- args[9]  # Replication timing file
# pdfdir <- args[10] # PDF directory
# output_qc <- args[11] # Output QC metrics
# output_seg <- args[12] # Output Segmentation metrics 

# # Perform data loading
# ascat.bc <- ascat.loadData(
#     Tumor_LogR_file = logrt,
#     Tumor_BAF_file = baft,
#     Germline_LogR_file = logrn,
#     Germline_BAF_file = bafn,
#     genomeVersion = "hg19",
#     gender = gender
# )

# ascat.gcrt <- ascat.correctLogR(
#   ascat.bc, 
#   GCcontentfile = gcfile, 
#   replictimingfile = rtfile)

# ascat.seg <- ascat.aspcf(
#   ascat.gcrt, 
#   out.dir = pdfdir,
#   seed = 42
#   )

# ascat.output <- ascat.runAscat(
#   ascat.seg, 
#   gamma = 1, 
#   pdfPlot = TRUE, 
#   img.dir = pdfdir
#   )

# QC <- ascat.metrics(
#   ascat.seg, 
#   ascat.output
#   )

# segments <- ascat.output[["segments"]]

# # Write output
# write.table(QC, file = output_qc, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(segments, file = output_seg, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# cat("ASCAT analysis completed successfully\n")