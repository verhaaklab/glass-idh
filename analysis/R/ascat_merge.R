#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
})

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 8) {
  stop("Usage: ascat_merge.R <tumor_baf_csv> <tumor_logr_csv> ",
       "<normal_baf_csv> <normal_logr_csv> ",
       "<out_tumor_baf> <out_tumor_logr> ",
       "<out_normal_baf> <out_normal_logr>\n")
}

tbaf_files  <- strsplit(args[1], ",")[[1]]
tlogr_files <- strsplit(args[2], ",")[[1]]
nbaf_files  <- strsplit(args[3], ",")[[1]]
nlogr_files <- strsplit(args[4], ",")[[1]]
out_tbaf    <- args[5]
out_tlogr   <- args[6]
out_nbaf    <- args[7]
out_nlogr   <- args[8]

# full-outer-join helper
full_join_all <- function(dfs, by_cols) {
  reduce(dfs, full_join, by = by_cols)
}

read_baf_tbl <- function(path) {
    df <- read.table(path,
                     header = TRUE, sep = "\t",
                     row.names = 1, check.names = FALSE,
                     stringsAsFactors = FALSE) %>%
          rownames_to_column("snp")

    metric_col <- names(df)[ncol(df)]     # <-- this is the aliquot_barcode
    df %>%
      select(snp, chr, pos, all_of(metric_col))   # keep it verbatim
}

read_logr_tbl <- function(path) {
    df <- read.table(path,
                     header = TRUE, sep = "\t",
                     row.names = 1, check.names = FALSE,
                     stringsAsFactors = FALSE) %>%
          rownames_to_column("snp")

    metric_col <- names(df)[ncol(df)]
    df %>%
      select(snp, chr, pos, all_of(metric_col))
}

# Merge & write function
merge_and_write <- function(files, reader, out_path) {
  merged <- map(files, reader) %>%
    full_join_all(by = c("snp","chr","pos")) %>%
    select(snp, chr, pos, everything()) %>%
    column_to_rownames("snp")

  write.table(merged, file      = out_path,
              sep       = "\t",
              quote     = FALSE,
              row.names = TRUE,
              col.names = NA)
}

# 1) Tumor BAF
merge_and_write(tbaf_files,  read_baf_tbl,  out_tbaf)

# 2) Tumor LogR
merge_and_write(tlogr_files, read_logr_tbl, out_tlogr)

# 3) Normal BAF
merge_and_write(nbaf_files,  read_baf_tbl,  out_nbaf)

# 4) Normal LogR
merge_and_write(nlogr_files, read_logr_tbl, out_nlogr)

cat("Merge complete: all files with snp, chr, pos + sample columns.\n")







# #!/usr/bin/env Rscript

# suppressPackageStartupMessages({
#   library(dplyr)
# })

# # ------------------------------------------------------------------------
# # 1) Read command-line arguments
# # ------------------------------------------------------------------------
# args <- commandArgs(trailingOnly = TRUE)

# if (length(args) < 8) {
#   stop("Usage: ascat_merge.R <tumor_baf_files_csv> <tumor_logr_files_csv> <normal_baf_file> <normal_logr_file> <out_tumor_baf> <out_tumor_logr> <out_normal_baf> <out_normal_logr>\n")
# }

# # 1) Tumor files (comma-separated lists)
# tumor_baf_files  <- unlist(strsplit(args[1], split = ","))
# tumor_logr_files <- unlist(strsplit(args[2], split = ","))

# # 2) Normal is a single file each
# normal_baf_file  <- args[3]
# normal_logr_file <- args[4]

# # 3) Output paths
# out_tumor_baf    <- args[5]
# out_tumor_logr   <- args[6]
# out_normal_baf   <- args[7]
# out_normal_logr  <- args[8]

# # ------------------------------------------------------------------------
# # Print summary
# # ------------------------------------------------------------------------
# cat("Merging tumor BAF:\n", tumor_baf_files, "\n")
# cat("Merging tumor LogR:\n", tumor_logr_files, "\n")
# cat("Using single normal BAF:", normal_baf_file, "\n")
# cat("Using single normal LogR:", normal_logr_file, "\n")
# cat("Output (tumor BAF):",   out_tumor_baf,   "\n")
# cat("Output (tumor LogR):",  out_tumor_logr,  "\n")
# cat("Output (normal BAF):",  out_normal_baf,  "\n")
# cat("Output (normal LogR):", out_normal_logr, "\n\n")

# # # ------------------------------------------------------------------------
# # # 4) Define helper functions
# # # ------------------------------------------------------------------------
# # # A helper to load one BAF/LogR file, using the first column as row names.
# # # The resulting dataframe has columns: "chromosome", "position", <last_col>.
# # # We rename that last column to `sample_name`.
# # load_file <- function(filepath, sample_name) {
# #   df <- read.table(
# #     filepath,
# #     header = TRUE,
# #     sep = "\t",
# #     row.names = 1,             # first column (e.g. "1_561156") becomes row names
# #     stringsAsFactors = FALSE
# #   )
# #   # Now 'df' has columns: "chromosome", "position", and one value column (index 3).
# #   df <- df[, c("chromosome", "position", colnames(df)[3])]  
# #   colnames(df)[3] <- sample_name
# #   return(df)
# # }

# # # Merges a list of files by (chromosome, position).
# # # Each file is loaded with load_file() (which assigns a sample column name),
# # # and all are combined via base::merge(..., by=c("chromosome","position")).
# # merge_files <- function(files) {
# #   Reduce(
# #     function(x, y) merge(x, y, by = c("chromosome", "position"), all = FALSE, sort = FALSE),
# #     lapply(seq_along(files), function(i) {
# #       # For the sample name, we remove known suffixes or just use the file's basename
# #       f      <- files[i]
# #       s_name <- sub("\\.BAF\\.tumor\\.txt$|\\.LogR\\.tumor\\.txt$|\\.BAF\\.normal\\.txt$|\\.LogR\\.normal\\.txt$", "", basename(f))
# #       load_file(f, s_name)
# #     })
# #   )
# # }

# # # ------------------------------------------------------------------------
# # # 5) Merge TUMOR BAF/LogR
# # # ------------------------------------------------------------------------
# # if (length(tumor_baf_files) > 1) {
# #   tbaf_merged <- merge_files(tumor_baf_files)
# # } else {
# #   f <- tumor_baf_files[1]
# #   s <- sub("\\.BAF\\.tumor\\.txt$|\\.LogR\\.tumor\\.txt$|\\.BAF\\.normal\\.txt$|\\.LogR\\.normal\\.txt$", "", basename(f))
# #   tbaf_merged <- load_file(f, s)
# # }
# # rownames(tbaf_merged) <- paste(tbaf_merged$chromosome, tbaf_merged$position, sep = "_")
# # write.table(tbaf_merged, file = out_tumor_baf, sep = "\t", quote = FALSE, row.names = TRUE)

# # if (length(tumor_logr_files) > 1) {
# #   tlogr_merged <- merge_files(tumor_logr_files)
# # } else {
# #   f <- tumor_logr_files[1]
# #   s <- sub("\\.BAF\\.tumor\\.txt$|\\.LogR\\.tumor\\.txt$|\\.BAF\\.normal\\.txt$|\\.LogR\\.normal\\.txt$", "", basename(f))
# #   tlogr_merged <- load_file(f, s)
# # }
# # rownames(tlogr_merged) <- paste(tlogr_merged$chromosome, tlogr_merged$position, sep = "_")
# # write.table(tlogr_merged, file = out_tumor_logr, sep = "\t", quote = FALSE, row.names = TRUE)

# # # ------------------------------------------------------------------------
# # # 6) Load SINGLE normal BAF/LogR
# # # ------------------------------------------------------------------------
# # # We do not merge multiple normals, but rename columns similarly
# # nbaf_df <- {
# #   s <- sub("\\.BAF\\.tumor\\.txt$|\\.LogR\\.tumor\\.txt$|\\.BAF\\.normal\\.txt$|\\.LogR\\.normal\\.txt$", "", basename(normal_baf_file))
# #   load_file(normal_baf_file, s)
# # }
# # rownames(nbaf_df) <- paste(nbaf_df$chromosome, nbaf_df$position, sep = "_")
# # write.table(nbaf_df, file = out_normal_baf, sep = "\t", quote = FALSE, row.names = TRUE)

# # nlogr_df <- {
# #   s <- sub("\\.BAF\\.tumor\\.txt$|\\.LogR\\.tumor\\.txt$|\\.BAF\\.normal\\.txt$|\\.LogR\\.normal\\.txt$", "", basename(normal_logr_file))
# #   load_file(normal_logr_file, s)
# # }
# # rownames(nlogr_df) <- paste(nlogr_df$chromosome, nlogr_df$position, sep = "_")
# # write.table(nlogr_df, file = out_normal_logr, sep = "\t", quote = FALSE, row.names = TRUE)

# # cat("\nMerging complete.\n")


# # ------------------------------------------------------------------------
# # 4) Define helper functions
# # ------------------------------------------------------------------------
# # A helper to load one BAF/LogR file, using the first column as row names.
# # The resulting dataframe has columns: "chromosome", "position", <last_col>.
# # We rename that last column to `sample_name`.
# load_file <- function(filepath, sample_name) {
#   df <- read.table(
#     filepath,
#     header = TRUE,
#     sep = "\t",
#     row.names = 1,
#     stringsAsFactors = FALSE
#   )
#   df <- df[, c("chromosome", "position", colnames(df)[3])]
#   colnames(df)[3] <- sample_name
#   return(df)
# }

# # Merge a list of files by (chromosome, position) with full join,
# # then filter out rows with more than max_missing missing sample values.
# merge_and_filter <- function(files, max_missing = 1) {
#   # 1) load all files with sample-named columns
#   dfs <- lapply(seq_along(files), function(i) {
#     f      <- files[i]
#     s_name <- sub("\\.BAF\\.tumor\\.txt$|\\.LogR\\.tumor\\.txt$|\\.BAF\\.normal\\.txt$|\\.LogR\\.normal\\.txt$", "", basename(f))
#     load_file(f, s_name)
#   })

#   # 2) full-join them all
#   merged <- Reduce(
#     function(x, y) merge(x, y, by = c("chromosome", "position"), all = TRUE, sort = FALSE),
#     dfs
#   )

#   # 3) count NAs per locus across sample columns
#   sample_cols <- setdiff(colnames(merged), c("chromosome", "position"))
#   na_counts   <- rowSums(is.na(merged[, sample_cols]))

#   total    <- nrow(merged)
#   kept     <- sum(na_counts <= max_missing)
#   dropped  <- total - kept
#   pct_kept <- 100 * kept / total

#   cat(sprintf("[QC] merged %d loci; dropping %d loci (>%d NAs); keeping %d (%.1f%%)\n",
#               total, dropped, max_missing, kept, pct_kept))

#   # 4) filter and return
#   merged[na_counts <= max_missing, ]
# }

# # ------------------------------------------------------------------------
# # 5) Merge TUMOR BAF
# # ------------------------------------------------------------------------
# if (length(tumor_baf_files) > 1) {
#   # allow up to 1 missing sample per locus
#   tbaf_merged <- merge_and_filter(tumor_baf_files, max_missing = 1)
# } else {
#   f <- tumor_baf_files[1]
#   s <- sub("\\.BAF\\.tumor\\.txt$|\\.LogR\\.tumor\\.txt$|\\.BAF\\.normal\\.txt$|\\.LogR\\.normal\\.txt$", "", basename(f))
#   tbaf_merged <- load_file(f, s)
# }
# rownames(tbaf_merged) <- paste(tbaf_merged$chromosome, tbaf_merged$position, sep = "_")
# write.table(tbaf_merged, file = out_tumor_baf, sep = "\t", quote = FALSE, row.names = TRUE)

# # ------------------------------------------------------------------------
# # 6) Merge TUMOR LogR
# # ------------------------------------------------------------------------
# if (length(tumor_logr_files) > 1) {
#   # apply the same full-join + NA-filter for LogR
#   tlogr_merged <- merge_and_filter(tumor_logr_files, max_missing = 1)
# } else {
#   f <- tumor_logr_files[1]
#   s <- sub("\\.BAF\\.tumor\\.txt$|\\.LogR\\.tumor\\.txt$|\\.BAF\\.normal\\.txt$|\\.LogR\\.normal\\.txt$", "", basename(f))
#   tlogr_merged <- load_file(f, s)
# }
# rownames(tlogr_merged) <- paste(tlogr_merged$chromosome, tlogr_merged$position, sep = "_")
# write.table(tlogr_merged, file = out_tumor_logr, sep = "\t", quote = FALSE, row.names = TRUE)

# # ------------------------------------------------------------------------
# # 7) Copy/rename normal BAF & LogR
# # ------------------------------------------------------------------------
# nbaf_df <- {
#   s <- sub("\\.BAF\\.tumor\\.txt$|\\.LogR\\.tumor\\.txt$|\\.BAF\\.normal\\.txt$|\\.LogR\\.normal\\.txt$", "", basename(normal_baf_file))
#   df <- load_file(normal_baf_file, s)
#   rownames(df) <- paste(df$chromosome, df$position, sep = "_")
#   df
# }
# write.table(nbaf_df, file = out_normal_baf, sep = "\t", quote = FALSE, row.names = TRUE)

# nlogr_df <- {
#   s <- sub("\\.BAF\\.tumor\\.txt$|\\.LogR\\.tumor\\.txt$|\\.BAF\\.normal\\.txt$|\\.LogR\\.normal\\.txt$", "", basename(normal_logr_file))
#   df <- load_file(normal_logr_file, s)
#   rownames(df) <- paste(df$chromosome, df$position, sep = "_")
#   df
# }
# write.table(nlogr_df, file = out_normal_logr, sep = "\t", quote = FALSE, row.names = TRUE)

# cat("\nMerging complete.\n")
