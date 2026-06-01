# GLASS-I: ShatterSeek
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(ShatterSeek)
library(tidyverse)
library(GenomicRanges)
library(parallel)
library(VariantAnnotation)
library(RPostgres)

#### Define Paths ####
sv_folder <- "~/svaba"                                                          ## data created with GLASS-I_SvABA.sbatch
outdir <- "~/shatterseek"

valid_chroms <- c(as.character(1:22), "X")

#### Import SCNA from GATK segment files ####
con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")               ## data available on Synapse (Synapse ID: syn17038081)

surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") 

gold.set <- dbGetQuery(con, "SELECT * FROM analysis.gold_set")%>%
  semi_join(surgeries, by = "case_barcode")%>%
  pivot_longer(
    cols = c(tumor_barcode_a, tumor_barcode_b),
    names_to = NULL,
    values_to = "aliquot_barcode"
  )

scna.data <- dbGetQuery(con, "SELECT * FROM variants.gatk_seg")%>%
  semi_join(gold.set, by = "aliquot_barcode")%>%
  filter(str_detect(aliquot_barcode, "-WGS-"))%>%
  mutate(
    pos = str_remove_all(pos, "[\\[\\) ]"),
    start = as.numeric(str_split_fixed(pos, ",", 2)[, 1]),
    end   = (as.numeric(str_split_fixed(pos, ",", 2)[, 2]) - 1)
  ) 

dbDisconnect(con)

# Transform log2 values to total CN numbers
log2_to_total_cn <- function(log2cr,  max_cn = 20) {
  cn <- round(2 * (2 ^ log2cr))
  cn[cn > max_cn] <- max_cn
  return(as.integer(cn))
}

merged_scna <- scna.data %>%
  mutate(
    chromosome = as.character(chrom),
    total_cn = log2_to_total_cn(log2_copy_ratio),
  ) %>%
  filter(chromosome %in% valid_chroms)%>%
  arrange(aliquot_barcode, chromosome, start) %>%
  group_by(aliquot_barcode, chromosome) %>%
  mutate(
    new_block = 
      is.na(lag(end)) |
      start != lag(end) + 1|
      total_cn != lag(total_cn)
  ) %>%
  mutate(block_id = cumsum(new_block)) %>%
  group_by(aliquot_barcode, chromosome, block_id, total_cn) %>%
  summarise(
    start = min(start),
    end   = max(end),
    .groups = "drop"
  )%>%
  dplyr::select(aliquot_barcode, chromosome, start, end, total_cn)

#### Import SVs from SvABA somatic VCFs ####
vcf_files <- list.files(sv_folder, pattern = "svaba\\.somatic\\.sv\\.vcf\\.gz$", full.names = T)

parse_bnd_alt <- function(alt) {
  alt <- as.character(alt)
  
  if (is.na(alt) || nchar(alt) == 0) return(rep(NA_character_, 3))
  
  # Extract mate chrom:pos
  coord <- str_extract(alt, "([0-9XY]+|chr[0-9XY]+):[0-9]+")
  if (is.na(coord)) return(rep(NA_character_, 3))
  
  chrom2 <- sub(":.*", "", coord)
  pos2 <- as.integer(sub(".*:", "", coord))
  
  # Extract bracket orientation
  brackets <- unlist(str_extract_all(alt, "\\[|\\]"))
  strand1 <- if (length(brackets) >= 1) ifelse(brackets[1] == "[", "-", "+") else NA_character_
  
  c(chrom2, as.character(pos2), strand1)
}

parse_svaba_vcf <- function(vcf_file) {
  vcf <- readVcf(vcf_file, genome = "hg19")
  rr <- rowRanges(vcf)
  alt_vec <- as.character(unlist(alt(vcf)))
  
  df <- tibble(
    aliquot_barcode = sub("\\.svaba\\.somatic\\.sv\\.vcf\\.gz$", "", basename(vcf_file)),
    chrom1 = as.character(seqnames(rr)),
    pos1 = start(rr),
    alt = alt_vec
  )
  
  # Extract mate info from ALT
  bnd_info <- t(vapply(alt_vec, parse_bnd_alt, character(3)))
  colnames(bnd_info) <- c("chrom2", "pos2", "strand1")
  df <- bind_cols(df, as.data.frame(bnd_info, stringsAsFactors = FALSE))
  df$pos2 <- as.integer(df$pos2)
  
  # Keep only valid chromosomes
  df <- df %>%
    mutate(
      chrom1 = sub("^chr", "", chrom1),
      chrom2 = sub("^chr", "", chrom2)
    ) %>%
    filter(chrom1 %in% valid_chroms, chrom2 %in% valid_chroms)
  
  # =========================
  # Pair mates
  # =========================
  # Assuming each SV has exactly two BNDs. Match by chrom2:pos2 -> chrom1:pos1
  mates <- df %>%
    dplyr::select(chrom1, pos1, strand1) %>%
    mutate(chrom2 = chrom1, 
           pos2 = pos1, 
           strand2 = strand1)%>%
    dplyr::select(chrom2, pos2, strand2)
  
  df <- df %>%
    left_join(mates, by = c("chrom2", "pos2"))
  
  # =========================
  # Classify SVs
  # =========================
  df <- df %>%
    mutate(
      svclass = case_when(
        chrom1 != chrom2 ~ "TRA",
        strand1 == "-" & strand2 == "+" ~ "DUP",
        strand1 == "+" & strand2 == "-" ~ "DEL",
        (strand1 == "+" & strand2 == "+" & pos2 > pos1) |
          (strand1 == "-" & strand2 == "-" & pos2 < pos1) ~ "h2hINV",
        (strand1 == "+" & strand2 == "+" & pos2 < pos1) |
          (strand1 == "-" & strand2 == "-" & pos2 > pos1) ~ "t2tINV",
        TRUE ~ NA_character_
      )
    )
  
  df
}

merged_sv <- bind_rows(lapply(vcf_files, parse_svaba_vcf))%>%
  unique(.) %>%
  filter(
    chrom1 %in% valid_chroms,
    chrom2 %in% valid_chroms,
    !str_detect(alt, "hs37d5|GL000|KI270|decoy|random|_alt")
  )%>%
  rowwise() %>%
  mutate(
    # Create an "unordered" SV key so both mates map to the same key
    sv_key = paste0(
      paste0(sort(c(chrom1, chrom2)), collapse = "_"),
      ":",
      paste0(sort(c(pos1, pos2)), collapse = "_")
    )
  ) %>%
  ungroup() %>%
  group_by(aliquot_barcode, sv_key) %>%
  slice_min(pos1) %>%   # pick one representative row
  ungroup() %>%
  dplyr::select(-sv_key)

#### Keep only SCNA aliquots with matching SVs ####
scna_aliquots <- unique(merged_scna$aliquot_barcode)
sv_aliquots <- unique(merged_sv$aliquot_barcode)
aliquots_to_process <- intersect(scna_aliquots, sv_aliquots)

saveRDS(aliquots_to_process, file = file.path(outdir, "shatterseek_aliquots_to_process.rds"))

##### ShatterSeek Function #####
run_shatterseek <- function(alq_id) {
  
  log <- character()
  t0 <- Sys.time()
  
  add_log <- function(msg) {
    log <<- c(log, paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg))
  }
  
  add_log(paste0("START aliquot: ", alq_id))
  
  ## ---- SVs ----
  ind_sv <- merged_sv %>% filter(aliquot_barcode == alq_id)
  
  ## ---- SCNAs ----
  ind_scna <- merged_scna %>%
    filter(aliquot_barcode == alq_id) %>%
    dplyr::select(chromosome, start, end, total_cn)
  
  dd <- ind_scna
  dd$total_cn[dd$total_cn == 0] <- 150000
  dd$total_cn[is.na(dd$total_cn)] <- 0
  
  dd <- as(dd,"GRanges")
  cov <- coverage(dd,weight = dd$total_cn)
  dd1 <- as(cov,"GRanges")
  dd1 <- as.data.frame(dd1)
  dd1 <- dd1[dd1$score !=0,]
  dd1 <- dd1[,c(1,2,3,6)]
  names(dd1) <- names(ind_scna)[1:4]
  dd1$total_cn[dd1$total_cn == 150000] <- 0
  ind_scna <- as_tibble(dd1)
  
  ## ---- Build objects ----
  res <- tryCatch({
    
    SV_data <- SVs(
      chrom1  = ind_sv$chrom1,
      pos1    = ind_sv$pos1,
      chrom2  = ind_sv$chrom2,
      pos2    = ind_sv$pos2,
      SVtype  = ind_sv$svclass,
      strand1 = ind_sv$strand1,
      strand2 = ind_sv$strand2
    )
    
    CN_data <- CNVsegs(
      chrom     = as.character(ind_scna$chromosome),
      start     = ind_scna$start,
      end       = ind_scna$end,
      total_cn  = ind_scna$total_cn
    )
    
    ct <- shatterseek(
      SV.sample  = SV_data,
      seg.sample = CN_data
    )
    
    add_log("ShatterSeek completed")
    
    ct
    
  }, error = function(e) {
    add_log(paste0("FAILED: ", conditionMessage(e)))
    return(NULL)
  })
  
  if (is.null(res)) {
    return(log)
  }
  
  ## ---- Save output ----
  ctpr_df <- as_tibble(res@chromSummary) %>%
    dplyr::mutate(aliquot_barcode = alq_id) %>%
    dplyr::select(aliquot_barcode, dplyr::everything())
  
  write_tsv(ctpr_df, file.path(outdir, paste0("shatterseek_", alq_id, ".tsv")))
  save(res, ctpr_df,
       file = file.path(outdir, paste0("shatterseek_", alq_id, ".RData")))
  
  return(log)
}

#### Run ShatterSeek in parallel ####
mclapply(aliquots_to_process, run_shatterseek, mc.preschedule = FALSE, mc.cores = 10)

message("Finished ShatterSeek run for GLASS cohort.")