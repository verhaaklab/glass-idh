# GLASS-I: Chromothripsis Annotation Table
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(tidyverse)
library(RPostgres)

# load data
con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")                    ## data available on Synapse (Synapse ID: syn17038081)

surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries")

dbDisconnect(con)

rdata_files <- list.files("~/shatterseek", pattern = "\\.RData$", full.names = TRUE) ## data created with GLASS-I_ShatterSeek.R

process_ctpr <- function(file_path) {
  load(file_path)
  
  sample_name <- substr(basename(file_path), 13, 42)
  
  ctpr_df <- ctpr_df %>%
    mutate(
      sample = sample_name,
      intrachrom_SVs = rowSums(dplyr::select(., number_DEL, number_DUP, number_h2hINV, number_t2tINV)),
      interchrom_SVs = rowSums(dplyr::select(., inter_number_DEL, inter_number_DUP, inter_number_h2hINV, inter_number_t2tINV)),
      max_number_oscillating_CN_segments_2_states = if_else(is.na(max_number_oscillating_CN_segments_2_states), 0, max_number_oscillating_CN_segments_2_states),
      high.confidence.test1 = if_else(intrachrom_SVs > 5 &
                                        max_number_oscillating_CN_segments_2_states > 6 &
                                        (chr_breakpoint_enrichment < 0.05 | pval_exp_cluster < 0.05) &
                                        pval_fragment_joins > 0.05 ,
                                      TRUE, FALSE),    
      high.confidence.test2 = if_else(intrachrom_SVs > 2 &
                                        (interchrom_SVs > 3 | number_TRA > 3) &
                                        max_number_oscillating_CN_segments_2_states > 6 &
                                        pval_fragment_joins > 0.05,
                                      TRUE, FALSE),
      low.confidence.test = if_else(intrachrom_SVs > 5 &
                                      max_number_oscillating_CN_segments_2_states > 3 &
                                      max_number_oscillating_CN_segments_2_states < 7 &
                                      (chr_breakpoint_enrichment < 0.05 | pval_exp_cluster < 0.05) &
                                      pval_fragment_joins > 0.05,
                                    TRUE, FALSE)
    )
  
  return(ctpr_df)
}

# Apply to all files and combine into one big dataframe
all_ctpr_df <- map_dfr(rdata_files, process_ctpr)

sample.data <- all_ctpr_df %>%
  group_by(aliquot_barcode)%>%
  summarise(
    sample_barcode = substr(aliquot_barcode[1], 1, 15),
    chromothripsis.high.conf = if_else(all(high.confidence.test1==F & high.confidence.test2==F), 0,1),
    chromothripsis.low.conf = if_else(all(low.confidence.test==F), 0,1),
    chromothripsis.any.conf = if_else(all(high.confidence.test1==F & high.confidence.test2==F & low.confidence.test==F), 0,1))%>%
  ungroup()%>%
  inner_join(surgeries, by = "sample_barcode")%>%
  select(aliquot_barcode, sample_barcode, case_barcode, surgery_number, idh_codel_subtype, chromothripsis.high.conf, chromothripsis.low.conf, chromothripsis.any.conf)

write_csv(sample.data, "~/chromothripsis_shatterseek_samples.csv")