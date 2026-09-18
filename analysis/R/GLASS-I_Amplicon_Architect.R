# GLASS-I: Amplicon Architect
# Author: E. Kocakavuk
# Date: 2026-06-18

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(tidyverse)
library(RPostgres)
library(pool)
library(stringr)
library(fs)
library(vroom)

### GLASS database connection
con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")

gold_set <- dbGetQuery(con, "SELECT * FROM analysis.gold_set")                     ## data available on Synapse (Synapse ID: syn17038081)

dbDisconnect(con)

## Read in AmpliconSuite results
aa_files <- dir_ls(
 c("~/ampsuite/run/",                                                             ## created with ampsuite.smk + ampsuite_manifest.py
  "~/ampsuite/run/"
  ),
  recurse = TRUE,
  glob    = "*/*/*_classification/*_result_table.tsv"
)

aa_query <- vroom::vroom(
  aa_files,
  delim = "\t",
  .name_repair = "minimal",
  progress = TRUE
)

aa_output <- aa_query %>%  
  filter(`Sample name` %in% c(gold_set$tumor_barcode_a, gold_set$tumor_barcode_b), 
          grepl("-WGS-", `Sample name`)) %>% 
  group_by(`Sample name`) %>% 
  mutate(ecDNA = ifelse(any(Classification == "ecDNA"), "ecDNA", "no_ecDNA")) %>% 
  mutate(ecDNA = ifelse(is.na(ecDNA), "no_ecDNA", ecDNA)) %>% #NA in this context means that there were no amplicons detected
  ungroup() %>% 
  select(aliquot_barcode = `Sample name`, ecDNA) %>% distinct()

aa_genes_files <- dir_ls(
 c("~/ampsuite/run/",
  "~/ampsuite/run/"
  ),
  recurse = TRUE,
  glob    = "*/*/*_classification/*_gene_list.tsv"
)

aa_genes_query <- vroom::vroom(
  aa_genes_files,
  delim = "\t",
  .name_repair = "minimal",
  progress = TRUE
)

aa_genes_output <- aa_genes_query %>%  
  filter(sample_name %in% c(gold_set$tumor_barcode_a, gold_set$tumor_barcode_b), 
          grepl("-WGS-", sample_name)) %>% 
  filter(grepl("ecDNA", feature))

#Export dataframe
write_tsv(aa_output, "~/glass5_ampsuite_20260317.tsv")
write_tsv(aa_genes_output, "~/glass5_ampsuite_genes_20260424.tsv")
