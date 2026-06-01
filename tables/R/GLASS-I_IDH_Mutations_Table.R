# GLASS-I: IDH Mutations Table
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(RPostgres)
library(tidyverse)

## Some housekeeping data
variant_list <- c("FRAME_SHIFT_DEL","FRAME_SHIFT_INS","IN_FRAME_DEL","IN_FRAME_INS","MISSENSE","NONSENSE","SPLICE_SITE")

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")             ## data available on Synapse (Synapse ID: syn17038081)
hm_info <- dbGetQuery(con, "SELECT * FROM analysis.tumor_clinical_comparison_2024")%>%
  mutate(HM_group = ifelse(sbs11_r_hypermutation==1,"HM","NHM"))%>%
  select(case_barcode, tumor_pair_barcode, HM_group)

df_gs <- dbGetQuery(con, "SELECT * FROM analysis.gold_set gs JOIN clinical.subtypes cs ON cs.case_barcode = gs.case_barcode WHERE cs.idh_codel_subtype NOT LIKE 'IDHwt'")%>%
  select(c(1,2,6))

mutation.types <- dbGetQuery(con, "select * from variants.passanno
where (gene_symbol = 'IDH1' and protein_change like 'p.R132%') OR (gene_symbol = 'IDH2' and protein_change like 'p.R172%')")[,c(1,17)]

## counting frequency table for IDH1
idh1_short <- dbGetQuery(con, "SELECT pg.* from variants.pgeno pg JOIN analysis.gold_set gs ON pg.tumor_pair_barcode = gs.tumor_pair_barcode WHERE pg.gene_symbol = 'IDH1'")%>%
  semi_join(df_gs, by = c("case_barcode","tumor_pair_barcode"))%>% 
  semi_join(mutation.types, by = c("variant_id"))%>% 
  mutate(driver_change = case_when(mutect2_call_a==1 & mutect2_call_b==0 ~ "P", mutect2_call_a==0 & mutect2_call_b==1 ~ "R", mutect2_call_a==1 & mutect2_call_b==1 ~ "S", TRUE ~ "none"))%>% 
  filter(!driver_change=="none" & variant_classification %in% variant_list)%>%
  full_join(df_gs, by = c("case_barcode","tumor_pair_barcode"))%>% 
  full_join(mutation.types, by = c("variant_id"))%>% 
  mutate(driver_status = ifelse(is.na(gene_symbol),"wildtype","mutant"))%>%
  left_join(hm_info, by = c("case_barcode","tumor_pair_barcode"))%>% 
  group_by(case_barcode)%>% 
  mutate(multihit = ifelse(n()>1, "yes","no"))%>% 
  slice_max(order_by = af_b, n = 1)%>%
  slice(1)%>% 
  ungroup()%>%
  mutate(gene_symbol = "IDH1")%>%
  select(case_barcode, tumor_pair_barcode, idh_codel_subtype, HM_group, gene_symbol, driver_status, driver_change, protein_change)%>%
  filter(!is.na(tumor_pair_barcode))

## counting frequency table for IDH2
idh2_short <- dbGetQuery(con, "SELECT pg.* from variants.pgeno pg JOIN analysis.gold_set gs ON pg.tumor_pair_barcode = gs.tumor_pair_barcode  WHERE pg.gene_symbol = 'IDH2'")%>%
  semi_join(df_gs, by = c("case_barcode","tumor_pair_barcode"))%>% 
  semi_join(mutation.types, by = c("variant_id"))%>% 
  mutate(driver_change = case_when(mutect2_call_a==1 & mutect2_call_b==0 ~ "P", mutect2_call_a==0 & mutect2_call_b==1 ~ "R", mutect2_call_a==1 & mutect2_call_b==1 ~ "S", TRUE ~ "none"))%>% 
  filter(!driver_change=="none" & variant_classification %in% variant_list)%>%
  full_join(df_gs, by = c("case_barcode","tumor_pair_barcode"))%>% 
  full_join(mutation.types, by = c("variant_id"))%>% 
  mutate(driver_status = ifelse(is.na(gene_symbol),"wildtype","mutant"))%>%
  left_join(hm_info, by = c("case_barcode","tumor_pair_barcode"))%>% 
  group_by(case_barcode)%>% 
  mutate(multihit = ifelse(n()>1, "yes","no"))%>% 
  slice_max(order_by = af_b, n = 1)%>%
  slice(1)%>% 
  ungroup()%>%
  mutate(gene_symbol = "IDH2")%>%
  select(case_barcode, tumor_pair_barcode, idh_codel_subtype, HM_group, gene_symbol, driver_status, driver_change, protein_change)%>%
  filter(!is.na(tumor_pair_barcode))

snv_data <- rbind(idh1_short, idh2_short)

## counting and combining frequency table for PI3K
idh1_idh2_short <- subset(snv_data, snv_data$gene_symbol=="IDH1"|snv_data$gene_symbol=="IDH2")%>%
  mutate(gene_symbol_old = gene_symbol,
         gene_symbol = "IDH1/2",
         driver_status_old = driver_status,
         driver_change_old = if_else(is.na(driver_change), "NA", driver_change),
         protein_change_old = if_else(is.na(protein_change), "NA", protein_change),)%>% 
  group_by(case_barcode) %>% 
  mutate(driver_status = if_else(any(driver_status_old == "mutant"), "mutant", "wildtype"),
         driver_change = if_else(all(driver_change_old == "NA"), NA,
                         if_else(any(driver_change_old == "S"), "S",
                         if_else(any(driver_change_old == "P") & any(driver_change_old == "R"), "S",
                         if_else(any(driver_change_old == "P"), "P", "R")))),
         protein_change= if_else(all(protein_change_old == "NA"), NA,
                                 if_else(any(protein_change_old == "p.R132H"), "R132H", "Other"))
           )%>%
  ungroup()%>%
  select(case_barcode, tumor_pair_barcode, idh_codel_subtype, HM_group, gene_symbol, driver_status, driver_change, protein_change)%>%
  distinct(.)

driver_changes <- rbind(snv_data, idh1_idh2_short)

write.table(driver_changes, file = "~/idh_mutations.txt", sep = "\t", row.names = FALSE)