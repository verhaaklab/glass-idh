# GLASS-I: Drivers Table
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

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")                  ## data available on Synapse (Synapse ID: syn17038081)
hm_info <- dbGetQuery(con, "SELECT * FROM analysis.tumor_clinical_comparison_2024")%>%
  mutate(HM_group = ifelse(sbs11_r_hypermutation==1,"HM","NHM"))%>%
  select(case_barcode, tumor_pair_barcode, HM_group)

df_gs <- dbGetQuery(con, "SELECT * FROM analysis.gold_set gs JOIN clinical.subtypes cs ON cs.case_barcode = gs.case_barcode WHERE cs.idh_codel_subtype NOT LIKE 'IDHwt'")%>%
  select(c(1,2,6))

cnv_data <- dbGetQuery(con, read_file("~/GLASS-I_CNA_Data.sql"))%>%              ## see SQL script of same name
  left_join(hm_info, by = c("case_barcode","tumor_pair_barcode"))%>%
  mutate(driver_status = if_else(is.na(cnv_state), "No_AMP/HD", cnv_state) ,
         driver_change = cnv_change)%>%
  select(case_barcode, tumor_pair_barcode, idh_codel_subtype, HM_group, gene_symbol, driver_status, driver_change)

## counting frequency table for PIK3CA
pik3ca_short <- dbGetQuery(con, "SELECT pg.* from variants.pgeno pg JOIN analysis.gold_set gs ON pg.tumor_pair_barcode = gs.tumor_pair_barcode WHERE pg.gene_symbol = 'PIK3CA'")%>%
  semi_join(df_gs, by = c("case_barcode","tumor_pair_barcode"))%>% 
  mutate(driver_change = case_when(mutect2_call_a==1 & mutect2_call_b==0 ~ "P", mutect2_call_a==0 & mutect2_call_b==1 ~ "R", mutect2_call_a==1 & mutect2_call_b==1 ~ "S", TRUE ~ "none"))%>% 
  filter(!driver_change=="none" & variant_classification %in% variant_list)%>%
  full_join(df_gs, by = c("case_barcode","tumor_pair_barcode"))%>% 
  mutate(driver_status = ifelse(is.na(gene_symbol),"wildtype","mutant"))%>%
  left_join(hm_info, by = c("case_barcode","tumor_pair_barcode"))%>% 
  group_by(case_barcode)%>% 
  mutate(multihit = ifelse(n()>1, "yes","no"))%>% 
  slice_max(order_by = af_b, n = 1)%>%
  slice(1)%>% 
  ungroup()%>%
  mutate(gene_symbol = "PIK3CA")%>%
  select(case_barcode, tumor_pair_barcode, idh_codel_subtype, HM_group, gene_symbol, driver_status, driver_change)

## counting frequency table for PIK3R1
pik3r1_short <- dbGetQuery(con, "SELECT pg.* from variants.pgeno pg JOIN analysis.gold_set gs ON pg.tumor_pair_barcode = gs.tumor_pair_barcode WHERE pg.gene_symbol = 'PIK3R1'")%>%
  semi_join(df_gs, by = c("case_barcode","tumor_pair_barcode"))%>% 
  mutate(driver_change = case_when(mutect2_call_a==1 & mutect2_call_b==0 ~ "P", mutect2_call_a==0 & mutect2_call_b==1 ~ "R", mutect2_call_a==1 & mutect2_call_b==1 ~ "S", TRUE ~ "none"))%>% 
  filter(!driver_change=="none" & variant_classification %in% variant_list)%>%
  full_join(df_gs, by = c("case_barcode","tumor_pair_barcode"))%>% 
  mutate(driver_status = ifelse(is.na(gene_symbol),"wildtype","mutant"))%>%
  left_join(hm_info, by = c("case_barcode","tumor_pair_barcode"))%>% 
  group_by(case_barcode)%>% 
  mutate(multihit = ifelse(n()>1, "yes","no"))%>% 
  slice_max(order_by = af_b, n = 1)%>%
  slice(1)%>% 
  ungroup()%>%
  mutate(gene_symbol = "PIK3R1")%>%
  select(case_barcode, tumor_pair_barcode, idh_codel_subtype, HM_group, gene_symbol, driver_status, driver_change)

## counting frequency table for NOTCH1
notch1_short <- dbGetQuery(con, "SELECT pg.* from variants.pgeno pg JOIN analysis.gold_set gs ON pg.tumor_pair_barcode = gs.tumor_pair_barcode WHERE pg.gene_symbol = 'NOTCH1'")%>%
  semi_join(df_gs, by = c("case_barcode","tumor_pair_barcode"))%>% 
  mutate(driver_change = case_when(mutect2_call_a==1 & mutect2_call_b==0 ~ "P", mutect2_call_a==0 & mutect2_call_b==1 ~ "R", mutect2_call_a==1 & mutect2_call_b==1 ~ "S", TRUE ~ "none"))%>% 
  filter(!driver_change=="none" & variant_classification %in% variant_list)%>%
  full_join(df_gs, by = c("case_barcode","tumor_pair_barcode"))%>% 
  mutate(driver_status = ifelse(is.na(gene_symbol),"wildtype","mutant"))%>%
  left_join(hm_info, by = c("case_barcode","tumor_pair_barcode"))%>% 
  group_by(case_barcode)%>% 
  mutate(multihit = ifelse(n()>1, "yes","no"))%>% 
  slice_max(order_by = af_b, n = 1)%>%
  slice(1)%>% 
  ungroup()%>%
  mutate(gene_symbol = "NOTCH1")%>%
  select(case_barcode, tumor_pair_barcode, idh_codel_subtype, HM_group, gene_symbol, driver_status, driver_change)

## counting frequency table for NOTCH1
notch2_short <- dbGetQuery(con, "SELECT pg.* from variants.pgeno pg JOIN analysis.gold_set gs ON pg.tumor_pair_barcode = gs.tumor_pair_barcode WHERE pg.gene_symbol = 'NOTCH2'")%>%
  semi_join(df_gs, by = c("case_barcode","tumor_pair_barcode"))%>% 
  mutate(driver_change = case_when(mutect2_call_a==1 & mutect2_call_b==0 ~ "P", mutect2_call_a==0 & mutect2_call_b==1 ~ "R", mutect2_call_a==1 & mutect2_call_b==1 ~ "S", TRUE ~ "none"))%>% 
  filter(!driver_change=="none" & variant_classification %in% variant_list)%>%
  full_join(df_gs, by = c("case_barcode","tumor_pair_barcode"))%>% 
  mutate(driver_status = ifelse(is.na(gene_symbol),"wildtype","mutant"))%>%
  left_join(hm_info, by = c("case_barcode","tumor_pair_barcode"))%>% 
  group_by(case_barcode)%>% 
  mutate(multihit = ifelse(n()>1, "yes","no"))%>% 
  slice_max(order_by = af_b, n = 1)%>%
  slice(1)%>% 
  ungroup()%>%
  mutate(gene_symbol = "NOTCH2")%>%
  select(case_barcode, tumor_pair_barcode, idh_codel_subtype, HM_group, gene_symbol, driver_status, driver_change)

## counting frequency table for NOTCH3
notch3_short <- dbGetQuery(con, "SELECT pg.* from variants.pgeno pg JOIN analysis.gold_set gs ON pg.tumor_pair_barcode = gs.tumor_pair_barcode WHERE pg.gene_symbol = 'NOTCH3'")%>%
  semi_join(df_gs, by = c("case_barcode","tumor_pair_barcode"))%>% 
  mutate(driver_change = case_when(mutect2_call_a==1 & mutect2_call_b==0 ~ "P", mutect2_call_a==0 & mutect2_call_b==1 ~ "R", mutect2_call_a==1 & mutect2_call_b==1 ~ "S", TRUE ~ "none"))%>% 
  filter(!driver_change=="none" & variant_classification %in% variant_list)%>%
  full_join(df_gs, by = c("case_barcode","tumor_pair_barcode"))%>% 
  mutate(driver_status = ifelse(is.na(gene_symbol),"wildtype","mutant"))%>%
  left_join(hm_info, by = c("case_barcode","tumor_pair_barcode"))%>% 
  group_by(case_barcode)%>% 
  mutate(multihit = ifelse(n()>1, "yes","no"))%>% 
  slice_max(order_by = af_b, n = 1)%>%
  slice(1)%>% 
  ungroup()%>%
  mutate(gene_symbol = "NOTCH3")%>%
  select(case_barcode, tumor_pair_barcode, idh_codel_subtype, HM_group, gene_symbol, driver_status, driver_change)

## counting frequency table for MSH6
msh6_short <- dbGetQuery(con, "SELECT pg.* from variants.pgeno pg JOIN analysis.gold_set gs ON pg.tumor_pair_barcode = gs.tumor_pair_barcode WHERE pg.gene_symbol = 'MSH6'")%>%
  semi_join(df_gs, by = c("case_barcode","tumor_pair_barcode"))%>% 
  mutate(driver_change = case_when(mutect2_call_a==1 & mutect2_call_b==0 ~ "P", mutect2_call_a==0 & mutect2_call_b==1 ~ "R", mutect2_call_a==1 & mutect2_call_b==1 ~ "S", TRUE ~ "none"))%>% 
  filter(!driver_change=="none" & variant_classification %in% variant_list)%>%
  full_join(df_gs, by = c("case_barcode","tumor_pair_barcode"))%>% 
  mutate(driver_status = ifelse(is.na(gene_symbol),"wildtype","mutant"))%>%
  left_join(hm_info, by = c("case_barcode","tumor_pair_barcode"))%>% 
  group_by(case_barcode)%>% 
  mutate(multihit = ifelse(n()>1, "yes","no"))%>% 
  slice_max(order_by = af_b, n = 1)%>%
  slice(1)%>% 
  ungroup()%>%
  mutate(gene_symbol = "MSH6")%>%
  select(case_barcode, tumor_pair_barcode, idh_codel_subtype, HM_group, gene_symbol, driver_status, driver_change)

snv_data <- rbind(msh6_short, notch1_short, notch2_short, notch3_short, pik3ca_short, pik3r1_short)

## counting frequency table for CNV drivers
## counting and combining frequency table for MET/PTPRZ1
met_ptprz1_short <- subset(cnv_data, cnv_data$gene_symbol=="MET"|cnv_data$gene_symbol=="PTPRZ1")%>%
  mutate(gene_symbol_old = gene_symbol,
         gene_symbol = "MET/PTPRZ1",
         driver_status_old = driver_status,
         driver_change_old = if_else(is.na(driver_change), "NA", driver_change))%>% 
  group_by(case_barcode) %>% 
  mutate(driver_status = if_else(any(driver_status_old == "HLAMP"), "HLAMP", "No_AMP/HD"),
         driver_change = if_else(all(driver_change_old == "NA"), NA,
                         if_else(any(driver_change_old == "S"), "S",
                         if_else(any(driver_change_old == "P") & any(driver_change_old == "R"), "S",
                                         if_else(any(driver_change_old == "P"), "P", "R")))))%>%
  ungroup()%>%
  select(case_barcode, tumor_pair_barcode, idh_codel_subtype, HM_group, gene_symbol, driver_status, driver_change)%>%
  distinct(.)

## counting and combining frequency table for CDK4/CDK6
cdk4_cdk6_short <- subset(cnv_data, cnv_data$gene_symbol=="CDK4"|cnv_data$gene_symbol=="CDK6")%>%
  mutate(gene_symbol_old = gene_symbol,
         gene_symbol = "CDK4/CDK6",
         driver_status_old = driver_status,
         driver_change_old = if_else(is.na(driver_change), "NA", driver_change))%>% 
  group_by(case_barcode) %>% 
  mutate(driver_status = if_else(any(driver_status_old == "HLAMP"), "HLAMP", "No_AMP/HD"),
         driver_change = if_else(all(driver_change_old == "NA"), NA,
                                 if_else(any(driver_change_old == "S"), "S",
                                         if_else(any(driver_change_old == "P") & any(driver_change_old == "R"), "S",
                                                 if_else(any(driver_change_old == "P"), "P", "R")))))%>%
  ungroup()%>%
  select(case_barcode, tumor_pair_barcode, idh_codel_subtype, HM_group, gene_symbol, driver_status, driver_change)%>%
  distinct(.)

## counting and combining frequency table for MYC/MYCN
myc_mycn_short <- subset(cnv_data, cnv_data$gene_symbol=="MYC"|cnv_data$gene_symbol=="MYCN")%>%
  mutate(gene_symbol_old = gene_symbol,
         gene_symbol = "MYC/MYCN",
         driver_status_old = driver_status,
         driver_change_old = if_else(is.na(driver_change), "NA", driver_change))%>% 
  group_by(case_barcode) %>% 
  mutate(driver_status = if_else(any(driver_status_old == "HLAMP"), "HLAMP", "No_AMP/HD"),
         driver_change = if_else(all(driver_change_old == "NA"), NA,
                                 if_else(any(driver_change_old == "S"), "S",
                                         if_else(any(driver_change_old == "P") & any(driver_change_old == "R"), "S",
                                                 if_else(any(driver_change_old == "P"), "P", "R")))))%>%
  ungroup()%>%
  select(case_barcode, tumor_pair_barcode, idh_codel_subtype, HM_group, gene_symbol, driver_status, driver_change)%>%
  distinct(.)

## counting and combining frequency table for all amplifications
any_amp_short <- subset(cnv_data, cnv_data$gene_symbol=="MYC"|cnv_data$gene_symbol=="MYCN"|
                          cnv_data$gene_symbol=="CDK4"|cnv_data$gene_symbol=="CDK6"|
                          cnv_data$gene_symbol=="MET"|cnv_data$gene_symbol=="PTPRZ1"|
                          cnv_data$gene_symbol=="CCND2"|cnv_data$gene_symbol=="EGFR"|
                          cnv_data$gene_symbol=="PDGFRA")%>%
  mutate(gene_symbol_old = gene_symbol,
         gene_symbol = "Any Amplification",
         driver_status_old = driver_status,
         driver_change_old = if_else(is.na(driver_change), "NA", driver_change))%>% 
  group_by(case_barcode) %>% 
  mutate(driver_status = if_else(any(driver_status_old == "HLAMP"), "HLAMP", "No_AMP/HD"),
         driver_change = if_else(all(driver_change_old == "NA"), NA,
                                 if_else(any(driver_change_old == "S"), "S",
                                         if_else(any(driver_change_old == "P") & any(driver_change_old == "R"), "S",
                                                 if_else(any(driver_change_old == "P"), "P", "R")))))%>%
  ungroup()%>%
  select(case_barcode, tumor_pair_barcode, idh_codel_subtype, HM_group, gene_symbol, driver_status, driver_change)%>%
  distinct(.)

## counting and combining frequency table for CDKN2A/CDKN2B
cdkn2a_cdkn2b_short <- subset(cnv_data, cnv_data$gene_symbol=="CDKN2A"|cnv_data$gene_symbol=="CDKN2B")%>%
  mutate(gene_symbol_old = gene_symbol,
         gene_symbol = "CDKN2A/CDKN2B",
         driver_status_old = driver_status,
         driver_change_old = if_else(is.na(driver_change), "NA", driver_change))%>% 
  group_by(case_barcode) %>% 
  mutate(driver_status = if_else(any(driver_status_old == "HLDEL"), "HLDEL", "No_AMP/HD"),
         driver_change = if_else(all(driver_change_old == "NA"), NA,
                                 if_else(any(driver_change_old == "S"), "S",
                                         if_else(any(driver_change_old == "P") & any(driver_change_old == "R"), "S",
                                                 if_else(any(driver_change_old == "P"), "P", "R")))))%>%
  ungroup()%>%
  select(case_barcode, tumor_pair_barcode, idh_codel_subtype, HM_group, gene_symbol, driver_status, driver_change)%>%
  distinct(.)

## counting and combining frequency table for all amplifications
any_del_short <- subset(cnv_data, cnv_data$gene_symbol=="CDKN2A"|cnv_data$gene_symbol=="CDKN2B"|
                          cnv_data$gene_symbol=="ATRX"|cnv_data$gene_symbol=="PTEN")%>%
  mutate(gene_symbol_old = gene_symbol,
         gene_symbol = "Any HD",
         driver_status_old = driver_status,
         driver_change_old = if_else(is.na(driver_change), "NA", driver_change))%>% 
  group_by(case_barcode) %>% 
  mutate(driver_status = if_else(any(driver_status_old == "HLDEL"), "HLDEL", "No_AMP/HD"),
         driver_change = if_else(all(driver_change_old == "NA"), NA,
                                 if_else(any(driver_change_old == "S"), "S",
                                         if_else(any(driver_change_old == "P") & any(driver_change_old == "R"), "S",
                                                 if_else(any(driver_change_old == "P"), "P", "R")))))%>%
  ungroup()%>%
  select(case_barcode, tumor_pair_barcode, idh_codel_subtype, HM_group, gene_symbol, driver_status, driver_change)%>%
  distinct(.)

## counting and combining frequency table for PI3K
any_snv_short <- subset(snv_data, snv_data$gene_symbol=="PIK3CA"|snv_data$gene_symbol=="PIK3R1"|
                          snv_data$gene_symbol=="NOTCH1"|snv_data$gene_symbol=="MSH6")%>%
  mutate(gene_symbol_old = gene_symbol,
         gene_symbol = "Any SNV",
         driver_status_old = driver_status,
         driver_change_old = if_else(is.na(driver_change), "NA", driver_change))%>% 
  group_by(case_barcode) %>% 
  mutate(driver_status = if_else(any(driver_status_old == "mutant"), "mutant", "wildtype"),
         driver_change = if_else(all(driver_change_old == "NA"), NA,
                                 if_else(any(driver_change_old == "S"), "S",
                                         if_else(any(driver_change_old == "P") & any(driver_change_old == "R"), "S",
                                                 if_else(any(driver_change_old == "P"), "P", "R")))))%>%
  ungroup()%>%
  select(case_barcode, tumor_pair_barcode, idh_codel_subtype, HM_group, gene_symbol, driver_status, driver_change)%>%
  distinct(.)

## counting and combining frequency table for PI3K
pik3ca_pik3r1_short <- subset(snv_data, snv_data$gene_symbol=="PIK3CA"|snv_data$gene_symbol=="PIK3R1")%>%
  mutate(gene_symbol_old = gene_symbol,
         gene_symbol = "PI3K",
         driver_status_old = driver_status,
         driver_change_old = if_else(is.na(driver_change), "NA", driver_change))%>% 
  group_by(case_barcode) %>% 
  mutate(driver_status = if_else(any(driver_status_old == "mutant"), "mutant", "wildtype"),
         driver_change = if_else(all(driver_change_old == "NA"), NA,
                                 if_else(any(driver_change_old == "S"), "S",
                                         if_else(any(driver_change_old == "P") & any(driver_change_old == "R"), "S",
                                                 if_else(any(driver_change_old == "P"), "P", "R")))))%>%
  ungroup()%>%
  select(case_barcode, tumor_pair_barcode, idh_codel_subtype, HM_group, gene_symbol, driver_status, driver_change)%>%
  distinct(.)

driver_changes <- rbind(snv_data, pik3ca_pik3r1_short, cnv_data, cdkn2a_cdkn2b_short, cdk4_cdk6_short, met_ptprz1_short, myc_mycn_short, any_amp_short, any_snv_short, any_del_short)
write.table(driver_changes, file = "~/driver_changes.txt", sep = "\t", row.names = FALSE)