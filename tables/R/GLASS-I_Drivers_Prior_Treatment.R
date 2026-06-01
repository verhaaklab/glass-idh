# GLASS-I: Correlations Driver sSNVs/CNAs To Prior Treatment 
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(RPostgres)
library(tidyverse)
library(gtsummary)
library(gt)

# load data
molecular.data <- read.delim("~/driver_changes.txt")%>%                       ## created with GLASS-I_Drivers_Table.R
  mutate(timing.recurrent = substr(tumor_pair_barcode, 20,21),
         timing.primary = substr(tumor_pair_barcode, 14,15) )%>%
  select(-idh_codel_subtype) 

prior.treatment <- read_csv("~/prior_treatment_gold_set.csv")%>%              ## created with GLASS-I_Prior_Treatment_Table.R
  select(-case_barcode, -surgery_number)

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")             ## data available on Synapse (Synapse ID: syn17038081)

surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
  mutate(timing.surgery = substr(sample_barcode, 14,15)) %>%
  filter(idh_codel_subtype != "IDHwt")%>%
  left_join(prior.treatment, by ="sample_barcode")

dbDisconnect(con)

driver.at.recurrence <- surgeries%>%
  left_join(molecular.data, by ="case_barcode")%>%
  filter(timing.recurrent == timing.surgery) %>%  
  group_by(case_barcode)%>%
  summarise(HM_group =         HM_group[1],
            RT =               if_else(RT[1] ==1, "Radiotherapy", "No radiotherapy"),
            TMZ =              TMZ[1],
            ALK_CHEMO =        if_else(ALK_CHEMO[1] ==1, "Alkylating chemotherapy", "No alkylating chemotherapy"),
            tumor_type =       if_else(idh_codel_subtype[1] != "IDHmut-codel", "IDH-mutant Astrocytomas", "IDH-mutant Oligodendrogliomas"),
            NOTCH1 =           if_else(any(gene_symbol == "NOTCH1" & driver_status == "mutant" & driver_change == "R"), 1, 
                                       if_else(any(gene_symbol == "NOTCH1" & driver_status == "mutant" & driver_change == "S"), 1, 0)),
            PI3K =             if_else(any(gene_symbol == "PI3K" & driver_status == "mutant" & driver_change == "R"), 1, 
                                       if_else(any(gene_symbol == "PI3K" & driver_status == "mutant" & driver_change == "S"), 1, 0)),
            MSH6 =             if_else(any(gene_symbol == "MSH6" & driver_status == "mutant" & driver_change == "R"), 1, 
                                       if_else(any(gene_symbol == "MSH6" & driver_status == "mutant" & driver_change == "S"), 1, 0)),
            any_amp =          if_else(any(gene_symbol == "Any Amplification" & driver_status == "HLAMP" & driver_change == "R"), 1, 
                                       if_else(any(gene_symbol == "Any Amplification" & driver_status == "HLAMP" & driver_change == "S"), 1, 0)),
            CDKN2AB =          if_else(any(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "R"), 1, 
                                       if_else(any(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "S"), 1, 0)),
            CDKN2AB.R.only =          if_else(any(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "R"), 1, 
                                              if_else(any(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "S"), NA, 0)),
  )%>%
  ungroup()%>%
  mutate(ANY_Rx = if_else(is.na(ALK_CHEMO) | is.na(RT), NA,
                          if_else(ALK_CHEMO ==  "Alkylating chemotherapy" | RT == "Radiotherapy", "Any treatment", "No treatment")))

astrocytomas.hm <- driver.at.recurrence%>%
  filter(tumor_type == "IDH-mutant Astrocytomas" & HM_group == "HM")

astrocytomas.nhm <- driver.at.recurrence%>%
  filter(tumor_type == "IDH-mutant Astrocytomas" & HM_group == "NHM")

oligodendrogliomas.hm <- driver.at.recurrence%>%
  filter(tumor_type == "IDH-mutant Oligodendrogliomas" & HM_group == "HM")

oligodendrogliomas.nhm <- driver.at.recurrence%>%
  filter(tumor_type == "IDH-mutant Oligodendrogliomas" & HM_group == "NHM")

oligodendrogliomas.hm %>%
  tbl_summary(include = c(
    RT,
    MSH6,
    NOTCH1,
    PI3K,
    any_amp,
    CDKN2AB
  ),
  by = RT,
  missing = "no",
  type = all_continuous() ~ "continuous2",
  statistic = list(
    all_categorical() ~ "{p}% ({n}/{N})"
  ))%>% 
  add_p()%>%
  add_q()

oligodendrogliomas.nhm %>%
  tbl_summary(include = c(
    RT,
    MSH6,
    NOTCH1,
    PI3K,
    any_amp,
    CDKN2AB
  ),
  by = RT,
  missing = "no",
  type = all_continuous() ~ "continuous2",
  statistic = list(
    all_categorical() ~ "{p}% ({n}/{N})"
  ))%>% 
  add_p()%>%
  add_q()

astrocytomas.hm %>%
  tbl_summary(include = c(
    RT,
    MSH6,
    NOTCH1,
    PI3K,
    any_amp,
    CDKN2AB
  ),
  by = RT,
  missing = "no",
  type = all_continuous() ~ "continuous2",
  statistic = list(
    all_categorical() ~ "{p}% ({n}/{N})"
  ))%>% 
  add_p()%>%
  add_q()

astrocytomas.nhm %>%
  tbl_summary(include = c(
    RT,
    MSH6,
    NOTCH1,
    PI3K,
    any_amp,
    CDKN2AB
  ),
  by = RT,
  missing = "no",
  type = all_continuous() ~ "continuous2",
  statistic = list(
    all_categorical() ~ "{p}% ({n}/{N})"
  ))%>% 
  add_p()%>%
  add_q()

oligodendrogliomas.hm %>%
  tbl_summary(include = c(
    ALK_CHEMO,
    MSH6,
    NOTCH1,
    PI3K,
    any_amp,
    CDKN2AB
  ),
  by = ALK_CHEMO,
  missing = "no",
  type = all_continuous() ~ "continuous2",
  statistic = list(
    all_categorical() ~ "{p}% ({n}/{N})"
  ))%>% 
  add_p()%>%
  add_q()

oligodendrogliomas.nhm %>%
  tbl_summary(include = c(
    ALK_CHEMO,
    MSH6,
    NOTCH1,
    PI3K,
    any_amp,
    CDKN2AB
  ),
  by = ALK_CHEMO,
  missing = "no",
  type = all_continuous() ~ "continuous2",
  statistic = list(
    all_categorical() ~ "{p}% ({n}/{N})"
  ))%>% 
  add_p()%>%
  add_q()

astrocytomas.hm %>%
  tbl_summary(include = c(
    ALK_CHEMO,
    MSH6,
    NOTCH1,
    PI3K,
    any_amp,
    CDKN2AB
  ),
  by = ALK_CHEMO,
  missing = "no",
  type = all_continuous() ~ "continuous2",
  statistic = list(
    all_categorical() ~ "{p}% ({n}/{N})"
  ))%>% 
  add_p()%>%
  add_q()

astrocytomas.nhm %>%
  tbl_summary(include = c(
    ALK_CHEMO,
    MSH6,
    NOTCH1,
    PI3K,
    any_amp,
    CDKN2AB
  ),
  by = ALK_CHEMO,
  missing = "no",
  type = all_continuous() ~ "continuous2",
  statistic = list(
    all_categorical() ~ "{p}% ({n}/{N})"
  ))%>% 
  add_p()%>%
  add_q()

oligodendrogliomas.hm %>%
  tbl_summary(include = c(
    ANY_Rx,
    MSH6,
    NOTCH1,
    PI3K,
    any_amp,
    CDKN2AB
  ),
  by = ANY_Rx,
  missing = "no",
  type = all_continuous() ~ "continuous2",
  statistic = list(
    all_categorical() ~ "{p}% ({n}/{N})"
  ))%>% 
  add_p()%>%
  add_q()

oligodendrogliomas.nhm %>%
  tbl_summary(include = c(
    ANY_Rx,
    MSH6,
    NOTCH1,
    PI3K,
    any_amp,
    CDKN2AB
  ),
  by = ANY_Rx,
  missing = "no",
  type = all_continuous() ~ "continuous2",
  statistic = list(
    all_categorical() ~ "{p}% ({n}/{N})"
  ))%>% 
  add_p()%>%
  add_q()

astrocytomas.hm %>%
  tbl_summary(include = c(
    ANY_Rx,
    MSH6,
    NOTCH1,
    PI3K,
    any_amp,
    CDKN2AB
  ),
  by = ANY_Rx,
  missing = "no",
  type = all_continuous() ~ "continuous2",
  statistic = list(
    all_categorical() ~ "{p}% ({n}/{N})"
  ))%>% 
  add_p()%>%
  add_q()

astrocytomas.nhm %>%
  tbl_summary(include = c(
    ANY_Rx,
    MSH6,
    NOTCH1,
    PI3K,
    any_amp,
    CDKN2AB
  ),
  by = ANY_Rx,
  missing = "no",
  type = all_continuous() ~ "continuous2",
  statistic = list(
    all_categorical() ~ "{p}% ({n}/{N})"
  ))%>% 
  add_p()%>%
  add_q()