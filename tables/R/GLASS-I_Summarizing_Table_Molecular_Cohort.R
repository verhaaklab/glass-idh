# GLASS-I: Summarizing Table Molecular Cohort
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(RPostgres)
library(tidyverse)

# load data
molecular.data <- read.delim("~/driver_changes.txt")%>%                         ## created with GLASS-I_Drivers_Table.R
  select(-idh_codel_subtype, -tumor_pair_barcode) 

signature.data <-  read_csv("~/sbs11_sbs119_id8_Rprop.csv")                     ## created with mutational_signature_pipeline.sh

sbs11.data <- signature.data %>%
  mutate(recurrence_fraction_SBS11 = SBS11_R)%>%
  select(case_barcode, recurrence_fraction_SBS11)

sbs119.data <- signature.data %>%
  mutate(recurrence_fraction_SBS119 = SBS119_R)%>%
  select(case_barcode, recurrence_fraction_SBS119)

id8.data <- signature.data %>%
  mutate(recurrence_fraction_ID8 = ID8_R)%>%
  select(case_barcode, recurrence_fraction_ID8)

chromothripsis.data <- read_csv("~/chromothripsis_shatterseek_samples.csv")%>%  ## created with GLASS-I_Chromothripsis.R
  group_by(case_barcode)%>%
  filter(n() == 2) %>%
  arrange(surgery_number)%>%
  summarise(chromothripsis = if_else(chromothripsis.any.conf[1] == 1 & chromothripsis.any.conf[2] ==0, "chrom-I",
                             if_else(chromothripsis.any.conf[1] == 1 & chromothripsis.any.conf[2] ==1, "chrom-S",
                             if_else(chromothripsis.any.conf[1] == 0 & chromothripsis.any.conf[2] ==1, "chrom-R",
                             if_else(chromothripsis.any.conf[1] == 0 & chromothripsis.any.conf[2] ==0, "non-chrom",
                                     NA))))
            )%>%
  ungroup()

prior.treatment <- read_csv("~/prior_treatment_gold_set.csv")%>%                ## created with GLASS-I_Prior_Treatment_Table.R
  group_by(case_barcode)%>%
  arrange(surgery_number)%>%
  summarise(radiotherapy_prior_to_recurrence = RT[2],
            temozolomide_prior_to_recurrence = TMZ[2],
            alkylating_chemotherapy_prior_to_recurrence = ALK_CHEMO[2],)%>%
  ungroup()

prior.publication.data <- read_csv("~/glass3_analysis_gold_set.csv")%>%         ## data available on Synapse (Synapse ID: syn17038081)
  mutate(GLASS_data_release = "2022")%>%
  select(tumor_pair_barcode, GLASS_data_release)

purity.coverage.data <- read_csv("~/idh_mut_gold_purity_coverage.csv")          ## created with create_idh_mut_gold_purity_coverage_from_files.R

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")               ## data available on Synapse (Synapse ID: syn17038081)

surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
  filter(idh_codel_subtype != "IDHwt")%>%
  mutate(timing.surgery = substr(sample_barcode, 14,15))

gold.set <- dbGetQuery(con, "SELECT * FROM analysis.gold_set")%>%
  mutate(timing.recurrent = substr(tumor_pair_barcode, 20,21),
         timing.primary = substr(tumor_pair_barcode, 14,15),
         sequencing_type = substr(tumor_pair_barcode, 27,29))%>%
  semi_join(surgeries, by = "case_barcode")

survival.data <- dbGetQuery(con, "SELECT * FROM clinical.cases")[,c(2,5:8)]%>%
  semi_join(surgeries, by="case_barcode")

tmb.data <-  dbGetQuery(con, "SELECT * FROM analysis.tumor_clinical_comparison_2024") %>%
  semi_join(gold.set, by = "tumor_pair_barcode")%>%
  mutate(initial_TMB = mutation_burden_a,
         recurrence_TMB = mutation_burden_b )%>%
  select(case_barcode, initial_TMB, recurrence_TMB)

dbDisconnect(con)

wgd.data <- read.delim("~/wgd_updated_20260509.tsv")%>%                         ## created with WGD_sequenza_script.R
  select(-case_barcode)%>%
  filter(aliquot_barcode %in% c(gold.set$tumor_barcode_a, gold.set$tumor_barcode_b))%>%
  mutate(sample_barcode = substr(aliquot_barcode, 1, 15))%>%
  left_join(surgeries, by = "sample_barcode")%>%
  group_by(case_barcode)%>%
  filter(n() == 2) %>%
  arrange(surgery_number)%>%
  summarise(WGD = if_else(WGD_final[1] == "Present" & WGD_final[2] == "Absent", "WGD-I",
                          if_else(WGD_final[1] == "Present" & WGD_final[2] == "Present", "WGD-S",
                                  if_else(WGD_final[1] == "Absent" & WGD_final[2] == "Present", "WGD-R",
                                          if_else(WGD_final[1] == "Absent" & WGD_final[2] == "Absent", "non-WGD",
                                                  NA))))
  )%>%
  ungroup()

ecDNA.data <- read.delim("~/glass5_ampsuite_20260317.tsv")%>%                   ## created with amplicon_architect_merge_glass5.R
  mutate(sample_barcode = substr(aliquot_barcode, 1, 15))%>%
  left_join(surgeries, by = "sample_barcode")%>%
  semi_join(gold.set, by = "case_barcode")%>%
  group_by(case_barcode)%>%
  filter(n() == 2) %>%
  arrange(surgery_number)%>%
  summarise(ecDNA = if_else(ecDNA[1] == "ecDNA" & ecDNA[2] == "no_ecDNA", "ecDNA-I",
                    if_else(ecDNA[1] == "ecDNA"  & ecDNA[2] == "ecDNA" , "ecDNA-S",
                    if_else(ecDNA[1] == "no_ecDNA" & ecDNA[2] == "ecDNA", "ecDNA-R",
                    if_else(ecDNA[1] == "no_ecDNA" & ecDNA[2] == "no_ecDNA", "non-ecDNA",
                            NA))))
  )%>%
  ungroup()

kataegis.data <- read_csv("~/glass5_kataegis_count_all_wgs.csv")%>%             ## created with kataegis_pipeline.sh
  mutate(sample_barcode = substr(aliquot_barcode, 1, 15))%>%
  left_join(surgeries, by = "sample_barcode")%>%
  semi_join(gold.set, by = "case_barcode")%>%
  group_by(case_barcode)%>%
  filter(n() == 2) %>%
  arrange(surgery_number)%>%
  summarise(kataegis = if_else(kataegis_event_count[1] != 0 & kataegis_event_count[2] ==0, "kataegis-I",
                       if_else(kataegis_event_count[1] != 0 & kataegis_event_count[2] != 0, "kataegis-S",
                       if_else(kataegis_event_count[1] == 0 & kataegis_event_count[2] != 0, "kataegis-R",
                       if_else(kataegis_event_count[1] == 0 & kataegis_event_count[2] ==0, "non-kataegis",
                               NA))))
  )%>%
  ungroup()

supp.data.molecular.cohort <- survival.data%>%
  right_join(gold.set, by ="case_barcode")%>%
  left_join(surgeries, by ="case_barcode")%>%
  left_join(molecular.data, by ="case_barcode")%>%
  left_join(sbs11.data, by ="case_barcode")%>%
  left_join(sbs119.data, by ="case_barcode")%>%
  left_join(id8.data, by ="case_barcode")%>%
  left_join(prior.treatment, by ="case_barcode")%>%
  left_join(prior.publication.data, by ="tumor_pair_barcode")%>%
  left_join(chromothripsis.data, by ="case_barcode")%>%
  left_join(wgd.data, by ="case_barcode")%>%
  left_join(ecDNA.data, by ="case_barcode")%>%
  left_join(kataegis.data, by ="case_barcode")%>%
  left_join(tmb.data, by ="case_barcode")%>%
  left_join(purity.coverage.data, by ="case_barcode")%>%
  group_by(case_barcode)%>%
  summarise(sequencing_platform =                                  sequencing_type[1],
            GLASS_data_release =                                   if_else(all(is.na(GLASS_data_release)), "2026", "2022"),
            glioma_type =                                          if_else(idh_codel_subtype[1] != "IDHmut-codel", "astrocytoma", "oligodendroglioma"),
            hypermutation_status =                                 if_else(HM_group[1] == "HM", "hypermutant",
                                                                           if_else(HM_group[1] == "NHM", "non-hypermutant", NA)),
            age_at_diagnosis =                                     case_age_diagnosis_years[1],
            patient_sex =                                          case_sex[1],
            os_status =                                            if_else(case_vital_status[1] == "dead", 1, 0),
            os_time =                                              round((case_overall_survival_mo[1])/12, digits = 6),
            surgery_number_initial =                               surgery_number[which(timing.surgery == timing.primary)][1],
            surgery_number_recurrence =                            surgery_number[which(timing.surgery == timing.recurrent)][1],
            radiotherapy_prior_to_recurrence =                     if_else(radiotherapy_prior_to_recurrence[1] == 1, "treated", "non-treated"),
            temozolomide_prior_to_recurrence =                     if_else(temozolomide_prior_to_recurrence[1] == 1, "treated", "non-treated"),
            alkylating_chemotherapy_prior_to_recurrence =          if_else(alkylating_chemotherapy_prior_to_recurrence[1] == 1, "treated", "non-treated"),
            
            initial_TMB =                                          initial_TMB[1],
            recurrence_TMB =                                       recurrence_TMB[1],
            driver_CNA_CDKN2AB =                                   if_else(any(gene_symbol == "CDKN2A/CDKN2B" & is.na(driver_change)),  "non-HD",
                                                                           if_else(any(gene_symbol == "CDKN2A/CDKN2B" & driver_change == "P"), "HD-I", 
                                                                                   if_else(any(gene_symbol == "CDKN2A/CDKN2B" & driver_change == "S"), "HD-S",
                                                                                           if_else(any(gene_symbol == "CDKN2A/CDKN2B" & driver_change == "R"), "HD-R",
                                                                                                   NA)))),
            driver_CNA_focal_oncogene_panel =                      if_else(any(gene_symbol == "Any Amplification" & is.na(driver_change)),  "non-amp",
                                                                           if_else(any(gene_symbol == "Any Amplification" & driver_change == "P"), "amp-I", 
                                                                                   if_else(any(gene_symbol == "Any Amplification" & driver_change == "S"), "amp-S",
                                                                                           if_else(any(gene_symbol == "Any Amplification" & driver_change == "R"), "amp-R",
                                                                                                   NA)))),
            driver_sSNV_PI3K =                                     if_else(any(gene_symbol == "PI3K" & is.na(driver_change)),  "wt",
                                                                           if_else(any(gene_symbol == "PI3K" & driver_change == "P"), "mut-I", 
                                                                                   if_else(any(gene_symbol == "PI3K" & driver_change == "S"), "mut-S",
                                                                                           if_else(any(gene_symbol == "PI3K" & driver_change == "R"), "mut-R",
                                                                                                   NA)))),
            driver_sSNV_NOTCH1 =                                   if_else(any(gene_symbol == "NOTCH1" & is.na(driver_change)),  "wt",
                                                                           if_else(any(gene_symbol == "NOTCH1" & driver_change == "P"), "mut-I", 
                                                                                   if_else(any(gene_symbol == "NOTCH1" & driver_change == "S"), "mut-S",
                                                                                           if_else(any(gene_symbol == "NOTCH1" & driver_change == "R"), "mut-R",
                                                                                                   NA)))),
            driver_sSNV_MSH6 =                                     if_else(any(gene_symbol == "MSH6" & is.na(driver_change)),  "wt",
                                                                           if_else(any(gene_symbol == "MSH6" & driver_change == "P"), "mut-I", 
                                                                                   if_else(any(gene_symbol == "MSH6" & driver_change == "S"), "mut-S",
                                                                                           if_else(any(gene_symbol == "MSH6" & driver_change == "R"), "mut-R",
                                                                                                   NA)))),
            recurrence_fraction_SBS11 =                            recurrence_fraction_SBS11[1],
            recurrence_fraction_SBS119 =                           recurrence_fraction_SBS119[1],
            recurrence_fraction_ID8 =                              recurrence_fraction_ID8[1],
            chromothripsis =                                       chromothripsis[1],
            WGD =                                                  WGD[1],
            kataegis =                                             kataegis[1],
            ecDNA =                                                ecDNA[1],
            tumor_purity_initial =                                 seqz_initial[1],
            tumor_purity_recurrence =                              seqz_recurrent[1],            
            coverage_mean_initial =                                round((cov_mean_initial[1]), digits = 6),
            coverage_mean_recurrence =                             round((cov_mean_recurrent[1]), digits = 6),
  )%>%   
  ungroup()%>% 
  arrange(sequencing_platform, glioma_type, hypermutation_status, GLASS_data_release)

write.csv(supp.data.molecular.cohort, "~/summary_molecular_cohort.csv")