#################################################################################################################
# Content: Multivariable overall survival stepwise-backward Cox PH models
# Author: Mircea Tesileanu
# Date: 2025.07.03
# R-version: 4.3.2
# Tables: Extended Data Table 2
#################################################################################################################

# Clean start and load libraries
rm(list = ls())
gc()

library(tidyverse)    # v2.0.0
library(RPostgres)    # v1.4.7
library(survival)     # v3.7-0
library(survminer)    # v0.5.0

#################################################################################################################

# Load relevant data
molecular.data <-  read.delim("/path/to/driver_changes_14042025.txt")%>%
  mutate(timing.recurrent = substr(tumor_pair_barcode, 20,21),
         timing.primary = substr(tumor_pair_barcode, 14,15) )%>%
  select(-idh_codel_subtype) 

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")

surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
  mutate(sample_barcode =   ifelse(surgery_number == 1, paste0(case_barcode, "-TP"),
                            ifelse(surgery_number == 2, paste0(case_barcode, "-R1"),
                            ifelse(surgery_number == 3, paste0(case_barcode, "-R2"),
                            ifelse(surgery_number == 4, paste0(case_barcode, "-R3"),
                            ifelse(surgery_number == 5, paste0(case_barcode, "-R4"),
                            ifelse(surgery_number == 6, paste0(case_barcode, "-R5"),
                            ifelse(surgery_number == 7, paste0(case_barcode, "-R6"),
                            ifelse(surgery_number == 8, paste0(case_barcode, "-R7"),
                            NA)))))))),
         timing.surgery = substr(sample_barcode, 14,15)) %>%
  filter(idh_codel_subtype != "IDHwt")

silver.set <- dbGetQuery(con, "SELECT * FROM analysis.silver_set")%>%
  semi_join(surgeries, by = "case_barcode")

survival.data <- dbGetQuery(con, "SELECT * FROM clinical.cases")[,c(2,5:8)]%>%
  semi_join(silver.set, by="case_barcode")

dbDisconnect(con)

full.data <- survival.data%>%
  left_join(surgeries, by ="case_barcode")%>%
  left_join(molecular.data, by ="case_barcode")%>%
  group_by(case_barcode)%>%
  mutate(age_at_diagnosis =               if_else(case_age_diagnosis_years >40, "40+", "40 or younger"),
         radiotherapy =                   if_else(surgery_number == 1,
                                                  case_when(treatment_radiotherapy == T ~ "Treated",
                                                            treatment_radiotherapy == F ~"Non-treated"),
                                                  NA),
         alkylating_chemotherapy=         if_else(surgery_number == 1,
                                                  case_when(treatment_alkylating_agent == T  &  (treatment_concurrent_tmz == T|treatment_concurrent_tmz == F|is.na(treatment_concurrent_tmz))~ "Treated",
                                                            treatment_concurrent_tmz == T  &  (treatment_alkylating_agent == T |treatment_alkylating_agent == F|is.na(treatment_alkylating_agent))~ "Treated",
                                                            treatment_alkylating_agent == F  &  (treatment_concurrent_tmz == F |  is.na(treatment_concurrent_tmz))~ "Non-treated",
                                                            is.na(treatment_alkylating_agent)  &  treatment_concurrent_tmz == F~ "Non-treated"),
                                                  NA),
         glioma.type =                    if_else(idh_codel_subtype != "IDHmut-codel", "IDH-mutant Astrocytomas", "IDH-mutant Oligodendrogliomas"),
         NOTCH1 =                         if_else(timing.primary == "TP" & surgery_number == 1,
                                                  if_else(gene_symbol == "NOTCH1" & driver_status == "mutant" & driver_change == "P", "present", 
                                                  if_else(gene_symbol == "NOTCH1" & driver_status == "mutant" & driver_change == "S", "present", "absent")),
                                                  NA),
         PIK3CA =                         if_else(timing.primary == "TP" & surgery_number == 1,
                                                  if_else(gene_symbol == "PIK3CA" & driver_status == "mutant" & driver_change == "P", "present", 
                                                  if_else(gene_symbol == "PIK3CA" & driver_status == "mutant" & driver_change == "S", "present", "absent")),
                                                  NA),
         CDKN2AB =                        if_else(timing.primary == "TP" & surgery_number == 1,
                                                  if_else(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "P", "present", 
                                                  if_else(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "S", "present", "absent")),
                                                  NA),
         CCND2 =                          if_else(timing.primary == "TP" & surgery_number == 1,
                                                  if_else(gene_symbol == "CCND2" & driver_status == "HLAMP" & driver_change == "P", "present", 
                                                  if_else(gene_symbol == "CCND2" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                  NA),
         PDGFRA =                         if_else(timing.primary == "TP" & surgery_number == 1,
                                                  if_else(gene_symbol == "PDGFRA" & driver_status == "HLAMP" & driver_change == "P", "present", 
                                                  if_else(gene_symbol == "PDGFRA" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                  NA),
         CDK46 =                          if_else(timing.primary == "TP" & surgery_number == 1,
                                                  if_else(gene_symbol == "CDK4/CDK6" & driver_status == "HLAMP" & driver_change == "P", "present", 
                                                  if_else(gene_symbol == "CDK4/CDK6" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                  NA),
         MYCMYCN =                        if_else(timing.primary == "TP" & surgery_number == 1,
                                                  if_else(gene_symbol == "MYC/MYCN" & driver_status == "HLAMP" & driver_change == "P", "present", 
                                                  if_else(gene_symbol == "MYC/MYCN" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                  NA),
         grade =                          if_else(surgery_number == 1,
                                                  case_when(CDKN2AB == "present"& idh_codel_subtype == "IDHmut-noncodel" & (grade == "II"|grade == "III"|grade == "IV"|is.na(grade)) ~ "4",
                                                            (CDKN2AB == "absent"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-noncodel" & grade == "II" ~ "2",
                                                            (CDKN2AB == "absent"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-noncodel" & grade == "III" ~ "3",
                                                            (CDKN2AB == "absent"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-noncodel" & grade == "IV" ~ "4",
                                                            (CDKN2AB == "absent"|CDKN2AB == "present"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-codel" & grade == "II" ~ "2",
                                                            (CDKN2AB == "absent"|CDKN2AB == "present"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-codel" & grade == "III" ~ "3"),
                                                  NA),
         os_status =                      if_else(case_vital_status == "dead", 1, 0),
         os_time =                        (case_overall_survival_mo)/12
  )%>%
  ungroup()

patient.lvl.data <- full.data%>%
  group_by(case_barcode)%>%
  arrange(surgery_number)%>%
  summarise(age_at_diagnosis =           case_when(any(age_at_diagnosis=="40+") ~ "40+", 
                                                   any(age_at_diagnosis=="40 or younger") ~"40 or younger"),
            radiotherapy =               case_when(any(radiotherapy=="Treated") ~ "Treated", 
                                                   any(radiotherapy=="Non-treated") ~"Non-treated"),
            alkylating_chemotherapy =    case_when(any(alkylating_chemotherapy=="Treated") ~ "Treated", 
                                                   any(alkylating_chemotherapy=="Non-treated") ~"Non-treated"),
            grade =                      case_when(any(grade=="2") ~ "2", 
                                                   any(grade=="3") ~"3",
                                                   any(grade=="4") ~"4"),
            glioma.type =                case_when(any(glioma.type=="IDH-mutant Oligodendrogliomas") ~ "IDH-mutant Oligodendrogliomas", 
                                                   any(glioma.type=="IDH-mutant Astrocytomas") ~"IDH-mutant Astrocytomas"),
            NOTCH1 =                     case_when(any(NOTCH1=="present") ~ "present", 
                                                   any(NOTCH1=="absent") ~"absent"),
            PIK3CA =                     case_when(any(PIK3CA=="present") ~ "present", 
                                                   any(PIK3CA=="absent") ~"absent"),
            CDKN2AB =                    case_when(any(CDKN2AB=="present") ~ "present", 
                                                   any(CDKN2AB=="absent") ~"absent"),
            CCND2 =                      case_when(any(CCND2=="present") ~ "present", 
                                                   any(CCND2=="absent") ~"absent"),
            PDGFRA =                     case_when(any(PDGFRA=="present") ~ "present", 
                                                   any(PDGFRA=="absent") ~"absent"),
            CDK46 =                      case_when(any(CDK46=="present") ~ "present", 
                                                   any(CDK46=="absent") ~"absent"),
            MYCMYCN =                    case_when(any(MYCMYCN=="present") ~ "present", 
                                                   any(MYCMYCN=="absent") ~"absent"),
            os_status =                  os_status[1],
            os_time =                    os_time[1]
  )%>%
  ungroup()%>%
  filter(!duplicated(case_barcode))

#################################################################################################################

# Multivariable overall survival analysis of patients with IDH-mutant astrocytomas
all.astros <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Astrocytomas")

coxph.astro<- coxph(Surv(os_time, os_status) ~ radiotherapy  + grade + NOTCH1 + PIK3CA+
                    CDKN2AB + CCND2 + PDGFRA + CDK46 + MYCMYCN, data = all.astros )
summary(coxph.astro)
cox.zph(coxph.astro)
concordance(coxph.astro)$concordance

coxph.astro.no.pdgfra <- coxph(Surv(os_time, os_status) ~ radiotherapy  + grade + NOTCH1 + PIK3CA+
                      CDKN2AB + CCND2 + CDK46 + MYCMYCN, data = all.astros )
step(coxph.astro.no.pdgfra, direction = c( "backward"), steps = 1000)
step.coxph.astro <- coxph(Surv(os_time, os_status) ~ radiotherapy + NOTCH1 + PIK3CA + CCND2+ CDK46, data = all.astros )
summary(step.coxph.astro)
cox.zph(step.coxph.astro)
concordance(step.coxph.astro)$concordance

#################################################################################################################

# Multivariable overall survival analysis of patients with IDH-mutant oligodendrogliomas
all.oligos <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Oligodendrogliomas") %>% 
  filter(!is.na(age_at_diagnosis) & !is.na(radiotherapy) & !is.na(alkylating_chemotherapy))

coxph.oligo <- coxph(Surv(os_time, os_status) ~ age_at_diagnosis + radiotherapy + alkylating_chemotherapy, data = all.oligos )
summary(coxph.oligo)
cox.zph(coxph.oligo)
concordance(coxph.oligo)$concordance

step(coxph.oligo, direction = c( "backward"), steps = 1000)
step.coxph.oligo <- coxph(Surv(os_time, os_status) ~ age_at_diagnosis + radiotherapy + alkylating_chemotherapy, data = all.oligos )
summary(step.coxph.oligo)
cox.zph(step.coxph.oligo)
concordance(step.coxph.oligo)$concordance

#################################################################################################################