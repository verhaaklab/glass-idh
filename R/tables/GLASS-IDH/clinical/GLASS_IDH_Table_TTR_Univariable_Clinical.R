#################################################################################################################
# Content: Univariable time-to-recurrence Cox PH models of clinical variables
# Author: Mircea Tesileanu
# Date: 2025.07.03
# R-version: 4.3.2
# Tables: Extended Data Table 3
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
         timing.primary =   substr(tumor_pair_barcode, 14,15) )%>%
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
         sex =                            if_else(case_sex == "female", "Female", "Male"),
         surgery =                        if_else(surgery_number == 1,
                                                  case_when(surgery_type == "Biopsy" & (surgery_extent_of_resection == "Biopsy"|is.na(surgery_extent_of_resection)) ~"Biopsy",
                                                            surgery_type == "Craniotomy" & (surgery_extent_of_resection == "Subtotal"|surgery_extent_of_resection == "Total"|is.na(surgery_extent_of_resection)) ~ "Resection",
                                                            surgery_type == "Craniotomy" & surgery_extent_of_resection == "Biopsy"  ~ "Biopsy",
                                                            is.na(surgery_type) & surgery_extent_of_resection == "Biopsy" ~ "Biopsy",
                                                            is.na(surgery_type) & (surgery_extent_of_resection == "Subtotal"|surgery_extent_of_resection == "Total") ~ "Resection"),
                                                  NA),
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
         CDKN2AB =                        if_else(timing.primary == "TP" & surgery_number == 1,
                                                  if_else(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "P", "present", 
                                                  if_else(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "S", "present", "absent")),
                                                  NA),
         grade =                          if_else(surgery_number == 1,
                                                  case_when(CDKN2AB == "present"& idh_codel_subtype == "IDHmut-noncodel" & (grade == "II"|grade == "III"|grade == "IV"|is.na(grade)) ~ "4",
                                                            (CDKN2AB == "absent"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-noncodel" & grade == "II" ~ "2",
                                                            (CDKN2AB == "absent"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-noncodel" & grade == "III" ~ "3",
                                                            (CDKN2AB == "absent"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-noncodel" & grade == "IV" ~ "4",
                                                            (CDKN2AB == "absent"|CDKN2AB == "present"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-codel" & grade == "II" ~ "2",
                                                            (CDKN2AB == "absent"|CDKN2AB == "present"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-codel" & grade == "III" ~ "3"),
                                                  NA),
         ttr_status =                     1,
         ttr_time =                       if_else(surgery_number == 2 & surgical_interval_mo !=0,
                                                  surgical_interval_mo/12,
                                                  NA)
  )%>%
  ungroup()

patient.lvl.data <- full.data%>%
  group_by(case_barcode)%>%
  arrange(surgery_number)%>%
  summarise(age_at_diagnosis =           age_at_diagnosis[1],
            glioma.type =                glioma.type[1],
            sex =                        sex[1],
            surgery =                    surgery[1],
            radiotherapy =               radiotherapy[1],
            alkylating_chemotherapy =    alkylating_chemotherapy[1],
            grade =                      grade[1],
            ttr_status =                 ttr_status[1],
            ttr_time =                   ttr_time[surgery_number == 2]
  )%>%
  ungroup()%>%
  filter(!duplicated(case_barcode))

#################################################################################################################

# Univariable time-to-recurrence analysis of patients with IDH-mutant astrocytomas
all.astros <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Astrocytomas")

summary(coxph(Surv(ttr_time, ttr_status) ~ age_at_diagnosis, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ age_at_diagnosis, data = all.astros))
survfit(Surv(ttr_time, ttr_status) ~ age_at_diagnosis, data = all.astros )

summary(coxph(Surv(ttr_time, ttr_status) ~ sex, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ sex, data = all.astros))
survfit(Surv(ttr_time, ttr_status) ~ sex, data = all.astros )

summary(coxph(Surv(ttr_time, ttr_status) ~ surgery, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ surgery, data = all.astros))
survfit(Surv(ttr_time, ttr_status) ~ surgery, data = all.astros )

summary(coxph(Surv(ttr_time, ttr_status) ~ radiotherapy, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ radiotherapy, data = all.astros))
survfit(Surv(ttr_time, ttr_status) ~ radiotherapy, data = all.astros )

summary(coxph(Surv(ttr_time, ttr_status) ~ alkylating_chemotherapy, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ alkylating_chemotherapy, data = all.astros))
survfit(Surv(ttr_time, ttr_status) ~ alkylating_chemotherapy, data = all.astros )

summary(coxph(Surv(ttr_time, ttr_status) ~ grade, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ grade, data = all.astros))
survfit(Surv(ttr_time, ttr_status) ~ grade, data = all.astros )

#################################################################################################################

# Univariable time-to-recurrence analysis of patients with IDH-mutant oligodendrogliomas
all.oligos <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Oligodendrogliomas")

summary(coxph(Surv(ttr_time, ttr_status) ~ age_at_diagnosis, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ age_at_diagnosis, data = all.oligos))
survfit(Surv(ttr_time, ttr_status) ~ age_at_diagnosis, data = all.oligos )

summary(coxph(Surv(ttr_time, ttr_status) ~ sex, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ sex, data = all.oligos))
survfit(Surv(ttr_time, ttr_status) ~ sex, data = all.oligos )

summary(coxph(Surv(ttr_time, ttr_status) ~ surgery, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ surgery, data = all.oligos))
survfit(Surv(ttr_time, ttr_status) ~ surgery, data = all.oligos )

summary(coxph(Surv(ttr_time, ttr_status) ~ radiotherapy, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ radiotherapy, data = all.oligos))
survfit(Surv(ttr_time, ttr_status) ~ radiotherapy, data = all.oligos )

summary(coxph(Surv(ttr_time, ttr_status) ~ alkylating_chemotherapy, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ alkylating_chemotherapy, data = all.oligos))
survfit(Surv(ttr_time, ttr_status) ~ alkylating_chemotherapy, data = all.oligos )

summary(coxph(Surv(ttr_time, ttr_status) ~ grade, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ grade, data = all.oligos))
survfit(Surv(ttr_time, ttr_status) ~ grade, data = all.oligos )

#################################################################################################################