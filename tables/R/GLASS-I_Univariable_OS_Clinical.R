# GLASS-I: Univariable Overall Survival Clinical Parameters
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(RPostgres)
library(tidyverse)
library(survival)
library(survminer)

# load data
molecular.data <- read.delim("~/driver_changes.txt")%>%                       ## created with GLASS-I_Drivers_Table.R
  select(-idh_codel_subtype) 

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")             ## data available on Synapse (Synapse ID: syn17038081)

#load surgeries and filter missing survival data samples
surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
  filter(idh_codel_subtype != "IDHwt")%>%
  filter(!is.na(surgical_interval_mo))%>%
  mutate(timing.surgery = substr(sample_barcode, 14,15)) 

gold.set <- dbGetQuery(con, "SELECT * FROM analysis.gold_set")%>%
  mutate(timing.recurrent = substr(tumor_pair_barcode, 20,21),
         timing.primary = substr(tumor_pair_barcode, 14,15) )%>%
  semi_join(surgeries, by = "case_barcode")

survival.data <- dbGetQuery(con, "SELECT * FROM clinical.cases")[,c(2,5:8)]%>%
  semi_join(surgeries, by="case_barcode")%>%
  filter(!is.na(case_overall_survival_mo))

dbDisconnect(con)

full.data <- survival.data%>%
  left_join(surgeries, by ="case_barcode")%>%
  left_join(molecular.data, by ="case_barcode")%>%
  left_join(gold.set, by ="case_barcode")%>%
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
         CDKN2AB =                        if_else(surgery_number == 1 & timing.primary == timing.surgery,
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
         os_status =                      if_else(case_vital_status == "dead", 1, 0),
         os_time =                        (case_overall_survival_mo)/12
  )%>%
  ungroup()

patient.lvl.data <- full.data%>%
  group_by(case_barcode)%>%
  summarise(age_at_diagnosis =           case_when(any(age_at_diagnosis=="40+") ~ "40+", 
                                                   any(age_at_diagnosis=="40 or younger") ~"40 or younger"),
            sex =                        case_when(any(sex=="Female") ~ "Female", 
                                                   any(sex=="Male") ~"Male"),
            os_status =                  os_status[1],
            os_time =                    os_time[1],
            surgery =                    case_when(any(surgery=="Biopsy") ~ "Biopsy", 
                                                   any(surgery=="Resection") ~"Resection"),
            radiotherapy =               case_when(any(radiotherapy=="Treated") ~ "Treated", 
                                                   any(radiotherapy=="Non-treated") ~"Non-treated"),
            alkylating_chemotherapy =    case_when(any(alkylating_chemotherapy=="Treated") ~ "Treated", 
                                                   any(alkylating_chemotherapy=="Non-treated") ~"Non-treated"),
            grade =                      case_when(any(grade=="2") ~ "2", 
                                                   any(grade=="3") ~"3",
                                                   any(grade=="4") ~"4"),
            glioma.type =                case_when(any(glioma.type=="IDH-mutant Oligodendrogliomas") ~ "IDH-mutant Oligodendrogliomas", 
                                                   any(glioma.type=="IDH-mutant Astrocytomas") ~"IDH-mutant Astrocytomas")
  )%>%
  ungroup()%>%
  filter(!duplicated(case_barcode))

#### IDH-mutant oligodendrogliomas ####
# All tumors
all.oligos <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Oligodendrogliomas")

summary(coxph(Surv(os_time, os_status) ~ age_at_diagnosis, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ age_at_diagnosis, data = all.oligos))
survfit(Surv(os_time, os_status) ~ age_at_diagnosis, data = all.oligos )

summary(coxph(Surv(os_time, os_status) ~ sex, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ sex, data = all.oligos)) 
survfit(Surv(os_time, os_status) ~ sex, data = all.oligos )

summary(coxph(Surv(os_time, os_status) ~ surgery, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ surgery, data = all.oligos))
survfit(Surv(os_time, os_status) ~ surgery, data = all.oligos )

summary(coxph(Surv(os_time, os_status) ~ radiotherapy, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ radiotherapy, data = all.oligos))
survfit(Surv(os_time, os_status) ~ radiotherapy, data = all.oligos )

summary(coxph(Surv(os_time, os_status) ~ alkylating_chemotherapy, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ alkylating_chemotherapy, data = all.oligos))
survfit(Surv(os_time, os_status) ~ alkylating_chemotherapy, data = all.oligos )

summary(coxph(Surv(os_time, os_status) ~ grade, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ grade, data = all.oligos))
survfit(Surv(os_time, os_status) ~ grade, data = all.oligos )

#### IDH-mutant astrocytomas ####
# All tumors
all.astros <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Astrocytomas")

summary(coxph(Surv(os_time, os_status) ~ age_at_diagnosis, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ age_at_diagnosis, data = all.astros))
survfit(Surv(os_time, os_status) ~ age_at_diagnosis, data = all.astros )

summary(coxph(Surv(os_time, os_status) ~ sex, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ sex, data = all.astros))
survfit(Surv(os_time, os_status) ~ sex, data = all.astros )

summary(coxph(Surv(os_time, os_status) ~ surgery, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ surgery, data = all.astros))
survfit(Surv(os_time, os_status) ~ surgery, data = all.astros )

summary(coxph(Surv(os_time, os_status) ~ radiotherapy, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ radiotherapy, data = all.astros))
survfit(Surv(os_time, os_status) ~ radiotherapy, data = all.astros )

summary(coxph(Surv(os_time, os_status) ~ alkylating_chemotherapy, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ alkylating_chemotherapy, data = all.astros))
survfit(Surv(os_time, os_status) ~ alkylating_chemotherapy, data = all.astros )

summary(coxph(Surv(os_time, os_status) ~ grade, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ grade, data = all.astros))
survfit(Surv(os_time, os_status) ~ grade, data = all.astros )