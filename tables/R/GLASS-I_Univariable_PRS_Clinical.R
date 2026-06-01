# GLASS-I: Univariable Post-Recurrence Survival Clinical Parameters
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
         surgery =                        if_else(timing.surgery == timing.recurrent,
                                                  case_when(surgery_type == "Biopsy" & (surgery_extent_of_resection == "Biopsy"|is.na(surgery_extent_of_resection)) ~"Biopsy",
                                                            surgery_type == "Craniotomy" & (surgery_extent_of_resection == "Subtotal"|surgery_extent_of_resection == "Total"|is.na(surgery_extent_of_resection)) ~ "Resection",
                                                            surgery_type == "Craniotomy" & surgery_extent_of_resection == "Biopsy"  ~ "Biopsy",
                                                            is.na(surgery_type) & surgery_extent_of_resection == "Biopsy" ~ "Biopsy",
                                                            is.na(surgery_type) & (surgery_extent_of_resection == "Subtotal"|surgery_extent_of_resection == "Total") ~ "Resection"),
                                                  NA),
         radiotherapy =                   if_else(timing.surgery == timing.recurrent,
                                                  case_when(treatment_radiotherapy == T ~ "Treated",
                                                            treatment_radiotherapy == F ~"Non-treated"),
                                                  NA),
         alkylating_chemotherapy=         if_else(timing.surgery == timing.recurrent,
                                                  case_when(treatment_alkylating_agent == T  &  (treatment_concurrent_tmz == T|treatment_concurrent_tmz == F|is.na(treatment_concurrent_tmz))~ "Treated",
                                                            treatment_concurrent_tmz == T  &  (treatment_alkylating_agent == T |treatment_alkylating_agent == F|is.na(treatment_alkylating_agent))~ "Treated",
                                                            treatment_alkylating_agent == F  &  (treatment_concurrent_tmz == F |  is.na(treatment_concurrent_tmz))~ "Non-treated",
                                                            is.na(treatment_alkylating_agent)  &  treatment_concurrent_tmz == F~ "Non-treated"),
                                                  NA),
         glioma.type =                    if_else(idh_codel_subtype != "IDHmut-codel", "IDH-mutant Astrocytomas", "IDH-mutant Oligodendrogliomas"),
         CDKN2AB =                        if_else(timing.surgery == timing.recurrent,
                                                  if_else(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "R", "present", 
                                                          if_else(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "S", "present", "absent")),
                                                  NA),
         Hypermutant =                    if_else(HM_group == "HM", "present", "absent"),
         grade =                          if_else(timing.surgery == timing.recurrent,
                                                  case_when(CDKN2AB == "present"& idh_codel_subtype == "IDHmut-noncodel" & (grade == "II"|grade == "III"|grade == "IV"|is.na(grade)) ~ "4",
                                                            (CDKN2AB == "absent"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-noncodel" & grade == "II" ~ "2",
                                                            (CDKN2AB == "absent"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-noncodel" & grade == "III" ~ "3",
                                                            (CDKN2AB == "absent"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-noncodel" & grade == "IV" ~ "4",
                                                            (CDKN2AB == "absent"|CDKN2AB == "present"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-codel" & grade == "II" ~ "2",
                                                            (CDKN2AB == "absent"|CDKN2AB == "present"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-codel" & grade == "III" ~ "3"),
                                                  NA),
         prs_status =                       if_else(timing.surgery == timing.recurrent,
                                                    if_else(case_vital_status == "dead", 1, 0),
                                                    NA),
         prs_time =                         if_else(timing.surgery == timing.recurrent,
                                                    ((case_overall_survival_mo)/12) - ((surgical_interval_mo)/12),
                                                    NA)
  )%>%
  ungroup()%>%
  filter(prs_time != 0)

patient.lvl.data <- full.data%>%
  group_by(case_barcode)%>%
  arrange(surgery_number)%>%
  summarise(age_at_diagnosis =           age_at_diagnosis[1],
            glioma.type =                glioma.type[1],
            sex =                        sex[1],
            Hypermutant =                case_when(any(Hypermutant=="absent") ~ "absent", 
                                                   any(Hypermutant=="present") ~"present"),
            HM_group =                   case_when(any(HM_group=="HM") ~ "HM", 
                                                   any(HM_group=="NHM") ~"NHM"),
            prs_status =      case_when(any(prs_status==1) ~ 1, 
                                        any(prs_status==0) ~0),
            prs_time =        prs_time[!is.na(prs_time)][1],
            surgery =                    case_when(any(surgery=="Biopsy") ~ "Biopsy", 
                                                   any(surgery=="Resection") ~"Resection"),
            radiotherapy =               case_when(any(radiotherapy=="Treated") ~ "Treated", 
                                                   any(radiotherapy=="Non-treated") ~"Non-treated"),
            alkylating_chemotherapy =    case_when(any(alkylating_chemotherapy=="Treated") ~ "Treated", 
                                                   any(alkylating_chemotherapy=="Non-treated") ~"Non-treated"),
            grade =                      case_when(any(grade=="2") ~ "2", 
                                                   any(grade=="3") ~"3",
                                                   any(grade=="4") ~"4"),
  )%>%
  ungroup()%>%
  filter(!duplicated(case_barcode))

#### IDH-mutant oligodendrogliomas ####
# All tumors
all.oligos <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Oligodendrogliomas")

summary(coxph(Surv(prs_time, prs_status) ~ age_at_diagnosis, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ age_at_diagnosis, data = all.oligos))
survfit(Surv(prs_time, prs_status) ~ age_at_diagnosis, data = all.oligos )

summary(coxph(Surv(prs_time, prs_status) ~ sex, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ sex, data = all.oligos)) ## PH assumption violated
survfit(Surv(prs_time, prs_status) ~ sex, data = all.oligos )

# summary(coxph(Surv(prs_time, prs_status) ~ surgery, data = all.oligos ), conf.int=0.9) ## only resections

summary(coxph(Surv(prs_time, prs_status) ~ radiotherapy, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ radiotherapy, data = all.oligos))
survfit(Surv(prs_time, prs_status) ~ radiotherapy, data = all.oligos )

summary(coxph(Surv(prs_time, prs_status) ~ alkylating_chemotherapy, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ alkylating_chemotherapy, data = all.oligos))
survfit(Surv(prs_time, prs_status) ~ alkylating_chemotherapy, data = all.oligos )

summary(coxph(Surv(prs_time, prs_status) ~ grade, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ grade, data = all.oligos))
survfit(Surv(prs_time, prs_status) ~ grade, data = all.oligos )

summary(coxph(Surv(prs_time, prs_status) ~ Hypermutant, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ Hypermutant, data = all.oligos))
survfit(Surv(prs_time, prs_status) ~ Hypermutant, data = all.oligos )

# Non-hypermutant tumors
nhm.oligos <- all.oligos %>% 
  filter(HM_group == "NHM")

summary(coxph(Surv(prs_time, prs_status) ~ age_at_diagnosis, data = nhm.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ age_at_diagnosis, data = nhm.oligos))
survfit(Surv(prs_time, prs_status) ~ age_at_diagnosis, data = nhm.oligos )

summary(coxph(Surv(prs_time, prs_status) ~ sex, data = nhm.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ sex, data = nhm.oligos)) ## PH assumption violated
survfit(Surv(prs_time, prs_status) ~ sex, data = nhm.oligos )

# summary(coxph(Surv(prs_time, prs_status) ~ surgery, data = nhm.oligos ), conf.int=0.9) ## only resections

summary(coxph(Surv(prs_time, prs_status) ~ radiotherapy, data = nhm.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ radiotherapy, data = nhm.oligos))
survfit(Surv(prs_time, prs_status) ~ radiotherapy, data = nhm.oligos )

summary(coxph(Surv(prs_time, prs_status) ~ alkylating_chemotherapy, data = nhm.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ alkylating_chemotherapy, data = nhm.oligos))
survfit(Surv(prs_time, prs_status) ~ alkylating_chemotherapy, data = nhm.oligos )

summary(coxph(Surv(prs_time, prs_status) ~ grade, data = nhm.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ grade, data = nhm.oligos))
survfit(Surv(prs_time, prs_status) ~ grade, data = nhm.oligos )

# Hypermutant tumors
hm.oligos <- all.oligos %>% 
  filter(HM_group == "HM")

summary(coxph(Surv(prs_time, prs_status) ~ age_at_diagnosis, data = hm.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ age_at_diagnosis, data = hm.oligos))
survfit(Surv(prs_time, prs_status) ~ age_at_diagnosis, data = hm.oligos )

summary(coxph(Surv(prs_time, prs_status) ~ sex, data = hm.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ sex, data = hm.oligos))
survfit(Surv(prs_time, prs_status) ~ sex, data = hm.oligos )

# summary(coxph(Surv(prs_time, prs_status) ~ surgery, data = hm.oligos ), conf.int=0.9) ## only resections

summary(coxph(Surv(prs_time, prs_status) ~ radiotherapy, data = hm.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ radiotherapy, data = hm.oligos))
survfit(Surv(prs_time, prs_status) ~ radiotherapy, data = hm.oligos )

summary(coxph(Surv(prs_time, prs_status) ~ alkylating_chemotherapy, data = hm.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ alkylating_chemotherapy, data = hm.oligos))
survfit(Surv(prs_time, prs_status) ~ alkylating_chemotherapy, data = hm.oligos )

summary(coxph(Surv(prs_time, prs_status) ~ grade, data = hm.oligos ), conf.int=0.9) 
cox.zph(coxph(Surv(prs_time, prs_status) ~ grade, data = hm.oligos))
survfit(Surv(prs_time, prs_status) ~ grade, data = hm.oligos )

#### IDH-mutant astrocytomas ####
# All tumors
all.astros <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Astrocytomas")

summary(coxph(Surv(prs_time, prs_status) ~ age_at_diagnosis, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ age_at_diagnosis, data = all.astros))
survfit(Surv(prs_time, prs_status) ~ age_at_diagnosis, data = all.astros )

summary(coxph(Surv(prs_time, prs_status) ~ sex, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ sex, data = all.astros))
survfit(Surv(prs_time, prs_status) ~ sex, data = all.astros )

# summary(coxph(Surv(prs_time, prs_status) ~ surgery, data = all.astros ), conf.int=0.9) ## only resections

summary(coxph(Surv(prs_time, prs_status) ~ radiotherapy, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ radiotherapy, data = all.astros))
survfit(Surv(prs_time, prs_status) ~ radiotherapy, data = all.astros )

summary(coxph(Surv(prs_time, prs_status) ~ alkylating_chemotherapy, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ alkylating_chemotherapy, data = all.astros))
survfit(Surv(prs_time, prs_status) ~ alkylating_chemotherapy, data = all.astros )

summary(coxph(Surv(prs_time, prs_status) ~ grade, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ grade, data = all.astros)) ## PH assumption violated
survfit(Surv(prs_time, prs_status) ~ grade, data = all.astros )

summary(coxph(Surv(prs_time, prs_status) ~ Hypermutant, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ Hypermutant, data = all.astros))
survfit(Surv(prs_time, prs_status) ~ Hypermutant, data = all.astros )

# Non-hypermutant tumors
nhm.astros <- all.astros %>% 
  filter(HM_group == "NHM")

summary(coxph(Surv(prs_time, prs_status) ~ age_at_diagnosis, data = nhm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ age_at_diagnosis, data = nhm.astros))
survfit(Surv(prs_time, prs_status) ~ age_at_diagnosis, data = nhm.astros )

summary(coxph(Surv(prs_time, prs_status) ~ sex, data = nhm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ sex, data = nhm.astros))
survfit(Surv(prs_time, prs_status) ~ sex, data = nhm.astros )

# summary(coxph(Surv(prs_time, prs_status) ~ surgery, data = nhm.astros ), conf.int=0.9) ## only resections

summary(coxph(Surv(prs_time, prs_status) ~ radiotherapy, data = nhm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ radiotherapy, data = nhm.astros))
survfit(Surv(prs_time, prs_status) ~ radiotherapy, data = nhm.astros )

summary(coxph(Surv(prs_time, prs_status) ~ alkylating_chemotherapy, data = nhm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ alkylating_chemotherapy, data = nhm.astros))
survfit(Surv(prs_time, prs_status) ~ alkylating_chemotherapy, data = nhm.astros )

summary(coxph(Surv(prs_time, prs_status) ~ grade, data = nhm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ grade, data = nhm.astros)) ## PH assumption violated
survfit(Surv(prs_time, prs_status) ~ grade, data = nhm.astros )

# Hypermutant tumors
hm.astros <- all.astros %>% 
  filter(HM_group == "HM")

summary(coxph(Surv(prs_time, prs_status) ~ age_at_diagnosis, data = hm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ age_at_diagnosis, data = hm.astros))
survfit(Surv(prs_time, prs_status) ~ age_at_diagnosis, data = hm.astros )

summary(coxph(Surv(prs_time, prs_status) ~ sex, data = hm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ sex, data = hm.astros))
survfit(Surv(prs_time, prs_status) ~ sex, data = hm.astros )

# summary(coxph(Surv(prs_time, prs_status) ~ surgery, data = hm.astros ), conf.int=0.9) ## only resections

summary(coxph(Surv(prs_time, prs_status) ~ radiotherapy, data = hm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ radiotherapy, data = hm.astros))
survfit(Surv(prs_time, prs_status) ~ radiotherapy, data = hm.astros )

summary(coxph(Surv(prs_time, prs_status) ~ alkylating_chemotherapy, data = hm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ alkylating_chemotherapy, data = hm.astros))
survfit(Surv(prs_time, prs_status) ~ alkylating_chemotherapy, data = hm.astros )

summary(coxph(Surv(prs_time, prs_status) ~ grade, data = hm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(prs_time, prs_status) ~ grade, data = hm.astros))
survfit(Surv(prs_time, prs_status) ~ grade, data = hm.astros )