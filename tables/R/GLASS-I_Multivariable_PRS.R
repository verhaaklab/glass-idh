# GLASS-I: Multivariable Post-Recurrence Survival
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
         MSH6 =                            if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "MSH6" & driver_status == "mutant" & driver_change == "R", "present", 
                                                           if_else(gene_symbol == "MSH6" & driver_status == "mutant" & driver_change == "S", "present", "absent")),
                                                   NA),
         NOTCH1 =                          if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "NOTCH1" & driver_status == "mutant" & driver_change == "R", "present", 
                                                           if_else(gene_symbol == "NOTCH1" & driver_status == "mutant" & driver_change == "S", "present", "absent")),
                                                   NA),
         PI3K =                            if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "PI3K" & driver_status == "mutant" & driver_change == "R", "present", 
                                                           if_else(gene_symbol == "PI3K" & driver_status == "mutant" & driver_change == "S", "present", "absent")),
                                                   NA),
         CDKN2AB =                         if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "R", "present", 
                                                           if_else(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "S", "present", "absent")),
                                                   NA),
         Any_Amplification =               if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "Any Amplification" & driver_status == "HLAMP" & driver_change == "R", "present", 
                                                           if_else(gene_symbol == "Any Amplification" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                   NA),
         Hypermutant =                     if_else(HM_group == "HM", "present", "absent"),
         grade =                           if_else(timing.surgery == timing.recurrent,
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
  ungroup()

patient.lvl.data <- full.data%>%
  group_by(case_barcode)%>%
  arrange(surgery_number)%>%
  summarise(MSH6 =               case_when(any(MSH6=="present") ~ "present", 
                                       any(MSH6=="absent") ~"absent"),
            NOTCH1 =             case_when(any(NOTCH1=="present") ~ "present", 
                                       any(NOTCH1=="absent") ~"absent"),
            PI3K =               case_when(any(PI3K=="present") ~ "present", 
                                       any(PI3K=="absent") ~"absent"),
            CDKN2AB =            case_when(any(CDKN2AB=="present") ~ "present", 
                                       any(CDKN2AB=="absent") ~"absent"),
            Any_Amplification =  case_when(any(Any_Amplification=="present") ~ "present", 
                                       any(Any_Amplification=="absent") ~"absent"),
            age_at_diagnosis =   age_at_diagnosis[1],
            glioma.type =        glioma.type[1],
            Hypermutant =        Hypermutant[1],
            HM_group =           case_when(any(HM_group=="HM") ~ "HM", 
                                       any(HM_group=="NHM") ~"NHM"),
            prs_status =         case_when(any(prs_status==1) ~ 1, 
                                        any(prs_status==0) ~0),
            prs_time =           prs_time[!is.na(prs_time)][1],
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
  filter(glioma.type == "IDH-mutant Oligodendrogliomas") %>% 
  filter(!is.na(age_at_diagnosis) & !is.na(grade) & !is.na(MSH6) &!is.na(PI3K) & !is.na(Hypermutant) ) 

coxph.oligo <- coxph(Surv(prs_time, prs_status) ~ age_at_diagnosis + grade + MSH6 + PI3K + Hypermutant , data = all.oligos )
summary(coxph.oligo)
cox.zph(coxph.oligo)
concordance(coxph.oligo)$concordance

step(coxph.oligo, direction = c( "backward"), steps = 1000)
step.coxph.oligo <- coxph(Surv(prs_time, prs_status) ~ age_at_diagnosis +  grade + Hypermutant, data = all.oligos )
summary(step.coxph.oligo)
cox.zph(step.coxph.oligo)
concordance(step.coxph.oligo)$concordance

# Non-Hypermutant tumors
nhm.oligos <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Oligodendrogliomas") %>% 
  filter(HM_group == "NHM")%>% 
  filter(!is.na(age_at_diagnosis) & !is.na(grade)) 

coxph.oligo.nhm <- coxph(Surv(prs_time, prs_status) ~ age_at_diagnosis + grade , data = nhm.oligos )  
summary(coxph.oligo.nhm)
cox.zph(coxph.oligo.nhm)
concordance(coxph.oligo.nhm)$concordance

step(coxph.oligo.nhm, direction = c( "backward"), steps = 1000)
step.coxph.oligo.nhm <- coxph(Surv(prs_time, prs_status) ~ age_at_diagnosis +  grade, data = nhm.oligos )
summary(step.coxph.oligo.nhm)
cox.zph(step.coxph.oligo.nhm)
concordance(step.coxph.oligo.nhm)$concordance

#### IDH-mutant astrocytomas ####
# All tumors
all.astros <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Astrocytomas")%>% 
  filter(!is.na(radiotherapy) &  !is.na(alkylating_chemotherapy) &!is.na(MSH6) & !is.na(NOTCH1)& !is.na(PI3K) & !is.na(Hypermutant)& !is.na(CDKN2AB) & !is.na(Any_Amplification))

coxph.astro<- coxph(Surv(prs_time, prs_status) ~ radiotherapy + alkylating_chemotherapy + MSH6 +NOTCH1 + PI3K + Hypermutant + CDKN2AB  + Any_Amplification, data = all.astros ) 
summary(coxph.astro)
cox.zph(coxph.astro)
concordance(coxph.astro)$concordance

step(coxph.astro, direction = c( "backward"), steps = 1000)
step.coxph.astro <- coxph(Surv(prs_time, prs_status) ~ radiotherapy + alkylating_chemotherapy + PI3K  + Hypermutant + Any_Amplification , data = all.astros )
summary(step.coxph.astro)
cox.zph(step.coxph.astro)
concordance(step.coxph.astro)$concordance

# Non-Hypermutant tumors
nhm.astros <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Astrocytomas")%>% 
  filter(HM_group == "NHM")%>% 
  filter(!is.na(radiotherapy) & !is.na(alkylating_chemotherapy) & !is.na(NOTCH1)& !is.na(CDKN2AB)& !is.na(Any_Amplification))

coxph.astro.nhm<- coxph(Surv(prs_time, prs_status) ~ radiotherapy + alkylating_chemotherapy + NOTCH1 + CDKN2AB + Any_Amplification , data = nhm.astros )
summary(coxph.astro.nhm)
cox.zph(coxph.astro.nhm)
concordance(coxph.astro.nhm)$concordance

step(coxph.astro.nhm, direction = c( "backward"), steps = 1000)
step.coxph.astro.nhm <- coxph(Surv(prs_time, prs_status) ~ radiotherapy + alkylating_chemotherapy + Any_Amplification, data = nhm.astros )
summary(step.coxph.astro.nhm)
cox.zph(step.coxph.astro.nhm)
concordance(step.coxph.astro.nhm)$concordance

# Hypermutant tumors
hm.astros <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Astrocytomas")%>% 
  filter(HM_group == "HM")%>% 
  filter(!is.na(CDKN2AB))

coxph.astro.hm<- coxph(Surv(prs_time, prs_status) ~  CDKN2AB, data = hm.astros )
summary(coxph.astro.hm)
cox.zph(coxph.astro.hm)
concordance(coxph.astro.hm)$concordance

step(coxph.astro.hm, direction = c( "backward"), steps = 1000)
step.coxph.astro.hm <- coxph(Surv(prs_time, prs_status) ~ CDKN2AB, data = hm.astros )
summary(step.coxph.astro.hm)
cox.zph(step.coxph.astro.hm)
concordance(step.coxph.astro.hm)$concordance