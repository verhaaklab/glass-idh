#################################################################################################################
# Content: Multivariable post-recurrence survival stepwise-backward Cox PH models
# Author: Mircea Tesileanu
# Date: 2025.07.03
# R-version: 4.3.2
# Tables: Extended Data Table 6
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
         sex =                            if_else(case_sex == "female", "Female", "Male"),
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
         PIK3CA =                          if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "PIK3CA" & driver_status == "mutant" & driver_change == "R", "present", 
                                                           if_else(gene_symbol == "PIK3CA" & driver_status == "mutant" & driver_change == "S", "present", "absent")),
                                                   NA),
         CDKN2AB =                         if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "R", "present", 
                                                           if_else(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "S", "present", "absent")),
                                                   NA),
         ATRX =                            if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "ATRX" & driver_status == "HLDEL" & driver_change == "R", "present", 
                                                           if_else(gene_symbol == "ATRX" & driver_status == "HLDEL" & driver_change == "S", "present", "absent")),
                                                   NA),
         CCND2 =                           if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "CCND2" & driver_status == "HLAMP" & driver_change == "R", "present", 
                                                           if_else(gene_symbol == "CCND2" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                   NA),
         CDK46 =                           if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "CDK4/CDK6" & driver_status == "HLAMP" & driver_change == "R", "present", 
                                                           if_else(gene_symbol == "CDK4/CDK6" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                   NA),
         MYCMYCN =                         if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "MYC/MYCN" & driver_status == "HLAMP" & driver_change == "R", "present", 
                                                           if_else(gene_symbol == "MYC/MYCN" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
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
         os_status =                      if_else(case_vital_status == "dead", 1, 0),
         os_time =                        (case_overall_survival_mo)/12
  )%>%
  ungroup()

patient.lvl.data <- full.data%>%
  group_by(case_barcode)%>%
  arrange(surgery_number)%>%
  summarise(MSH6 =           case_when(any(MSH6=="present") ~ "present", 
                                       any(MSH6=="absent") ~"absent"),
            NOTCH1 =         case_when(any(NOTCH1=="present") ~ "present", 
                                       any(NOTCH1=="absent") ~"absent"),
            PIK3CA =         case_when(any(PIK3CA=="present") ~ "present", 
                                       any(PIK3CA=="absent") ~"absent"),
            CDKN2AB =        case_when(any(CDKN2AB=="present") ~ "present", 
                                       any(CDKN2AB=="absent") ~"absent"),
            ATRX =           case_when(any(ATRX=="present") ~ "present", 
                                       any(ATRX=="absent") ~"absent"),
            CCND2 =          case_when(any(CCND2=="present") ~ "present", 
                                       any(CCND2=="absent") ~"absent"),
            CDK46 =          case_when(any(CDK46=="present") ~ "present", 
                                       any(CDK46=="absent") ~"absent"),
            MYCMYCN =        case_when(any(MYCMYCN=="present") ~ "present", 
                                       any(MYCMYCN=="absent") ~"absent"),
            age_at_diagnosis =           age_at_diagnosis[1],
            glioma.type =                glioma.type[1],
            sex =                        sex[1],
            Hypermutant =                Hypermutant[1],
            HM_group =       case_when(any(HM_group=="HM") ~ "HM", 
                                       any(HM_group=="NHM") ~"NHM"),
            os_status =                  os_status[1],
            os_time =                    os_time[1],
            alkylating_chemotherapy =    case_when(any(alkylating_chemotherapy=="Treated") ~ "Treated", 
                                                   any(alkylating_chemotherapy=="Non-treated") ~"Non-treated"),
            grade =                      case_when(any(grade=="2") ~ "2", 
                                                   any(grade=="3") ~"3",
                                                   any(grade=="4") ~"4"),
  )%>%
  ungroup()%>%
  filter(!duplicated(case_barcode))

#################################################################################################################

# Multivariable post-recurrence survival analysis of patients with IDH-mutant astrocytomas
# All tumors
all.astros <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Astrocytomas")%>% 
  filter(!is.na(alkylating_chemotherapy) & !is.na(grade) & !is.na(NOTCH1)& !is.na(Hypermutant)&
           !is.na(CDKN2AB)& !is.na(CDK46))

coxph.astro<- coxph(Surv(os_time, os_status) ~ alkylating_chemotherapy + grade +NOTCH1 + Hypermutant + 
                    CDKN2AB  + CDK46 , data = all.astros ) 
summary(coxph.astro)
cox.zph(coxph.astro)
concordance(coxph.astro)$concord

step(coxph.astro, direction = c( "backward"), steps = 1000)
step.coxph.astro <- coxph(Surv(os_time, os_status) ~ alkylating_chemotherapy + Hypermutant + CDK46 , data = all.astros )
summary(step.coxph.astro)
cox.zph(step.coxph.astro)
concordance(step.coxph.astro)$concordance

# Hypermutant tumors
hm.astros <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Astrocytomas")%>% 
  filter(HM_group == "HM")

coxph.astro.hm<- coxph(Surv(os_time, os_status) ~ ATRX  + MYCMYCN, data = hm.astros )
summary(coxph.astro.hm)
cox.zph(coxph.astro.hm)
concordance(coxph.astro.hm)$concordance

step(coxph.astro.hm, direction = c( "backward"), steps = 1000)
step.coxph.astro.hm <- coxph(Surv(os_time, os_status) ~ ATRX, data = hm.astros )
summary(step.coxph.astro.hm)
cox.zph(step.coxph.astro.hm)
concordance(step.coxph.astro.hm)$concordance

# Non-Hypermutant tumors
nhm.astros <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Astrocytomas")%>% 
  filter(HM_group == "NHM")%>% 
  filter(!is.na(age_at_diagnosis) & !is.na(sex) & !is.na(grade)& !is.na(CCND2))

coxph.astro.nhm<- coxph(Surv(os_time, os_status) ~ age_at_diagnosis + sex + grade + CCND2 , data = nhm.astros )
summary(coxph.astro.nhm)
cox.zph(coxph.astro.nhm)
concordance(coxph.astro.nhm)$concordance

step(coxph.astro.nhm, direction = c( "backward"), steps = 1000)
step.coxph.astro.nhm <- coxph(Surv(os_time, os_status) ~ CCND2, data = nhm.astros )
summary(step.coxph.astro.nhm)
cox.zph(step.coxph.astro.nhm)
concordance(step.coxph.astro.nhm)$concordance

#################################################################################################################

# Multivariable post-recurrence survival analysis of patients with IDH-mutant oligodendrogliomas
# All tumors
all.oligos <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Oligodendrogliomas") %>% 
  filter(!is.na(age_at_diagnosis) & !is.na(alkylating_chemotherapy) & !is.na(grade) & 
           !is.na(MSH6) & !is.na(NOTCH1) &!is.na(PIK3CA) & !is.na(Hypermutant) ) 

coxph.oligo <- coxph(Surv(os_time, os_status) ~ age_at_diagnosis + alkylating_chemotherapy  + grade + 
                       MSH6 + NOTCH1 + PIK3CA + Hypermutant , data = all.oligos )
summary(coxph.oligo)
cox.zph(coxph.oligo)
concordance(coxph.oligo)$concordance

step(coxph.oligo, direction = c( "backward"), steps = 1000)
step.coxph.oligo <- coxph(Surv(os_time, os_status) ~ age_at_diagnosis +  alkylating_chemotherapy + Hypermutant, data = all.oligos )
summary(step.coxph.oligo)
cox.zph(step.coxph.oligo)
concordance(step.coxph.oligo)$concordance

# Non-Hypermutant tumors
nhm.oligos <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Oligodendrogliomas") %>% 
  filter(HM_group == "NHM")%>% 
  filter(!is.na(age_at_diagnosis)  & !is.na(alkylating_chemotherapy)) 

coxph.oligo.nhm <- coxph(Surv(os_time, os_status) ~ age_at_diagnosis + alkylating_chemotherapy , data = nhm.oligos )  
summary(coxph.oligo.nhm)
cox.zph(coxph.oligo.nhm)
concordance(coxph.oligo.nhm)$concordance

step(coxph.oligo.nhm, direction = c( "backward"), steps = 1000)
step.coxph.oligo.nhm <- coxph(Surv(os_time, os_status) ~ age_at_diagnosis +  alkylating_chemotherapy, data = nhm.oligos )
summary(step.coxph.oligo.nhm)
cox.zph(step.coxph.oligo.nhm)
concordance(step.coxph.oligo.nhm)$concordance

#################################################################################################################