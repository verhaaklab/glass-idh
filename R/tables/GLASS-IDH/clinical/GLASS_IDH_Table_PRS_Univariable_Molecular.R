#################################################################################################################
# Content: Univariable post-recurrence survival Cox PH models of molecular variables
# Author: Mircea Tesileanu
# Date: 2025.07.03
# R-version: 4.3.2
# Tables: Extended Data Table 5
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
  mutate(MSH6 =                            if_else(timing.surgery == timing.recurrent,
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
         PIK3R1 =                          if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "PIK3R1" & driver_status == "mutant" & driver_change == "R", "present", 
                                                   if_else(gene_symbol == "PIK3R1" & driver_status == "mutant" & driver_change == "S", "present", "absent")),
                                                   NA),
         CDKN2AB =                         if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "R", "present", 
                                                   if_else(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "S", "present", "absent")),
                                                   NA),
         PTEN =                            if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "PTEN" & driver_status == "HLDEL" & driver_change == "R", "present", 
                                                   if_else(gene_symbol == "PTEN" & driver_status == "HLDEL" & driver_change == "S", "present", "absent")),
                                                   NA),
         ATRX =                            if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "ATRX" & driver_status == "HLDEL" & driver_change == "R", "present", 
                                                   if_else(gene_symbol == "ATRX" & driver_status == "HLDEL" & driver_change == "S", "present", "absent")),
                                                   NA),
         CCND2 =                           if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "CCND2" & driver_status == "HLAMP" & driver_change == "R", "present", 
                                                   if_else(gene_symbol == "CCND2" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                   NA),
         EGFR =                            if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "EGFR" & driver_status == "HLAMP" & driver_change == "R", "present", 
                                                   if_else(gene_symbol == "EGFR" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                   NA),
         PDGFRA =                          if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "PDGFRA" & driver_status == "HLAMP" & driver_change == "R", "present", 
                                                   if_else(gene_symbol == "PDGFRA" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                   NA),
         CDK46 =                           if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "CDK4/CDK6" & driver_status == "HLAMP" & driver_change == "R", "present", 
                                                   if_else(gene_symbol == "CDK4/CDK6" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                   NA),
         METPTPRZ1 =                       if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "MET/PTPRZ1" & driver_status == "HLAMP" & driver_change == "R", "present", 
                                                   if_else(gene_symbol == "MET/PTPRZ1" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                   NA),
         MYCMYCN =                         if_else(timing.surgery == timing.recurrent,
                                                   if_else(gene_symbol == "MYC/MYCN" & driver_status == "HLAMP" & driver_change == "R", "present", 
                                                   if_else(gene_symbol == "MYC/MYCN" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                   NA),
         MGMT =                            if_else(timing.surgery == timing.recurrent,
                                                   if_else(mgmt_methylation == "U", "absent", "present"),
                                                   NA),
         glioma.type =                     if_else(idh_codel_subtype != "IDHmut-codel", "IDH-mutant Astrocytomas", "IDH-mutant Oligodendrogliomas"),
         os_status =                       if_else(case_vital_status == "dead", 1, 0),
         os_time =                         (case_overall_survival_mo)/12
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
            PIK3R1 =         case_when(any(PIK3R1=="present") ~ "present", 
                                       any(PIK3R1=="absent") ~"absent"),
            CDKN2AB =        case_when(any(CDKN2AB=="present") ~ "present", 
                                       any(CDKN2AB=="absent") ~"absent"),
            PTEN =           case_when(any(PTEN=="present") ~ "present", 
                                       any(PTEN=="absent") ~"absent"),
            ATRX =           case_when(any(ATRX=="present") ~ "present", 
                                       any(ATRX=="absent") ~"absent"),
            CCND2 =          case_when(any(CCND2=="present") ~ "present", 
                                       any(CCND2=="absent") ~"absent"),
            EGFR =           case_when(any(EGFR=="present") ~ "present", 
                                       any(EGFR=="absent") ~"absent"),
            PDGFRA =         case_when(any(PDGFRA=="present") ~ "present", 
                                       any(PDGFRA=="absent") ~"absent"),
            CDK46 =          case_when(any(CDK46=="present") ~ "present", 
                                       any(CDK46=="absent") ~"absent"),
            METPTPRZ1 =      case_when(any(METPTPRZ1=="present") ~ "present", 
                                       any(METPTPRZ1=="absent") ~"absent"),
            MYCMYCN =        case_when(any(MYCMYCN=="present") ~ "present", 
                                       any(MYCMYCN=="absent") ~"absent"),
            MGMT =           case_when(any(MGMT=="present") ~ "present", 
                                       any(MGMT=="absent") ~"absent"),
            glioma.type =    glioma.type[1],
            HM_group =       case_when(any(HM_group=="HM") ~ "HM", 
                                       any(HM_group=="NHM") ~"NHM"),
            os_status =      os_status[1],
            os_time =        os_time[1]
  )%>%
  ungroup()%>%
  filter(!duplicated(case_barcode))

#################################################################################################################

# Univariable post-recurrence survival analysis of patients with IDH-mutant astrocytomas
# All tumors
all.astros <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Astrocytomas")

summary(coxph(Surv(os_time, os_status) ~ MSH6, data = all.astros ), conf.int=0.9) 
cox.zph(coxph(Surv(os_time, os_status) ~ MSH6, data = all.astros))
survfit(Surv(os_time, os_status) ~ MSH6, data = all.astros )

summary(coxph(Surv(os_time, os_status) ~ NOTCH1, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ NOTCH1, data = all.astros))
survfit(Surv(os_time, os_status) ~ NOTCH1, data = all.astros )

summary(coxph(Surv(os_time, os_status) ~ PIK3CA, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ PIK3CA, data = all.astros))
survfit(Surv(os_time, os_status) ~ PIK3CA, data = all.astros )

summary(coxph(Surv(os_time, os_status) ~ PIK3R1, data = all.astros ), conf.int=0.9) 
cox.zph(coxph(Surv(os_time, os_status) ~ PIK3R1, data = all.astros))
survfit(Surv(os_time, os_status) ~ PIK3R1, data = all.astros )

summary(coxph(Surv(os_time, os_status) ~ CDKN2AB, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ CDKN2AB, data = all.astros))
survfit(Surv(os_time, os_status) ~ CDKN2AB, data = all.astros)

# summary(coxph(Surv(os_time, os_status) ~ PTEN, data = all.astros), conf.int=0.9) ## only non-HD

summary(coxph(Surv(os_time, os_status) ~ ATRX, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ ATRX, data = all.astros))
survfit(Surv(os_time, os_status) ~ ATRX, data = all.astros)

summary(coxph(Surv(os_time, os_status) ~ CCND2, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ CCND2, data = all.astros))
survfit(Surv(os_time, os_status) ~ CCND2, data = all.astros)

summary(coxph(Surv(os_time, os_status) ~ EGFR, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ EGFR, data = all.astros))
survfit(Surv(os_time, os_status) ~ EGFR, data = all.astros)

summary(coxph(Surv(os_time, os_status) ~ PDGFRA, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ PDGFRA, data = all.astros)) 
survfit(Surv(os_time, os_status) ~ PDGFRA, data = all.astros)

summary(coxph(Surv(os_time, os_status) ~ CDK46, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ CDK46, data = all.astros))
survfit(Surv(os_time, os_status) ~ CDK46, data = all.astros)

summary(coxph(Surv(os_time, os_status) ~ METPTPRZ1, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ METPTPRZ1, data = all.astros))
survfit(Surv(os_time, os_status) ~ METPTPRZ1, data = all.astros)

summary(coxph(Surv(os_time, os_status) ~ MYCMYCN, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ MYCMYCN, data = all.astros)) ## PH assumption violated
survfit(Surv(os_time, os_status) ~ MYCMYCN, data = all.astros)

summary(coxph(Surv(os_time, os_status) ~ MGMT, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ MGMT, data = all.astros))
survfit(Surv(os_time, os_status) ~ MGMT, data = all.astros)

# Hypermutant tumors
hm.astros <- all.astros %>% 
  filter(HM_group == "HM")

summary(coxph(Surv(os_time, os_status) ~ MSH6, data = hm.astros ), conf.int=0.9) 
cox.zph(coxph(Surv(os_time, os_status) ~ MSH6, data = hm.astros))
survfit(Surv(os_time, os_status) ~ MSH6, data = hm.astros )

summary(coxph(Surv(os_time, os_status) ~ NOTCH1, data = hm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ NOTCH1, data = hm.astros))
survfit(Surv(os_time, os_status) ~ NOTCH1, data = hm.astros )

summary(coxph(Surv(os_time, os_status) ~ PIK3CA, data = hm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ PIK3CA, data = hm.astros))
survfit(Surv(os_time, os_status) ~ PIK3CA, data = hm.astros )

summary(coxph(Surv(os_time, os_status) ~ PIK3R1, data = hm.astros ), conf.int=0.9) 
cox.zph(coxph(Surv(os_time, os_status) ~ PIK3R1, data = hm.astros))
survfit(Surv(os_time, os_status) ~ PIK3R1, data = hm.astros )

summary(coxph(Surv(os_time, os_status) ~ CDKN2AB, data = hm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ CDKN2AB, data = hm.astros))
survfit(Surv(os_time, os_status) ~ CDKN2AB, data = hm.astros)

# summary(coxph(Surv(os_time, os_status) ~ PTEN, data = hm.astros ), conf.int=0.9) ## only non-HD

summary(coxph(Surv(os_time, os_status) ~ ATRX, data = hm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ ATRX, data = hm.astros))
survfit(Surv(os_time, os_status) ~ ATRX, data = hm.astros)

summary(coxph(Surv(os_time, os_status) ~ CCND2, data = hm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ CCND2, data = hm.astros))
survfit(Surv(os_time, os_status) ~ CCND2, data = hm.astros)

# summary(coxph(Surv(os_time, os_status) ~ EGFR, data = hm.astros ), conf.int=0.9) ## only non-amp

summary(coxph(Surv(os_time, os_status) ~ PDGFRA, data = hm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ PDGFRA, data = hm.astros))
survfit(Surv(os_time, os_status) ~ PDGFRA, data = hm.astros)

summary(coxph(Surv(os_time, os_status) ~ CDK46, data = hm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ CDK46, data = hm.astros))
survfit(Surv(os_time, os_status) ~ CDK46, data = hm.astros)

summary(coxph(Surv(os_time, os_status) ~ METPTPRZ1, data = hm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ METPTPRZ1, data = hm.astros))
survfit(Surv(os_time, os_status) ~ METPTPRZ1, data = hm.astros)

summary(coxph(Surv(os_time, os_status) ~ MYCMYCN, data = hm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ MYCMYCN, data = hm.astros))
survfit(Surv(os_time, os_status) ~ MYCMYCN, data = hm.astros)

summary(coxph(Surv(os_time, os_status) ~ MGMT, data = hm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ MGMT, data = hm.astros))
survfit(Surv(os_time, os_status) ~ MGMT, data = hm.astros)

# Non-Hypermutant tumors
nhm.astros <- all.astros %>% 
  filter(HM_group == "NHM")

# summary(coxph(Surv(os_time, os_status) ~ MSH6, data = nhm.astros ), conf.int=0.9) ## only wildtype

summary(coxph(Surv(os_time, os_status) ~ NOTCH1, data = nhm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ NOTCH1, data = nhm.astros)) 
survfit(Surv(os_time, os_status) ~ NOTCH1, data = nhm.astros )

summary(coxph(Surv(os_time, os_status) ~ PIK3CA, data = nhm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ PIK3CA, data = nhm.astros))
survfit(Surv(os_time, os_status) ~ PIK3CA, data = nhm.astros )

summary(coxph(Surv(os_time, os_status) ~ PIK3R1, data = nhm.astros ), conf.int=0.9) 
cox.zph(coxph(Surv(os_time, os_status) ~ PIK3R1, data = nhm.astros))
survfit(Surv(os_time, os_status) ~ PIK3R1, data = nhm.astros )

summary(coxph(Surv(os_time, os_status) ~ CDKN2AB, data = nhm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ CDKN2AB, data = nhm.astros)) ## PH assumption violated
survfit(Surv(os_time, os_status) ~ CDKN2AB, data = nhm.astros)

# summary(coxph(Surv(os_time, os_status) ~ PTEN, data = nhm.astros ), conf.int=0.9) ## only non-HD

summary(coxph(Surv(os_time, os_status) ~ ATRX, data = nhm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ ATRX, data = nhm.astros))
survfit(Surv(os_time, os_status) ~ ATRX, data = nhm.astros)

summary(coxph(Surv(os_time, os_status) ~ CCND2, data = nhm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ CCND2, data = nhm.astros)) 
survfit(Surv(os_time, os_status) ~ CCND2, data = nhm.astros)

summary(coxph(Surv(os_time, os_status) ~ EGFR, data = nhm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ EGFR, data = nhm.astros))
survfit(Surv(os_time, os_status) ~ EGFR, data = nhm.astros)

summary(coxph(Surv(os_time, os_status) ~ PDGFRA, data = nhm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ PDGFRA, data = nhm.astros))
survfit(Surv(os_time, os_status) ~ PDGFRA, data = nhm.astros)

summary(coxph(Surv(os_time, os_status) ~ CDK46, data = nhm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ CDK46, data = nhm.astros))
survfit(Surv(os_time, os_status) ~ CDK46, data = nhm.astros)

summary(coxph(Surv(os_time, os_status) ~ METPTPRZ1, data = nhm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ METPTPRZ1, data = nhm.astros))
survfit(Surv(os_time, os_status) ~ METPTPRZ1, data = nhm.astros)

summary(coxph(Surv(os_time, os_status) ~ MYCMYCN, data = nhm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ MYCMYCN, data = nhm.astros))
survfit(Surv(os_time, os_status) ~ MYCMYCN, data = nhm.astros)

summary(coxph(Surv(os_time, os_status) ~ MGMT, data = nhm.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ MGMT, data = nhm.astros))
survfit(Surv(os_time, os_status) ~ MGMT, data = nhm.astros)

#################################################################################################################

# Univariable post-recurrence survival analysis of patients with IDH-mutant oligodendrogliomas
# All tumors
all.oligos <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Oligodendrogliomas")

summary(coxph(Surv(os_time, os_status) ~ MSH6, data = all.oligos ), conf.int=0.9) 
cox.zph(coxph(Surv(os_time, os_status) ~ MSH6, data = all.oligos))
survfit(Surv(os_time, os_status) ~ MSH6, data = all.oligos )

summary(coxph(Surv(os_time, os_status) ~ NOTCH1, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ NOTCH1, data = all.oligos))
survfit(Surv(os_time, os_status) ~ NOTCH1, data = all.oligos )

summary(coxph(Surv(os_time, os_status) ~ PIK3CA, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ PIK3CA, data = all.oligos))
survfit(Surv(os_time, os_status) ~ PIK3CA, data = all.oligos )

summary(coxph(Surv(os_time, os_status) ~ PIK3R1, data = all.oligos ), conf.int=0.9) 
cox.zph(coxph(Surv(os_time, os_status) ~ PIK3R1, data = all.oligos))
survfit(Surv(os_time, os_status) ~ PIK3R1, data = all.oligos )

summary(coxph(Surv(os_time, os_status) ~ CDKN2AB, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ CDKN2AB, data = all.oligos)) ## PH assumption violated
survfit(Surv(os_time, os_status) ~ CDKN2AB, data = all.oligos)

# summary(coxph(Surv(os_time, os_status) ~ PTEN, data = all.oligos ), conf.int=0.9) ## only non-HD

# summary(coxph(Surv(os_time, os_status) ~ ATRX, data = all.oligos ), conf.int=0.9)## only non-HD

# summary(coxph(Surv(os_time, os_status) ~ CCND2, data = all.oligos ), conf.int=0.9) ## only non-amp

# summary(coxph(Surv(os_time, os_status) ~ EGFR, data = all.oligos ), conf.int=0.9) ## only non-amp

# summary(coxph(Surv(os_time, os_status) ~ PDGFRA, data = all.oligos ), conf.int=0.9) ## only non-amp

summary(coxph(Surv(os_time, os_status) ~ CDK46, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ CDK46, data = all.oligos))
survfit(Surv(os_time, os_status) ~ CDK46, data = all.oligos)

summary(coxph(Surv(os_time, os_status) ~ METPTPRZ1, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ METPTPRZ1, data = all.oligos))
survfit(Surv(os_time, os_status) ~ METPTPRZ1, data = all.oligos)

# summary(coxph(Surv(os_time, os_status) ~ MYCMYCN, data = all.oligos ), conf.int=0.9) ## only non-amp

summary(coxph(Surv(os_time, os_status) ~ MGMT, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ MGMT, data = all.oligos))
survfit(Surv(os_time, os_status) ~ MGMT, data = all.oligos)

# Hypermutant tumors
hm.oligos <- all.oligos %>% 
  filter(HM_group == "HM")

summary(coxph(Surv(os_time, os_status) ~ MSH6, data = hm.oligos ), conf.int=0.9) 
cox.zph(coxph(Surv(os_time, os_status) ~ MSH6, data = hm.oligos))
survfit(Surv(os_time, os_status) ~ MSH6, data = hm.oligos )

summary(coxph(Surv(os_time, os_status) ~ NOTCH1, data = hm.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ NOTCH1, data = hm.oligos))
survfit(Surv(os_time, os_status) ~ NOTCH1, data = hm.oligos )

summary(coxph(Surv(os_time, os_status) ~ PIK3CA, data = hm.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ PIK3CA, data = hm.oligos))
survfit(Surv(os_time, os_status) ~ PIK3CA, data = hm.oligos )

summary(coxph(Surv(os_time, os_status) ~ PIK3R1, data = hm.oligos ), conf.int=0.9) 
cox.zph(coxph(Surv(os_time, os_status) ~ PIK3R1, data = hm.oligos))
survfit(Surv(os_time, os_status) ~ PIK3R1, data = hm.oligos )

# summary(coxph(Surv(os_time, os_status) ~ CDKN2AB, data = hm.oligos ), conf.int=0.9) ## only non-HD

# summary(coxph(Surv(os_time, os_status) ~ PTEN, data = hm.oligos ), conf.int=0.9) ## only non-HD

# summary(coxph(Surv(os_time, os_status) ~ ATRX, data = hm.oligos ), conf.int=0.9) ## only non-HD

# summary(coxph(Surv(os_time, os_status) ~ CCND2, data = hm.oligos ), conf.int=0.9) ## only non-amp

# summary(coxph(Surv(os_time, os_status) ~ EGFR, data = hm.oligos ), conf.int=0.9)## only non-amp

# summary(coxph(Surv(os_time, os_status) ~ PDGFRA, data = hm.oligos ), conf.int=0.9) ## only non-amp

summary(coxph(Surv(os_time, os_status) ~ CDK46, data = hm.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ CDK46, data = hm.oligos))
survfit(Surv(os_time, os_status) ~ CDK46, data = hm.oligos)

# summary(coxph(Surv(os_time, os_status) ~ METPTPRZ1, data = hm.oligos ), conf.int=0.9)## only non-amp

# summary(coxph(Surv(os_time, os_status) ~ MYCMYCN, data = hm.oligos ), conf.int=0.9) ## only non-amp

# summary(coxph(Surv(os_time, os_status) ~ MGMT, data = hm.oligos ), conf.int=0.9) ## only methylated

# Non-Hypermutant tumors
nhm.oligos <- all.oligos %>% 
  filter(HM_group == "NHM")

# summary(coxph(Surv(os_time, os_status) ~ MSH6, data = nhm.oligos ), conf.int=0.9) ## only wildtype

summary(coxph(Surv(os_time, os_status) ~ NOTCH1, data = nhm.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ NOTCH1, data = nhm.oligos))
survfit(Surv(os_time, os_status) ~ NOTCH1, data = nhm.oligos )

summary(coxph(Surv(os_time, os_status) ~ PIK3CA, data = nhm.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ PIK3CA, data = nhm.oligos))
survfit(Surv(os_time, os_status) ~ PIK3CA, data = nhm.oligos )

summary(coxph(Surv(os_time, os_status) ~ PIK3R1, data = nhm.oligos ), conf.int=0.9) 
cox.zph(coxph(Surv(os_time, os_status) ~ PIK3R1, data = nhm.oligos))
survfit(Surv(os_time, os_status) ~ PIK3R1, data = nhm.oligos )

summary(coxph(Surv(os_time, os_status) ~ CDKN2AB, data = nhm.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ CDKN2AB, data = nhm.oligos)) ## PH assumption violated
survfit(Surv(os_time, os_status) ~ CDKN2AB, data = nhm.oligos)

# summary(coxph(Surv(os_time, os_status) ~ PTEN, data = nhm.oligos ), conf.int=0.9) ## only non-HD

# summary(coxph(Surv(os_time, os_status) ~ ATRX, data = nhm.oligos ), conf.int=0.9)## only non-HD

# summary(coxph(Surv(os_time, os_status) ~ CCND2, data = nhm.oligos ), conf.int=0.9)## only non-amp

# summary(coxph(Surv(os_time, os_status) ~ EGFR, data = nhm.oligos ), conf.int=0.9)## only non-amp

# summary(coxph(Surv(os_time, os_status) ~ PDGFRA, data = nhm.oligos ), conf.int=0.9)## only non-amp

summary(coxph(Surv(os_time, os_status) ~ CDK46, data = nhm.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ CDK46, data = nhm.oligos))
survfit(Surv(os_time, os_status) ~ CDK46, data = nhm.oligos)

summary(coxph(Surv(os_time, os_status) ~ METPTPRZ1, data = nhm.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ METPTPRZ1, data = nhm.oligos))
survfit(Surv(os_time, os_status) ~ METPTPRZ1, data = nhm.oligos)

# summary(coxph(Surv(os_time, os_status) ~ MYCMYCN, data = nhm.oligos ), conf.int=0.9)## only non-amp

summary(coxph(Surv(os_time, os_status) ~ MGMT, data = nhm.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ MGMT, data = nhm.oligos))
survfit(Surv(os_time, os_status) ~ MGMT, data = nhm.oligos)

#################################################################################################################