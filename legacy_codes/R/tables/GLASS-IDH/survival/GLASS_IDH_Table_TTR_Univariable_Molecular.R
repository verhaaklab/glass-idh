#################################################################################################################
# Content: Univariable time-to-recurrence Cox PH models of molecular variables
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
  mutate(MSH6 =                            if_else(timing.primary == "TP" & surgery_number == 1,
                                                   if_else(gene_symbol == "MSH6" & driver_status == "mutant" & driver_change == "P", "present", 
                                                           if_else(gene_symbol == "MSH6" & driver_status == "mutant" & driver_change == "S", "present", "absent")),
                                                   NA),
         NOTCH1 =                          if_else(timing.primary == "TP" & surgery_number == 1,
                                                   if_else(gene_symbol == "NOTCH1" & driver_status == "mutant" & driver_change == "P", "present", 
                                                           if_else(gene_symbol == "NOTCH1" & driver_status == "mutant" & driver_change == "S", "present", "absent")),
                                                   NA),
         PIK3CA =                          if_else(timing.primary == "TP" & surgery_number == 1,
                                                   if_else(gene_symbol == "PIK3CA" & driver_status == "mutant" & driver_change == "P", "present", 
                                                           if_else(gene_symbol == "PIK3CA" & driver_status == "mutant" & driver_change == "S", "present", "absent")),
                                                   NA),
         PIK3R1 =                          if_else(timing.primary == "TP" & surgery_number == 1,
                                                   if_else(gene_symbol == "PIK3R1" & driver_status == "mutant" & driver_change == "P", "present", 
                                                           if_else(gene_symbol == "PIK3R1" & driver_status == "mutant" & driver_change == "S", "present", "absent")),
                                                   NA),
         CDKN2AB =                         if_else(timing.primary == "TP" & surgery_number == 1,
                                                   if_else(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "P", "present", 
                                                           if_else(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "S", "present", "absent")),
                                                   NA),
         PTEN =                            if_else(timing.primary == "TP" & surgery_number == 1,
                                                   if_else(gene_symbol == "PTEN" & driver_status == "HLDEL" & driver_change == "P", "present", 
                                                           if_else(gene_symbol == "PTEN" & driver_status == "HLDEL" & driver_change == "S", "present", "absent")),
                                                   NA),
         ATRX =                            if_else(timing.primary == "TP" & surgery_number == 1,
                                                   if_else(gene_symbol == "ATRX" & driver_status == "HLDEL" & driver_change == "P", "present", 
                                                           if_else(gene_symbol == "ATRX" & driver_status == "HLDEL" & driver_change == "S", "present", "absent")),
                                                   NA),
         CCND2 =                           if_else(timing.primary == "TP" & surgery_number == 1,
                                                   if_else(gene_symbol == "CCND2" & driver_status == "HLAMP" & driver_change == "P", "present", 
                                                           if_else(gene_symbol == "CCND2" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                   NA),
         EGFR =                            if_else(timing.primary == "TP" & surgery_number == 1,
                                                   if_else(gene_symbol == "EGFR" & driver_status == "HLAMP" & driver_change == "P", "present", 
                                                           if_else(gene_symbol == "EGFR" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                   NA),
         PDGFRA =                          if_else(timing.primary == "TP" & surgery_number == 1,
                                                   if_else(gene_symbol == "PDGFRA" & driver_status == "HLAMP" & driver_change == "P", "present", 
                                                           if_else(gene_symbol == "PDGFRA" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                   NA),
         CDK46 =                           if_else(timing.primary == "TP" & surgery_number == 1,
                                                   if_else(gene_symbol == "CDK4/CDK6" & driver_status == "HLAMP" & driver_change == "P", "present", 
                                                           if_else(gene_symbol == "CDK4/CDK6" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                   NA),
         METPTPRZ1 =                       if_else(timing.primary == "TP" & surgery_number == 1,
                                                   if_else(gene_symbol == "MET/PTPRZ1" & driver_status == "HLAMP" & driver_change == "P", "present", 
                                                           if_else(gene_symbol == "MET/PTPRZ1" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                   NA),
         MYCMYCN =                         if_else(timing.primary == "TP" & surgery_number == 1,
                                                   if_else(gene_symbol == "MYC/MYCN" & driver_status == "HLAMP" & driver_change == "P", "present", 
                                                           if_else(gene_symbol == "MYC/MYCN" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                   NA),
         MGMT =                            if_else(surgery_number == 1,
                                                   if_else(mgmt_methylation == "U", "absent", "present"),
                                                   NA),
         glioma.type =                     if_else(idh_codel_subtype != "IDHmut-codel", "IDH-mutant Astrocytomas", "IDH-mutant Oligodendrogliomas"),
         ttr_status =                     1,
         ttr_time =                       if_else(surgery_number == 2 & surgical_interval_mo !=0,
                                                  surgical_interval_mo/12,
                                                  NA)
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
            ttr_status =     ttr_status[1],
            ttr_time =       ttr_time[surgery_number == 2]
  )%>%
  ungroup()%>%
  filter(!duplicated(case_barcode))

#################################################################################################################

# Univariable time-to-recurrence analysis of patients with IDH-mutant astrocytomas
all.astros <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Astrocytomas")

# summary(coxph(Surv(ttr_time, ttr_status) ~ MSH6, data = all.astros ), conf.int=0.9) ## only wildtype

summary(coxph(Surv(ttr_time, ttr_status) ~ NOTCH1, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ NOTCH1, data = all.astros))
survfit(Surv(ttr_time, ttr_status) ~ NOTCH1, data = all.astros )

summary(coxph(Surv(ttr_time, ttr_status) ~ PIK3CA, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ PIK3CA, data = all.astros))
survfit(Surv(ttr_time, ttr_status) ~ PIK3CA, data = all.astros )

# summary(coxph(Surv(ttr_time, ttr_status) ~ PIK3R1, data = all.astros ), conf.int=0.9) ## only wildtype

summary(coxph(Surv(ttr_time, ttr_status) ~ CDKN2AB, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ CDKN2AB, data = all.astros)) 
survfit(Surv(ttr_time, ttr_status) ~ CDKN2AB, data = all.astros)

summary(coxph(Surv(ttr_time, ttr_status) ~ PTEN, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ PTEN, data = all.astros)) 
survfit(Surv(ttr_time, ttr_status) ~ PTEN, data = all.astros)

summary(coxph(Surv(ttr_time, ttr_status) ~ ATRX, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ ATRX, data = all.astros))
survfit(Surv(ttr_time, ttr_status) ~ ATRX, data = all.astros)

summary(coxph(Surv(ttr_time, ttr_status) ~ CCND2, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ CCND2, data = all.astros))
survfit(Surv(ttr_time, ttr_status) ~ CCND2, data = all.astros)

summary(coxph(Surv(ttr_time, ttr_status) ~ EGFR, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ EGFR, data = all.astros)) 
survfit(Surv(ttr_time, ttr_status) ~ EGFR, data = all.astros)

summary(coxph(Surv(ttr_time, ttr_status) ~ PDGFRA, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ PDGFRA, data = all.astros))
survfit(Surv(ttr_time, ttr_status) ~ PDGFRA, data = all.astros)

summary(coxph(Surv(ttr_time, ttr_status) ~ CDK46, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ CDK46, data = all.astros))
survfit(Surv(ttr_time, ttr_status) ~ CDK46, data = all.astros)

summary(coxph(Surv(ttr_time, ttr_status) ~ METPTPRZ1, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ METPTPRZ1, data = all.astros)) ## PH assumption violated
survfit(Surv(ttr_time, ttr_status) ~ METPTPRZ1, data = all.astros)

summary(coxph(Surv(ttr_time, ttr_status) ~ MYCMYCN, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ MYCMYCN, data = all.astros)) 
survfit(Surv(ttr_time, ttr_status) ~ MYCMYCN, data = all.astros)

summary(coxph(Surv(ttr_time, ttr_status) ~ MGMT, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ MGMT, data = all.astros)) 
survfit(Surv(ttr_time, ttr_status) ~ MGMT, data = all.astros)

#################################################################################################################

# Univariable time-to-recurrence analysis of patients with IDH-mutant oligodendrogliomas
all.oligos <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Oligodendrogliomas")

# summary(coxph(Surv(ttr_time, ttr_status) ~ MSH6, data = all.oligos ), conf.int=0.9) ## only wildtype

summary(coxph(Surv(ttr_time, ttr_status) ~ NOTCH1, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ NOTCH1, data = all.oligos)) 
survfit(Surv(ttr_time, ttr_status) ~ NOTCH1, data = all.oligos )

summary(coxph(Surv(ttr_time, ttr_status) ~ PIK3CA, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ PIK3CA, data = all.oligos)) 
survfit(Surv(ttr_time, ttr_status) ~ PIK3CA, data = all.oligos )

summary(coxph(Surv(ttr_time, ttr_status) ~ PIK3R1, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ PIK3R1, data = all.oligos)) 
survfit(Surv(ttr_time, ttr_status) ~ PIK3R1, data = all.oligos )

# summary(coxph(Surv(ttr_time, ttr_status) ~ CDKN2AB, data = all.oligos ), conf.int=0.9) ## all non-HD

# summary(coxph(Surv(ttr_time, ttr_status) ~ PTEN, data = all.oligos ), conf.int=0.9) ## all non-HD

# summary(coxph(Surv(ttr_time, ttr_status) ~ ATRX, data = all.oligos ), conf.int=0.9) ## all non-HD

# summary(coxph(Surv(ttr_time, ttr_status) ~ CCND2, data = all.oligos ), conf.int=0.9) ## all non-amp

summary(coxph(Surv(ttr_time, ttr_status) ~ EGFR, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ EGFR, data = all.oligos))
survfit(Surv(ttr_time, ttr_status) ~ EGFR, data = all.oligos)

# summary(coxph(Surv(ttr_time, ttr_status) ~ PDGFRA, data = all.oligos ), conf.int=0.9)  ## all non-amp

# summary(coxph(Surv(ttr_time, ttr_status) ~ CDK46, data = all.oligos ), conf.int=0.9) ## all non-amp

summary(coxph(Surv(ttr_time, ttr_status) ~ METPTPRZ1, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(ttr_time, ttr_status) ~ METPTPRZ1, data = all.oligos))  
survfit(Surv(ttr_time, ttr_status) ~ METPTPRZ1, data = all.oligos)

# summary(coxph(Surv(ttr_time, ttr_status) ~ MYCMYCN, data = all.oligos ), conf.int=0.9)  ## all non-amp

# summary(coxph(Surv(ttr_time, ttr_status) ~ MGMT, data = all.oligos ), conf.int=0.9)  ## methylated

#################################################################################################################