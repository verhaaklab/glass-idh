# GLASS-I: Univariable Overall Survival Molecular Parameters
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
  mutate(MSH6 =                            if_else(surgery_number == 1 & timing.primary == timing.surgery,
                                                   if_else(gene_symbol == "MSH6" & driver_status == "mutant" & driver_change == "P", "present", 
                                                   if_else(gene_symbol == "MSH6" & driver_status == "mutant" & driver_change == "S", "present", "absent")),
                                                   NA),
         NOTCH1 =                          if_else(surgery_number == 1 & timing.primary == timing.surgery,
                                                   if_else(gene_symbol == "NOTCH1" & driver_status == "mutant" & driver_change == "P", "present", 
                                                   if_else(gene_symbol == "NOTCH1" & driver_status == "mutant" & driver_change == "S", "present", "absent")),
                                                   NA),
         PI3K =                          if_else(surgery_number == 1 & timing.primary == timing.surgery,
                                                   if_else(gene_symbol == "PI3K" & driver_status == "mutant" & driver_change == "P", "present", 
                                                   if_else(gene_symbol == "PI3K" & driver_status == "mutant" & driver_change == "S", "present", "absent")),
                                                   NA),
         CDKN2AB =                         if_else(surgery_number == 1 & timing.primary == timing.surgery,
                                                   if_else(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "P", "present", 
                                                   if_else(gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change == "S", "present", "absent")),
                                                   NA),
         Any_Amplification =                            if_else(surgery_number == 1 & timing.primary == timing.surgery,
                                                  if_else(gene_symbol == "Any Amplification" & driver_status == "HLAMP" & driver_change == "P", "present", 
                                                  if_else(gene_symbol == "Any Amplification" & driver_status == "HLAMP" & driver_change == "S", "present", "absent")),
                                                  NA),
         MGMT =                            if_else(surgery_number == 1 & timing.primary == timing.surgery,
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
            PI3K =         case_when(any(PI3K=="present") ~ "present", 
                                       any(PI3K=="absent") ~"absent"),
            CDKN2AB =        case_when(any(CDKN2AB=="present") ~ "present", 
                                       any(CDKN2AB=="absent") ~"absent"),
            Any_Amplification =           case_when(any(Any_Amplification=="present") ~ "present", 
                                       any(Any_Amplification=="absent") ~"absent"),
            MGMT =           case_when(any(MGMT=="present") ~ "present", 
                                       any(MGMT=="absent") ~"absent"),
            glioma.type =    glioma.type[1],
            os_status =      os_status[1],
            os_time =        os_time[1]
            )%>%
  ungroup()%>%
  filter(!duplicated(case_barcode))

#### IDH-mutant oligodendrogliomas ####
# All tumors
all.oligos <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Oligodendrogliomas")

# summary(coxph(Surv(os_time, os_status) ~ MSH6, data = all.oligos ), conf.int=0.9) ## only wildtype

summary(coxph(Surv(os_time, os_status) ~ NOTCH1, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ NOTCH1, data = all.oligos))
survfit(Surv(os_time, os_status) ~ NOTCH1, data = all.oligos )

summary(coxph(Surv(os_time, os_status) ~ PI3K, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ PI3K, data = all.oligos))
survfit(Surv(os_time, os_status) ~ PI3K, data = all.oligos )

# summary(coxph(Surv(os_time, os_status) ~ CDKN2AB, data = all.oligos ), conf.int=0.9) ## only non-HD

summary(coxph(Surv(os_time, os_status) ~ Any_Amplification, data = all.oligos ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ Any_Amplification, data = all.oligos))
survfit(Surv(os_time, os_status) ~ Any_Amplification, data = all.oligos)

# summary(coxph(Surv(os_time, os_status) ~ MGMT, data = all.oligos ), conf.int=0.9)  ## only methylated

#### IDH-mutant astrocytomas ####
# All tumors
all.astros <- patient.lvl.data %>% 
  filter(glioma.type == "IDH-mutant Astrocytomas")

# summary(coxph(Surv(os_time, os_status) ~ MSH6, data = all.astros ), conf.int=0.9) ## only wildtype

summary(coxph(Surv(os_time, os_status) ~ NOTCH1, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ NOTCH1, data = all.astros))
survfit(Surv(os_time, os_status) ~ NOTCH1, data = all.astros )

summary(coxph(Surv(os_time, os_status) ~ PI3K, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ PI3K, data = all.astros))
survfit(Surv(os_time, os_status) ~ PI3K, data = all.astros )

summary(coxph(Surv(os_time, os_status) ~ CDKN2AB, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ CDKN2AB, data = all.astros))
survfit(Surv(os_time, os_status) ~ CDKN2AB, data = all.astros)

summary(coxph(Surv(os_time, os_status) ~ Any_Amplification, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ Any_Amplification, data = all.astros))
survfit(Surv(os_time, os_status) ~ Any_Amplification, data = all.astros)

summary(coxph(Surv(os_time, os_status) ~ MGMT, data = all.astros ), conf.int=0.9)
cox.zph(coxph(Surv(os_time, os_status) ~ MGMT, data = all.astros)) ## PH assumption violated
survfit(Surv(os_time, os_status) ~ MGMT, data = all.astros)