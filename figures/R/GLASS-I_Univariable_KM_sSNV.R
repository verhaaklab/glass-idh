# GLASS-I: Post-Recurrence Survival Driver sSNVs
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(tidyverse)
library(RPostgres)
library(survival)
library(survminer)
library(patchwork)

# load data
molecular.data <- read.delim("~/driver_changes.txt")%>%                       ## created with GLASS-I_Drivers_Table.R
  select(-idh_codel_subtype) 

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")             ## data available on Synapse (Synapse ID: syn17038081)

surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
  mutate(timing.surgery = substr(sample_barcode, 14,15)) %>%
  filter(idh_codel_subtype != "IDHwt")

gold.set <- dbGetQuery(con, "SELECT * FROM analysis.gold_set")%>%
  mutate(timing.recurrent = substr(tumor_pair_barcode, 20,21),
         timing.primary = substr(tumor_pair_barcode, 14,15) )%>%
  semi_join(surgeries, by = "case_barcode")

survival.data <- dbGetQuery(con, "SELECT * FROM clinical.cases")[,c(2,5:8)]%>%
  right_join(gold.set, by="case_barcode")%>%
  mutate(os_status =         if_else(case_vital_status == "dead", 1, 0),
         os_time =           (case_overall_survival_mo)/12
  )

dbDisconnect(con)

# PI3K or NOTCH1 mutations in NHM IDH-mutant oligodendrogliomas
pi3k.notch1.data <- survival.data%>%
  left_join(surgeries, by ="case_barcode")%>%
  left_join(molecular.data, by ="case_barcode")%>%
  filter(timing.recurrent == timing.surgery) %>% 
  filter(surgery_number != 1) %>% 
  filter((gene_symbol == "PI3K" |gene_symbol == "NOTCH1") & HM_group == "NHM")%>%
  mutate(variable_status = if_else(driver_status == "wildtype", "wildtype",
                                   if_else(driver_change == "P", "wildtype","mutant"))
  )%>%
  group_by(case_barcode)%>%
  summarise(prs_time =      os_time[1] - (surgical_interval_mo[1]/12),
            os_status =     os_status[1],
            tumor_type =    if_else(idh_codel_subtype[1] == "IDHmut-noncodel", "IDH-mutant Astrocytoma", "IDH-mutant Oligodendroglioma"),
            either =        if_else(any(variable_status == "mutant"), "mutant",
                                    if_else(any(variable_status == "wildtype"),"wildtype", NA))
  )%>%
  ungroup() %>%
  filter(tumor_type == "IDH-mutant Oligodendroglioma",)

ggsurvplot(survfit(Surv(prs_time, os_status) ~ either, data = pi3k.notch1.data),
           pi3k.notch1.data, 
           palette = c("#377EB8", "gray60"), 
           short.panel.labs = T,
           ggtheme = theme_linedraw(),
           legend.position= "top",
           legend.title = "",
           legend.labs = c("SNV+", "SNV-"),
           pval.coord = c(7.5, 0.85),
           surv.median.line = "v",
           break.x.by = 2.5,
           font.legend = c(13, "plain", "black"),
           censor.shape = "|",
           censor.size = 4,
           xlab ="Post-recurrence survival (years)",
           xlim = c(0,17.5),
           pval = T,
           risk.table.col = "strata",
           risk.table.y.text =F,
           risk.table = "nrisk_cumcensor")