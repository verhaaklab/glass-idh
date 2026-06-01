# GLASS-I: Post-Recurrence Survival Hypermutation Status
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

# load data
molecular.data <- read.delim("~/driver_changes.txt")%>%                       ## created with GLASS-I_Drivers_Table.R
  mutate(timing.recurrent = substr(tumor_pair_barcode, 20,21),
         timing.primary = substr(tumor_pair_barcode, 14,15))%>%
  select(case_barcode, HM_group,timing.recurrent, timing.primary)%>%
  distinct(.)

prior.treatment <- read_csv("~/prior_treatment_gold_set.csv")%>%              ## created with GLASS-I_Prior_Treatment_Table.R
  select(-case_barcode, -surgery_number)

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")             ## data available on Synapse (Synapse ID: syn17038081)

surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
  mutate(timing.surgery = substr(sample_barcode, 14,15),
         tumor_type = if_else(idh_codel_subtype == "IDHmut-noncodel", "IDH-mutant Astrocytoma", "IDH-mutant Oligodendroglioma")  ) %>%
  filter(idh_codel_subtype != "IDHwt")%>%
  left_join(prior.treatment, by ="sample_barcode")

gold.set <- dbGetQuery(con, "SELECT * FROM analysis.gold_set")%>%
  mutate(timing.recurrent = substr(tumor_pair_barcode, 20,21),
         timing.primary = substr(tumor_pair_barcode, 14,15) )%>%
  semi_join(surgeries, by = "case_barcode")

survival.data <- dbGetQuery(con, "SELECT * FROM clinical.cases")[,c(2,5:8)]%>%
  semi_join(gold.set, by="case_barcode")%>%
  mutate(os_status =         if_else(case_vital_status == "dead", 1, 0),
         os_time =           (case_overall_survival_mo)/12
  )

dbDisconnect(con)

hm.data <- survival.data%>%
  left_join(surgeries, by ="case_barcode")%>%
  left_join(molecular.data, by ="case_barcode")%>%
  filter(timing.surgery == timing.recurrent)%>%
  filter(ALK_CHEMO == 1)%>%
  mutate(prs_time = os_time - (surgical_interval_mo/12))

hm.data.oligo <- hm.data %>%
  filter(tumor_type == "IDH-mutant Oligodendroglioma",)

ggsurvplot(survfit(Surv(prs_time, os_status) ~ HM_group, data = hm.data.oligo),
           hm.data.oligo, 
           palette = c("#377EB8", "gray60"), 
           short.panel.labs = T,
           ggtheme = theme_linedraw(),
           legend.position= "top",
           legend.title = "",
           legend.labs = c("Hypermutant", "Non-Hypermutant"),
           pval.coord = c(7.5, 0.85),
           surv.median.line = "v",
           break.x.by = 2.5,
           font.legend = c(13, "plain", "black"),
           censor.shape = "|",
           censor.size = 4,
           xlab ="Post-recurrence survival (years)",
           xlim = c(0,13),
           pval = T,
           risk.table.col = "strata",
           risk.table.y.text =F,
           risk.table = "nrisk_cumcensor")

hm.data.astro <- hm.data %>%
  filter(tumor_type == "IDH-mutant Astrocytoma",)

ggsurvplot(survfit(Surv(prs_time, os_status) ~ HM_group, data = hm.data.astro),
           hm.data.astro, 
           palette = c("#377EB8", "gray60"), 
           short.panel.labs = T,
           ggtheme = theme_linedraw(),
           legend.position= "top",
           legend.title = "",
           legend.labs = c("Hypermutant", "Non-Hypermutant"),
           pval.coord = c(7.5, 0.85),
           surv.median.line = "v",
           break.x.by = 2.5,
           font.legend = c(13, "plain", "black"),
           censor.shape = "|",
           censor.size = 4,
           xlab ="Post-recurrence survival (years)",
           xlim = c(0,12.5),
           pval = T,
           risk.table.col = "strata",
           risk.table.y.text =F,
           risk.table = "nrisk_cumcensor")