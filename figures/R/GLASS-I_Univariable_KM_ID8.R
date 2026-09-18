# GLASS-I: Post-Recurrence Survival ID8
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
id8.data <-  read_csv("~/sbs11_sbs119_id8_Rprop.csv")%>%                      ## created with mutational_signature_pipeline.sh
  mutate(sig_prop_id8 = ID8_R)%>%
  select(case_barcode, sig_prop_id8)

prior.treatment <- read_csv("~/prior_treatment_gold_set.csv")%>%              ## created with GLASS-I_Prior_Treatment_Table.R
  select(-case_barcode, -surgery_number)

hm.data <- read.delim("~/driver_changes.txt")%>%                              ## created with GLASS-I_Drivers_Table.R
  select(case_barcode, HM_group)%>%
  distinct(.)

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")             ## data available on Synapse (Synapse ID: syn17038081)

surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
  mutate(timing.surgery = substr(sample_barcode, 14,15)) %>%
  filter(idh_codel_subtype != "IDHwt")%>%
  left_join(prior.treatment, by ="sample_barcode")

gold.set <- dbGetQuery(con, "SELECT * FROM analysis.gold_set")%>%
  mutate(timing.recurrent = substr(tumor_pair_barcode, 20,21),
         timing.primary = substr(tumor_pair_barcode, 14,15) )%>%
  filter(str_detect(tumor_pair_barcode, "-WGS"))%>%
  semi_join(surgeries, by = "case_barcode")%>%
  select(case_barcode, timing.recurrent, timing.primary, tumor_pair_barcode)

survival.data <- dbGetQuery(con, "SELECT * FROM clinical.cases")[,c(2,5:8)]%>%
  right_join(gold.set, by="case_barcode")%>%
  mutate(os_status =         if_else(case_vital_status == "dead", 1, 0),
         os_time =           (case_overall_survival_mo)/12,
  )

initial.tmb.data <-  dbGetQuery(con, "SELECT * FROM analysis.tumor_clinical_comparison_2024") %>%
  filter(str_detect(tumor_pair_barcode, "-WGS"))%>%
  filter((sbs11_r_hypermutation == F & mutation_burden_a <10) |sbs11_r_hypermutation == T)%>%
  semi_join(gold.set, by = "case_barcode")

dbDisconnect(con)

id8 <- survival.data%>%
  semi_join(initial.tmb.data, by = "case_barcode")%>%
  left_join(surgeries, by ="case_barcode")%>%
  left_join(id8.data, by ="case_barcode")%>%
  left_join(hm.data, by ="case_barcode")%>%
  filter(timing.surgery == timing.recurrent)%>%
  filter(HM_group == "NHM")%>%
  mutate(id8 = if_else(sig_prop_id8 > mean(sig_prop_id8), T, F),
         prs_time =          os_time - (surgical_interval_mo/12))

id8.astro <- id8%>%
  filter(idh_codel_subtype == "IDHmut-noncodel")

ggsurvplot(survfit(Surv(prs_time, os_status) ~ id8, data = id8.astro ),
           id8.astro, 
           ggtheme = theme_linedraw(),
           palette = c("gray60", "#377EB8"),
           legend.labs = c("ID8 low", "ID8 high"),
           legend.title = "NHM Astro.",
           short.panel.labs = T,
           pval.coord = c(10, 0.85),
           surv.median.line = "v",
           break.x.by = 2.5,
           font.legend = c(13, "plain", "black"),
           censor.shape = "|",
           censor.size = 4,
           xlab ="Post-recurrence survival (years)",
           xlim = c(0,15),
           pval = T,
           risk.table.col = "strata",
           risk.table.y.text =F,
           risk.table = "nrisk_cumcensor")

id8.oligo <- id8%>%
  filter(idh_codel_subtype == "IDHmut-codel")

ggsurvplot(survfit(Surv(prs_time, os_status) ~ id8, data = id8.oligo ),
           id8.oligo, 
           ggtheme = theme_linedraw(),
           palette = c("gray60", "#377EB8"),
           legend.labs = c("ID8 low", "ID8 high"),
           legend.title = "NHM Oligo.",
           short.panel.labs = T,
           pval.coord = c(10, 0.85),
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