#################################################################################################################
# Content: Kaplan-Meier estimates of IDH-mutant gliomas divided by hypermutation status
# Author: Mircea Tesileanu
# Date: 2025.07.03
# R-version: 4.3.2
# Figures: 3A and 3B
#################################################################################################################

# Clean start and load libraries
rm(list = ls())
gc()

library(tidyverse)    # v2.0.0
library(RPostgres)    # v1.4.7
library(survival)     # v3.7-0
library(survminer)    # v0.5.0

#################################################################################################################

# Function to change Alk_Chemo column to any previous treatment with alkylating chemotherapy
modify_alk_chemo <- function(df) {
  df <- df %>%
    arrange(case_barcode, surgery_number) %>%
    group_by(case_barcode) %>%
    mutate(
      Alk_Chemo = if_else(cummax(replace_na(alk_chemo, 0)) == 1, 1, alk_chemo)
    )
  return(df)
}

#################################################################################################################

# Load relevant data 
molecular.data <-  read.delim("/path/to/driver_changes_14042025.txt")%>%
  mutate(timing.recurrent = substr(tumor_pair_barcode, 20,21),
         timing.primary = substr(tumor_pair_barcode, 14,15))%>%
  select(case_barcode, HM_group,timing.recurrent, timing.primary)%>%
  distinct(.)

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
         timing.surgery =   substr(sample_barcode, 14,15),
         tumor_type =       if_else(idh_codel_subtype == "IDHmut-noncodel", "IDH-mutant Astrocytoma", "IDH-mutant Oligodendroglioma"),
         alk_chemo=         case_when(treatment_alkylating_agent == T  &  (treatment_concurrent_tmz == T|treatment_concurrent_tmz == F|is.na(treatment_concurrent_tmz))~ 1,
                                      treatment_concurrent_tmz == T  &  (treatment_alkylating_agent == T |treatment_alkylating_agent == F|is.na(treatment_alkylating_agent))~ 1,
                                      treatment_alkylating_agent == F  &  (treatment_concurrent_tmz == F |  is.na(treatment_concurrent_tmz))~ 0,
                                      is.na(treatment_alkylating_agent)  &  treatment_concurrent_tmz == F~ 0)
         )%>%
  filter(idh_codel_subtype != "IDHwt")

silver.set <- dbGetQuery(con, "SELECT * FROM analysis.silver_set")%>%
  semi_join(surgeries, by = "case_barcode")

survival.data <- dbGetQuery(con, "SELECT * FROM clinical.cases")[,c(2,5:8)]%>%
  semi_join(silver.set, by="case_barcode")%>%
  mutate(os_status =         if_else(case_vital_status == "dead", 1, 0),
         os_time =           (case_overall_survival_mo)/12
  )

dbDisconnect(con)

surgeries <- modify_alk_chemo(surgeries)

hm.data <- survival.data%>%
  left_join(surgeries, by ="case_barcode")%>%
  left_join(molecular.data, by ="case_barcode")%>%
  filter(timing.surgery == timing.recurrent)%>%
  filter(timing.primary == "TP")%>%
  filter(Alk_Chemo == 1)%>%
  mutate(prs_time = os_time - (surgical_interval_mo/12))

#################################################################################################################

# Plot Figure 3A: Kaplan-Meier estimate of IDH-mutant astrocytomas divided by hypermutation status
hm.data.astro <- hm.data %>%
  filter(tumor_type == "IDH-mutant Astrocytoma")

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
           xlim = c(0,10),
           pval = T,
           risk.table.col = "strata",
           risk.table.y.text =F,
           risk.table = "nrisk_cumcensor")

#################################################################################################################

# Plot Figure 3B: Kaplan-Meier estimate of IDH-mutant oligodendrogliomas divided by hypermutation status
hm.data.oligo <- hm.data %>%
  filter(tumor_type == "IDH-mutant Oligodendroglioma")

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

#################################################################################################################