# GLASS-I: Post-Recurrence Survival Driver CNAs
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

cnv.data.recurrent.nhm <-    survival.data%>%
  left_join(surgeries, by ="case_barcode")%>%
  left_join(molecular.data, by ="case_barcode")%>%
  filter(surgery_number != 1) %>% 
  filter(timing.recurrent == timing.surgery) %>% 
  filter((gene_symbol == "CDKN2A/CDKN2B" |gene_symbol == "Any Amplification") & HM_group == "NHM" & idh_codel_subtype == "IDHmut-noncodel")%>%
  mutate(os_time =         os_time - (surgical_interval_mo/12),
         variable_status = if_else(driver_status == "No_AMP/HD", "No Amplification or Homozygous Deletion",
                                   if_else(driver_change == "P", "No Amplification or Homozygous Deletion","Amplification or Homozygous Deletion")),
         gene_symbol =     if_else(gene_symbol== "CDKN2A/CDKN2B" , "CDKN2A/B", "Any Focal Amplifications"),
         tumor_type =      "IDH-mutant \nAstrocytoma")

p1 <- ggsurvplot(survfit( Surv(os_time, os_status) ~ variable_status, data = cnv.data.recurrent.nhm ),
                 cnv.data.recurrent.nhm, 
                 facet.by = c("gene_symbol") ,
                 palette = c("#377EB8", "gray60"), 
                 ggtheme = theme_linedraw(),
                 legend.title = "CNAs at recurrence",
                 legend.labs = c("Amp. or HD", "No Amp. or HD"),
                 short.panel.labs = T,
                 pval.coord = c(10, 0.95),
                 surv.median.line = "v",
                 break.x.by = 5,
                 font.legend = c(13, "plain", "black"),
                 censor.shape = "|",
                 censor.size = 4,
                 panel.labs.font = list(face = "bold", size =13),
                 xlab ="Post-recurrence survival (years)",
                 xlim = c(0,15),
                 pval = T)

table(cnv.data.recurrent.nhm$gene_symbol, cnv.data.recurrent.nhm$variable_status)
coxph(Surv(os_time, os_status) ~ variable_status + gene_symbol, data = cnv.data.recurrent.nhm )
survfit(Surv(os_time, os_status) ~ variable_status + gene_symbol, data = cnv.data.recurrent.nhm )

cnv.data.recurrent.nhm.grade <- cnv.data.recurrent.nhm%>%
  group_by(case_barcode)%>%
  summarise(
    Grade =           if_else(any(gene_symbol == "CDKN2A/B" & variable_status  == "Amplification or Homozygous Deletion"), "Grade 4",
                              if_else(any(grade == "II"),"Grade 2/3",
                                      if_else(any(grade == "III"),"Grade 2/3",
                                              if_else(any(grade == "IV"),"Grade 4", 
                                                      NA))))
  )%>%
  ungroup()%>%
  select(case_barcode, Grade)

cnv.data.recurrent.nhm.full <- cnv.data.recurrent.nhm%>% 
  left_join(cnv.data.recurrent.nhm.grade, by = "case_barcode")%>%
  filter(!is.na(Grade))

p2 <- ggsurvplot(survfit( Surv(os_time, os_status) ~ variable_status, data = cnv.data.recurrent.nhm.full ),
           cnv.data.recurrent.nhm.full, 
           facet.by = c("Grade", "gene_symbol") ,
           palette = c("#377EB8", "gray60"), 
           ggtheme = theme_linedraw(),
           legend.title = "CNAs at recurrence",
           legend.labs = c("Amp. or HD", "No Amp. or HD"),
           short.panel.labs = T,
           pval.coord = c(10, 0.95),
           surv.median.line = "v",
           break.x.by = 5,
           font.legend = c(13, "plain", "black"),
           censor.shape = "|",
           censor.size = 4,
           panel.labs.font = list(face = "bold", size =13),
           xlab ="Post-recurrence survival (years)",
           xlim = c(0,15),
           pval = T)

# Focal Oncogene Amplification or HD-CDKN2A/B in NHM IDH-mutant gliomas
amp.cdkn2ab.data.astro <- survival.data%>%
  left_join(surgeries, by ="case_barcode")%>%
  left_join(molecular.data, by ="case_barcode")%>%
  filter(timing.recurrent == timing.surgery) %>% 
  filter(surgery_number != 1) %>% 
  filter((gene_symbol == "CDKN2A/CDKN2B" |gene_symbol == "Any Amplification") & HM_group == "NHM")%>%
  mutate(variable_status = if_else(driver_status == "No_AMP/HD", "No Amplification or Homozygous Deletion",
                                   if_else(driver_change == "P", "No Amplification or Homozygous Deletion","Amplification or Homozygous Deletion")),
         gene_symbol =     if_else(gene_symbol== "CDKN2A/CDKN2B" , "CDKN2A/B", "Any Focal Amplifications")
  )%>%
  group_by(case_barcode)%>%
  summarise(prs_time =      os_time[1] - (surgical_interval_mo[1]/12),
            os_status =     os_status[1],
            tumor_type =    if_else(idh_codel_subtype[1] == "IDHmut-noncodel", "IDH-mutant Astrocytoma", "IDH-mutant Oligodendroglioma"),
            either =        if_else(any(variable_status == "Amplification or Homozygous Deletion"), "Amplification or Homozygous Deletion",
                                    if_else(any(variable_status == "No Amplification or Homozygous Deletion"),"No Amplification or Homozygous Deletion", NA)),
            CDKN2AB =       if_else(any(driver_status == "HLDEL" & (driver_change == "R"| driver_change == "S")), "present", "absent"),
            Focal_amp =     if_else(any(driver_status == "HLAMP" & (driver_change == "R"| driver_change == "S")), "present", "absent"),
  )%>%
  ungroup() %>%
  filter(tumor_type == "IDH-mutant Astrocytoma")

p3 <- ggsurvplot(survfit(Surv(prs_time, os_status) ~ either, data = amp.cdkn2ab.data.astro),
                 amp.cdkn2ab.data.astro, 
                 palette = c("#377EB8", "gray60"), 
                 short.panel.labs = T,
                 ggtheme = theme_linedraw(),
                 legend.position= "top",
                 legend.title = "",
                 legend.labs = c("CNA+", "CNA-"),
                 pval.coord = c(7.5, 0.85),
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

survfit(Surv(prs_time, os_status) ~ either, data = amp.cdkn2ab.data.astro)
coxph(Surv(prs_time, os_status) ~ either, data = amp.cdkn2ab.data.astro)

amp.cdkn2ab.data.astro.no.mixed <- amp.cdkn2ab.data.astro%>%
  filter(!(CDKN2AB == "present" & Focal_amp== "present"))

p4 <- ggsurvplot(survfit(Surv(prs_time, os_status) ~ CDKN2AB + Focal_amp, data = amp.cdkn2ab.data.astro.no.mixed),
                 amp.cdkn2ab.data.astro.no.mixed, 
                 palette = c("gray60", "#377EB8", "#c46666"),
                 short.panel.labs = T,
                 ggtheme = theme_linedraw(),
                 legend.position= "top",
                 legend.title = "",
                 # legend.labs = c("CNA+", "CNA-"),
                 pval.coord = c(7.5, 0.85),
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

survfit(Surv(prs_time, os_status) ~ CDKN2AB + Focal_amp, data = amp.cdkn2ab.data.astro.no.mixed)
coxph(Surv(prs_time, os_status) ~ CDKN2AB + Focal_amp, data = amp.cdkn2ab.data.astro.no.mixed)

p1
p2
p3
p4