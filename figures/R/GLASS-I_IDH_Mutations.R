# GLASS-I: Canonical IDH mutations
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
library(patchwork)
library(ggbeeswarm)
library(rstatix)

# load data
molecular.data <- read.delim("~/idh_mutations.txt")%>%                       ## created with GLASS-I_IDH_Mutations_Table.R
  filter(gene_symbol!="IDH1/2")%>%
  select(-idh_codel_subtype, -tumor_pair_barcode)

hypermut.data <- molecular.data%>%
  dplyr::select(case_barcode, HM_group)%>%
  filter(!duplicated(case_barcode))

prior.treatment <- read_csv("~/prior_treatment_gold_set.csv")%>%              ## created with GLASS-I_Prior_Treatment_Table.R
  select(-case_barcode, -surgery_number)

driver.data <- read.delim("~/driver_changes.txt")%>%                          ## created with GLASS-I_Drivers_Table.R
  filter(gene_symbol == "Any Amplification" | gene_symbol == "CDKN2A/CDKN2B" | gene_symbol == "PI3K" |gene_symbol == "NOTCH1"|gene_symbol == "MSH6")%>%
  select(case_barcode, gene_symbol, driver_status, driver_change) 

fga.data <-  read.delim("~/fga_calc_annotated.txt")                           ## created with script_fga_calculation.R

values.tumor.type = c("Oligodendroglioma"="#29b4b4",
                      "Astrocytoma"="#800074")

values.mutations = c("IDH1 p.R132H"="#999999",
                     "Other IDH mutation"="#377EB8")

values.timing= c("Primary"="#cb2f66",
                 "Recurrent"="#30B3CA")

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")             ## data available on Synapse (Synapse ID: syn17038081)

surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
  mutate(timing.surgery = substr(sample_barcode, 14,15)) %>%
  filter(idh_codel_subtype != "IDHwt")%>%
  left_join(prior.treatment, by ="sample_barcode")

gold.set <- dbGetQuery(con, "SELECT * FROM analysis.gold_set")%>%
  mutate(timing.recurrent = substr(tumor_pair_barcode, 20,21),
         timing.primary = substr(tumor_pair_barcode, 14,15) )%>%
  semi_join(surgeries, by = "case_barcode")

initial.tmb.data <-  dbGetQuery(con, "SELECT * FROM analysis.tumor_clinical_comparison_2024") %>%
  filter(str_detect(tumor_pair_barcode, "-WGS"))%>%
  filter((sbs11_r_hypermutation == F & mutation_burden_a <10) |sbs11_r_hypermutation == T)%>%
  semi_join(gold.set, by = "case_barcode")

gold.set.samples<- dbGetQuery(con, "SELECT * FROM analysis.gold_set")%>%
  semi_join(surgeries, by = "case_barcode")%>%
  pivot_longer(
    cols = c(tumor_barcode_a, tumor_barcode_b),
    names_to = NULL,
    values_to = "aliquot_barcode"
  )%>%
  mutate(sample_barcode = substr(aliquot_barcode, 1, 15))

survival.data <- dbGetQuery(con, "SELECT * FROM clinical.cases")[,c(2,5:8)]%>%
  right_join(gold.set, by="case_barcode")%>%
  mutate(os_status =         if_else(case_vital_status == "dead", 1, 0),
         os_time =           (case_overall_survival_mo)/12
  )

dbDisconnect(con)

signature.data <-  read_csv("~/sbs11_sbs119_id8_Rprop.csv")%>%                ## created with mutational_signature_pipeline.sh
  semi_join(initial.tmb.data, by = "case_barcode")

id8.data <- signature.data %>%
  mutate(sig_prop_id8 = ID8_R)%>%
  select(case_barcode, sig_prop_id8)

sbs11.data <- signature.data %>%
  mutate(sig_prop_sbs11 = SBS11_R)%>%
  select(case_barcode, sig_prop_sbs11)

sbs119.data <- signature.data %>%
  mutate(sig_prop_sbs119 = SBS119_R)%>%
  select(case_barcode, sig_prop_sbs119)

surgeries <- surgeries%>%
  semi_join(gold.set.samples, by ="sample_barcode")

idh.data <-   survival.data%>%
  left_join(surgeries, by ="case_barcode")%>%
  left_join(molecular.data, by ="case_barcode")%>%
  group_by(case_barcode)%>%
  summarise(protein_change = protein_change[!is.na(protein_change)][1],
            idh1_mutation =     case_when(any(str_detect(protein_change, "p.R132")) ~ 1,
                                          (!all(str_detect(protein_change, "p.R132"))) ~ 0),
            idh2_mutation =     case_when(any(str_detect(protein_change, "p.R172")) ~ 1,
                                          (!all(str_detect(protein_change, "p.R172"))) ~ 0),
            idh1_r132h = case_when(any(str_detect(protein_change, "p.R132H")) ~ 1,
                                   (!all(str_detect(protein_change, "p.R132H"))) ~ 0),
            os_status   = os_status[1],
            os_time   = os_time[1],
            tumor_type = case_when(any(idh_codel_subtype == "IDHmut-codel") ~ "Oligodendroglioma",
                                   any(idh_codel_subtype == "IDHmut-noncodel") ~ "Astrocytoma"),
            idh_codel_subtype = idh_codel_subtype[1],
            wgs = if_else(str_detect(tumor_pair_barcode[1], "-WGS"), 1,0)
  )%>%
  ungroup()

table(idh.data$protein_change, idh.data$tumor_type)
fisher.test(table(idh.data$protein_change, idh.data$tumor_type))

wgs.idh.data <- idh.data%>%
  left_join(surgeries, by ="case_barcode")%>%
  left_join(id8.data, by ="case_barcode")%>%
  left_join(sbs11.data, by ="case_barcode")%>%
  left_join(sbs119.data, by ="case_barcode")%>%
  left_join(hypermut.data, by ="case_barcode")%>%
  filter(wgs == 1)%>%
  group_by(case_barcode)%>%
  arrange(surgery_number)%>%
  summarise(protein_change = protein_change[1],
            idh1_mutation =  idh1_mutation[1],
            idh2_mutation =  idh2_mutation[1],
            idh1_r132h =    case_when(idh1_r132h[1] == 1 ~ "IDH1 p.R132H",
                                      idh1_r132h[1] == 0 ~ "Other IDH mutation"),
            tumor_type = tumor_type[1],
            sig_prop_id8 = sig_prop_id8[1]*100,
            sig_prop_sbs11 = sig_prop_sbs11[1]*100,
            sig_prop_sbs119 = sig_prop_sbs119[1]*100,
            TMZ = TMZ[2],
            RT = RT[2],
            HM_group = HM_group[1]
  )%>%
  ungroup()%>%
  mutate(sig_prop_id8 =    if_else(RT==1 & !is.na(sig_prop_id8)& HM_group == "NHM", sig_prop_id8, NA),
         sig_prop_sbs11 =  if_else(TMZ==1 & !is.na(sig_prop_sbs11) & HM_group == "HM", sig_prop_sbs11, NA),
         sig_prop_sbs119 = if_else(TMZ==1 & !is.na(sig_prop_sbs119) & HM_group == "NHM", sig_prop_sbs119, NA)
  )

astro <-  idh.data%>%
  filter(idh_codel_subtype == "IDHmut-noncodel")

p1 <- ggsurvplot(survfit( Surv(os_time, os_status) ~ idh1_r132h, data = astro ),
                 astro, 
                 palette = c("#377EB8", "gray60"),
                 ggtheme = theme_linedraw(),
                 legend.title = "IDH1/2 mutation",
                 legend.labs = c("Other", "IDH1 p.R132H"),
                 short.panel.labs = T,
                 pval.coord = c(0, 0.14),
                 surv.median.line = "v",
                 break.x.by = 5,
                 font.legend = c(13, "plain", "black"),
                 censor.shape = "|",
                 censor.size = 4,
                 panel.labs.font = list(face = "bold", size =13),
                 xlab ="Overall survival (years)",
                 xlim = c(0,30),
                 pval = T,
                 risk.table.col = "strata",
                 risk.table.y.text =F,
                 risk.table = "nrisk_cumcensor")

oligo <- idh.data%>%
  filter(idh_codel_subtype == "IDHmut-codel")

p2 <- ggsurvplot(survfit( Surv(os_time, os_status) ~ idh1_r132h, data = oligo ),
                 oligo, 
                 palette = c("#377EB8", "gray60"),
                 ggtheme = theme_linedraw(),
                 legend.title = "IDH1/2 mutation",
                 legend.labs = c("Other", "IDH1 p.R132H"),
                 short.panel.labs = T,
                 pval.coord = c(0, 0.14),
                 surv.median.line = "v",
                 break.x.by = 5,
                 font.legend = c(13, "plain", "black"),
                 censor.shape = "|",
                 censor.size = 4,
                 panel.labs.font = list(face = "bold", size =13),
                 xlab ="Overall survival (years)",
                 xlim = c(0,30),
                 pval = T,
                 risk.table.col = "strata",
                 risk.table.y.text =F,
                 risk.table = "nrisk_cumcensor")

plot.idh.vs.sbs11 <- ggplot(wgs.idh.data)+
  geom_violin(aes(x=idh1_r132h, y= sig_prop_sbs11,  fill = idh1_r132h), color = "gray10", alpha = 0.4)+
  geom_beeswarm(aes(x=idh1_r132h, y= sig_prop_sbs11)) +
  facet_grid(~tumor_type)+
  theme_linedraw()+
  scale_color_manual(values = values.mutations)+
  scale_fill_manual(values = values.mutations)+
  ylab("SBS11 recurrence proportion (%)")+
  ylim(0,105)+
  xlab("")+
  ggtitle("Hypermutant/TMZ-treated") + 
  theme(legend.position = "none") +
  stat_compare_means(label.x = 1.05,
                     label.y = 102,  
                     method = "wilcox.test",
                     paired = F,
                     aes(x=idh1_r132h, y= sig_prop_sbs11, label = sprintf("Wilcoxon: p = %5.2f", as.numeric(..p.format..))))

plot.idh.vs.sbs119 <- ggplot(wgs.idh.data)+
  geom_violin(aes(x=idh1_r132h, y= sig_prop_sbs119,  fill = idh1_r132h), color = "gray10", alpha = 0.4)+
  geom_beeswarm(aes(x=idh1_r132h, y= sig_prop_sbs119)) +
  facet_grid(~tumor_type)+
  theme_linedraw()+
  scale_color_manual(values = values.mutations)+
  scale_fill_manual(values = values.mutations)+
  ylab("SBS119 recurrence proportion (%)")+
  ylim(0,105)+
  xlab("")+
  ggtitle("Non-Hypermutant/TMZ-treated") + 
  theme(legend.position = "none") +
  stat_compare_means(label.x = 1.05,
                     label.y = 102,  
                     method = "wilcox.test",
                     paired = F,
                     aes(x=idh1_r132h, y= sig_prop_sbs119, label = sprintf("Wilcoxon: p = %5.2f", as.numeric(..p.format..))))

plot.idh.vs.id8 <- ggplot(wgs.idh.data)+
  geom_violin(aes(x=idh1_r132h, y= sig_prop_id8,  fill = idh1_r132h), color = "gray10", alpha = 0.4)+
  geom_beeswarm(aes(x=idh1_r132h, y= sig_prop_id8)) +
  facet_grid(~tumor_type)+
  theme_linedraw()+
  scale_color_manual(values = values.mutations)+
  scale_fill_manual(values = values.mutations)+
  ylab("ID8 recurrence proportion (%)")+
  ylim(0,105)+
  xlab("")+
  ggtitle("Non-Hypermutant/RT-treated") + 
  theme(legend.position = "none") +
  stat_compare_means(label.x = 1.05,
                     label.y = 102,  
                     method = "wilcox.test",
                     paired = F,
                     aes(x=idh1_r132h, y= sig_prop_id8, label = sprintf("Wilcoxon: p = %5.2f", as.numeric(..p.format..))))

plot.idh.vs.sbs11 + plot.idh.vs.sbs119 + plot.idh.vs.id8+ 
  plot_annotation() +
  plot_layout(widths = c(2, 2))&
  theme(plot.tag = element_text(face = 'bold'),
        plot.title = element_text(face = 'bold', size =12, hjust = 0.5))

p1
p2

gold.set.aliquots <- gold.set%>%
  select(tumor_barcode_a, tumor_barcode_b) %>%
  pivot_longer(cols = everything(), values_to = "aliquot_barcode") %>%
  select(aliquot_barcode) 

idh.fga.data <- gold.set.aliquots%>%
  left_join(fga.data, by = "aliquot_barcode")%>%
  left_join(idh.data, by = "case_barcode")%>%
  mutate(idh1_r132h = if_else(idh1_r132h == 1, "IDH1 p.R132H", "Other IDH mutation"))%>%
  filter(!is.na(sample_type))%>%
  filter(HM_group == "NHM")%>%
  group_by(case_barcode)%>%
  summarise(idh1_r132h = idh1_r132h[1],
            tumor_type = tumor_type[1],
            delta.fga_norm = fga_norm[sample_type == "Recurrent"]- fga_norm[sample_type == "Primary"])%>%
  ungroup()

stat.test <- idh.fga.data %>%
  group_by(tumor_type) %>%
  rstatix::wilcox_test(delta.fga_norm ~ idh1_r132h, paired = F) %>%
  mutate(
    label = ifelse(
      p < 0.001,
      "Wilcoxon: p < 0.001",
      sprintf("Wilcoxon: p = %.3f", p)
    ),
    group1 = 1, 
    group2 = 2
  )

ggplot(idh.fga.data) +
  geom_violin(aes(x=idh1_r132h, y=delta.fga_norm, fill=idh1_r132h),
              color="gray10", alpha=0.4) +
  geom_beeswarm(aes(x=idh1_r132h, y=delta.fga_norm)) +
  facet_grid(~tumor_type) +
  theme_linedraw() +
  scale_color_manual(values=values.mutations) +
  scale_fill_manual(values=values.mutations) +
  ylab("FGA (Recurrence - Primary)") +
  ylim(-0.8, 0.8)+
  xlab("") +
  theme(legend.position = "none") +
  stat_pvalue_manual(stat.test, label = "label", y.position = 0.75)

idh.driver.data<- driver.data%>%
  left_join(idh.data, by = "case_barcode")%>%
  left_join(hypermut.data, by = "case_barcode")%>%
  mutate(idh1_r132h = if_else(idh1_r132h == 1, "IDH1 p.R132H", "Other IDH mutation"))

oligo.driver.data <- idh.driver.data%>%
  filter(tumor_type == "Oligodendroglioma")

oligo.amplification.data <- oligo.driver.data %>%
  filter(gene_symbol == "Any Amplification" & HM_group =="NHM")

table(oligo.amplification.data$idh1_r132h, oligo.amplification.data$driver_status)
fisher.test(table(oligo.amplification.data$idh1_r132h, oligo.amplification.data$driver_status)) 

oligo.CDKN2AB.data <- oligo.driver.data %>%
  filter(gene_symbol == "CDKN2A/CDKN2B" & HM_group =="NHM")

table(oligo.CDKN2AB.data$idh1_r132h, oligo.CDKN2AB.data$driver_status)
fisher.test(table(oligo.CDKN2AB.data$idh1_r132h, oligo.CDKN2AB.data$driver_status)) 

oligo.NOTCH1.data <- oligo.driver.data %>%
  filter(gene_symbol == "NOTCH1" & HM_group =="NHM")

table(oligo.NOTCH1.data$idh1_r132h, oligo.NOTCH1.data$driver_status)
fisher.test(table(oligo.NOTCH1.data$idh1_r132h, oligo.NOTCH1.data$driver_status)) 

oligo.PI3K.data <- oligo.driver.data %>%
  filter(gene_symbol == "PI3K"& HM_group =="NHM")

table(oligo.PI3K.data$idh1_r132h, oligo.PI3K.data$driver_status)
fisher.test(table(oligo.PI3K.data$idh1_r132h, oligo.PI3K.data$driver_status)) 

oligo.MSH6.data <- oligo.driver.data %>%
  filter(gene_symbol == "MSH6" & HM_group =="HM")

table(oligo.MSH6.data$idh1_r132h, oligo.MSH6.data$driver_status)
fisher.test(table(oligo.MSH6.data$idh1_r132h, oligo.MSH6.data$driver_status)) 

astro.driver.data <- idh.driver.data%>%
  filter(tumor_type == "Astrocytoma")

astro.amplification.data <- astro.driver.data %>%
  filter(gene_symbol == "Any Amplification" & HM_group =="NHM")

table(astro.amplification.data$idh1_r132h, astro.amplification.data$driver_status)
fisher.test(table(astro.amplification.data$idh1_r132h, astro.amplification.data$driver_status)) 

astro.CDKN2AB.data <- astro.driver.data %>%
  filter(gene_symbol == "CDKN2A/CDKN2B"& HM_group =="NHM")

table(astro.CDKN2AB.data$idh1_r132h, astro.CDKN2AB.data$driver_status)
fisher.test(table(astro.CDKN2AB.data$idh1_r132h, astro.CDKN2AB.data$driver_status)) 

astro.NOTCH1.data <- astro.driver.data %>%
  filter(gene_symbol == "NOTCH1" & HM_group =="NHM")

table(astro.NOTCH1.data$idh1_r132h, astro.NOTCH1.data$driver_status)
fisher.test(table(astro.NOTCH1.data$idh1_r132h, astro.NOTCH1.data$driver_status)) 

astro.PI3K.data <- astro.driver.data %>%
  filter(gene_symbol == "PI3K" & HM_group =="NHM")

table(astro.PI3K.data$idh1_r132h, astro.PI3K.data$driver_status)
fisher.test(table(astro.PI3K.data$idh1_r132h, astro.PI3K.data$driver_status)) 

astro.MSH6.data <- astro.driver.data %>%
  filter(gene_symbol == "MSH6" & HM_group =="HM")

table(astro.MSH6.data$idh1_r132h, astro.MSH6.data$driver_status)
fisher.test(table(astro.MSH6.data$idh1_r132h, astro.MSH6.data$driver_status)) 