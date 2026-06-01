# GLASS-I: Post-Recurrence Survival Drivers GENIE
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# Load libraries
library(tidyverse)
library(survival)
library(survminer)

# Load data
data_clinical_patient <- read.delim("~/genie_v18.5/data_clinical_patient.txt", skip=4) ## data available on Synapse (Synapse ID: syn7222066)
data_clinical_sample <- read.delim("~/genie_v18.5/data_clinical_sample.txt", skip=4)   ## data available on Synapse (Synapse ID: syn7222066)
data_cna <- read.csv("~/genie_v18.5/filtered_CNA.csv")                                 ## data available on Synapse (Synapse ID: syn7222066)
data_gene_matrix <- read.delim("~/genie_v18.5/data_gene_matrix.txt")                   ## data available on Synapse (Synapse ID: syn7222066)
data_mutations_extended <- read.delim("~/genie_v18.5/data_mutations_extended.txt")     ## data available on Synapse (Synapse ID: syn7222066)
genomic_information <- read.delim("~/genie_v18.5/genomic_information.txt")             ## data available on Synapse (Synapse ID: syn7222066)
data_cna_full <- read.delim("~/genie_v18.5/data_cna_hg19.seg")                         ## data available on Synapse (Synapse ID: syn7222066)

# Define IDH-mutant glioma cases
canonical_idh_mutations <- data_mutations_extended%>%
  filter(Hugo_Symbol == "IDH1" & Protein_position == 132 | Hugo_Symbol == "IDH2" & Protein_position == 172)%>%
  filter(t_alt_count !=0 & Variant_Classification != "Silent")%>%
  add_count(HGVSp_Short, name = "mutation_count")%>%
  filter(mutation_count >5)%>%
  mutate(SAMPLE_ID = Tumor_Sample_Barcode)

idh_mutant_gliomas_samples <- data_clinical_sample  %>% 
  filter(CANCER_TYPE_DETAILED %in% c("Anaplastic Astrocytoma", "Anaplastic Oligodendroglioma",  "Anaplastic Oligoastrocytoma", "Astrocytoma",
                                     "Diffuse Astrocytoma", "Diffuse Glioma",  "Glioblastoma", "Glioblastoma Multiforme", "Glioma, NOS",
                                     "Gliosarcoma", "High-Grade Glioma, NOS", "Low-Grade Glioma, NOS","Oligoastrocytoma", 
                                     "Oligodendroglioma"
  )) %>%
  semi_join(canonical_idh_mutations, by = "SAMPLE_ID")

idh_mutant_gliomas_patients<- data_clinical_patient  %>% 
  semi_join(idh_mutant_gliomas_samples, by = "PATIENT_ID")

idh_mutant_gliomas_mutations <- data_mutations_extended %>%
  mutate(SAMPLE_ID = Tumor_Sample_Barcode,
         vaf = t_alt_count / (t_ref_count + t_alt_count))%>% 
  semi_join(idh_mutant_gliomas_samples, by = "SAMPLE_ID")%>%
  filter(vaf > 0.05 & t_alt_count >5)

pi3k_notch1_mutations <- idh_mutant_gliomas_mutations %>%
  filter(Hugo_Symbol == "PIK3CA"|Hugo_Symbol == "PIK3R1"|Hugo_Symbol == "NOTCH1")%>%
  filter(Variant_Classification != "Intron" & Variant_Classification != "Silent")%>%
  select(SAMPLE_ID)%>%
  distinct()%>%
  mutate(pi3k.notch1 = "mutant")

panel_size <- genomic_information %>%
  mutate(region_length = End_Position - Start_Position + 1) %>% 
  group_by(SEQ_ASSAY_ID) %>% 
  summarise(panel_size_bp = sum(region_length, na.rm = TRUE)) %>% 
  mutate(panel_mb = panel_size_bp / 1e6) %>% 
  filter(panel_size_bp > 499999)

# Calculate TMB for glioma cases
tmb_glioma <- idh_mutant_gliomas_mutations %>%
  left_join(select(idh_mutant_gliomas_samples, SAMPLE_ID, SEQ_ASSAY_ID), by = "SAMPLE_ID") %>%
  inner_join(panel_size, by = "SEQ_ASSAY_ID") %>%
  add_count(Tumor_Sample_Barcode, name = "mutation_count") %>%
  mutate(TMB = mutation_count / panel_mb) %>%
  select(SAMPLE_ID = Tumor_Sample_Barcode, TMB) %>% 
  distinct() %>% 
  filter(!is.na(TMB))

# Mutations
mut_list_glioma <- idh_mutant_gliomas_mutations %>%
  group_by(SAMPLE_ID) %>%
  summarise(mut_genes = paste(unique(Hugo_Symbol), collapse = ", ")) %>%
  ungroup() %>% 
  select(SAMPLE_ID, mut_genes)

# CNVs
cnv_glioma <- data_cna %>% 
  semi_join(idh_mutant_gliomas_samples, by = "SAMPLE_ID")

amp_or_hdcdkn2a_glioma <- cnv_glioma %>%
  mutate(amp_oncogene = ifelse(CCND2 == 2 |CDK4 ==2 |CDK6==2|EGFR==2|
                                 MET==2|MYC==2|MYCN==2|PDGFRA==2, "yes", "no"))%>%
  select(SAMPLE_ID, amp_oncogene, CDKN2A)

# Surv 
surv_glioma <- idh_mutant_gliomas_samples %>% 
  left_join(idh_mutant_gliomas_patients, by = "PATIENT_ID") %>% 
  filter(INT_CONTACT != "Unknown" & INT_CONTACT != "Not Collected" & INT_CONTACT != "" & INT_CONTACT != "<6570" & INT_CONTACT != ">32485" & INT_DOD != "Unknown") %>% 
  mutate(surv_time = if_else(INT_DOD == "Not Applicable"  ,  
                             (as.numeric(INT_CONTACT) - as.numeric(AGE_AT_SEQ_REPORT_DAYS))/365.25,
                             (as.numeric(INT_DOD) - as.numeric(AGE_AT_SEQ_REPORT_DAYS))/365.25),
         vital_status = ifelse(DEAD == "True", 1, 0)) %>% 
  filter(surv_time > 0) %>%
  select(PATIENT_ID, surv_time, vital_status)


#### Oligo's
oligo_data <- idh_mutant_gliomas_samples %>% 
  filter(CANCER_TYPE_DETAILED %in% c("Anaplastic Oligodendroglioma", "Oligodendroglioma", "Anaplastic Oligoastrocytoma", "Oligoastrocytoma")) %>%
  left_join(surv_glioma, by = "PATIENT_ID") %>% 
  left_join(tmb_glioma, by = "SAMPLE_ID") %>% 
  left_join(pi3k_notch1_mutations, by = "SAMPLE_ID") %>%
  left_join(mut_list_glioma, by = "SAMPLE_ID")%>%
  mutate(hypermutation = ifelse(TMB > 10, "hypermutated", "non-hypermutated"),
         group = ifelse(hypermutation == "hypermutated", "Hypermutated", 
                        ifelse(is.na(pi3k.notch1), "All WT", "PIK3CA/PIK3R1/NOTCH1")))%>%
  filter(!is.na(surv_time))

oligo_nontp53 <- oligo_data %>% 
  filter(!str_detect(mut_genes, "TP53")) %>% 
  filter(!str_detect(mut_genes, "ATRX"))%>%   
  filter(group == "PIK3CA/PIK3R1/NOTCH1"| group == "All WT")

oligo_nontp53_os <- oligo_nontp53%>%
  filter(SAMPLE_TYPE != "Metastasis")%>%
  group_by(PATIENT_ID) %>%
  slice_min(surv_time, with_ties = FALSE)

oligo_nontp53_prs <- oligo_nontp53%>%
  filter(SAMPLE_TYPE_DETAILED == "Local recurrence")%>%
  group_by(PATIENT_ID) %>%
  slice_min(surv_time, with_ties = FALSE)

ggsurvplot(survfit(Surv(surv_time, vital_status) ~ group, data = oligo_nontp53_os),
           oligo_nontp53_os, 
           palette = c( "gray60", "#377EB8"), 
           short.panel.labs = T,
           ggtheme = theme_linedraw(),
           legend.position= "top",
           legend.title = "",
           legend.labs = c("Neither sSNVs\n(n=306)", "PI3K and/or NOTCH1\nmutations (n=159)"),
           pval.coord = c(7.5, 0.85),
           surv.median.line = "v",
           break.x.by = 2.5,
           font.legend = c(13, "plain", "black"),
           censor.shape = "|",
           censor.size = 4,
           xlab ="Overall survival (years)",
           pval = T,
           risk.table.col = "strata",
           risk.table.y.text =F,
           risk.table = "nrisk_cumcensor")

ggsurvplot(survfit(Surv(surv_time, vital_status) ~ group, data = oligo_nontp53_prs),
           oligo_nontp53_prs, 
           palette = c( "gray60", "#377EB8"), 
           short.panel.labs = T,
           ggtheme = theme_linedraw(),
           legend.position= "top",
           legend.title = "",
           legend.labs = c("Neither sSNVs\n(n=77)", "PI3K and/or NOTCH1\nmutations (n=49)"),
           pval.coord = c(7.5, 0.85),
           surv.median.line = "v",
           break.x.by = 2.5,
           font.legend = c(13, "plain", "black"),
           censor.shape = "|",
           censor.size = 4,
           xlab ="Post-recurrence survival (years)",
           pval = T,
           risk.table.col = "strata",
           risk.table.y.text =F,
           risk.table = "nrisk_cumcensor")

#### Astro's
astro_data <- idh_mutant_gliomas_samples %>% 
  filter(CANCER_TYPE_DETAILED %in% c("Anaplastic Astrocytoma", "Astrocytoma", "Anaplastic Oligoastrocytoma", "Oligoastrocytoma", "Diffuse Astrocytoma")) %>%
  left_join(surv_glioma, by = "PATIENT_ID") %>% 
  left_join(tmb_glioma, by = "SAMPLE_ID") %>% 
  left_join(amp_or_hdcdkn2a_glioma, by = "SAMPLE_ID") %>% 
  left_join(mut_list_glioma, by = "SAMPLE_ID")%>%
  mutate(hypermutation = ifelse(TMB > 10, "hypermutated", "non-hypermutated"),
         group = ifelse(hypermutation == "hypermutated", "Hypermutated", 
                        ifelse(CDKN2A == -2 & amp_oncogene != "yes", "HD-CDKN2A/B", 
                               ifelse(CDKN2A == 0 & amp_oncogene == "yes", "Focal amp",  
                                      ifelse(CDKN2A != -2 & amp_oncogene != "yes", "no CNA", NA))))
  )%>%
  filter(!is.na(surv_time))

astro_tp53_or_atrx <- astro_data %>% 
  filter(str_detect(mut_genes, "TP53") | str_detect(mut_genes, "ATRX"))%>%   
  filter(group != "Hypermutated")%>%  
  filter(!(amp_oncogene == "yes" & CDKN2A == -2))

astro_tp53_or_atrx_os <- astro_tp53_or_atrx%>%
  filter(SAMPLE_TYPE != "Metastasis")%>%
  group_by(PATIENT_ID) %>%
  slice_min(surv_time, with_ties = FALSE)

astro_tp53_or_atrx_prs <- astro_tp53_or_atrx%>%
  filter(SAMPLE_TYPE_DETAILED == "Local recurrence")%>%
  group_by(PATIENT_ID) %>%
  slice_min(surv_time, with_ties = FALSE)

ggsurvplot(survfit(Surv(surv_time, vital_status) ~ group, data = astro_tp53_or_atrx_os),
           astro_tp53_or_atrx_os, 
           palette = c(  "#377EB8", "#c46666", "gray60"), 
           short.panel.labs = T,
           ggtheme = theme_linedraw(),
           legend.position= "top",
           legend.title = "",
           legend.labs = c("Focal Onc. Amp.only\n(n=57)", "HD-CDKN2A/B only\n(n=29)", "Neither CNAs\n(n=366)"),
           pval.coord = c(7.5, 0.85),
           surv.median.line = "v",
           break.x.by = 2.5,
           font.legend = c(13, "plain", "black"),
           censor.shape = "|",
           censor.size = 4,
           xlab ="Overall survival (years)",
           pval = T,
           risk.table.col = "strata",
           risk.table.y.text =F,
           risk.table = "nrisk_cumcensor")

ggsurvplot(survfit(Surv(surv_time, vital_status) ~ group, data = astro_tp53_or_atrx_prs),
           astro_tp53_or_atrx_prs, 
           palette = c(  "#377EB8", "#c46666", "gray60"), 
           short.panel.labs = T,
           ggtheme = theme_linedraw(),
           legend.position= "top",
           legend.title = "",
           legend.labs = c("Focal Onc. Amp.only\n(n=8)", "HD-CDKN2A/B only\n(n=10)", "Neither CNAs\n(n=49)"),
           pval.coord = c(7.5, 0.85),
           surv.median.line = "v",
           break.x.by = 2.5,
           font.legend = c(13, "plain", "black"),
           censor.shape = "|",
           censor.size = 4,
           xlab ="Post-recurrence survival (years)",
           pval = T,
           risk.table.col = "strata",
           risk.table.y.text =F,
           risk.table = "nrisk_cumcensor")