# GLASS-I: Sankey Diagram Quality Control Pipeline
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(tidyverse)
library(RPostgres)
library(ComplexHeatmap)
library(ggsankey)

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")             ## data available on Synapse (Synapse ID: syn17038081)

aliquots <- dbGetQuery(con, "SELECT * FROM biospecimen.aliquots") %>%
  select(aliquot_barcode, sample_barcode, aliquot_analysis_type) %>%
  filter(!str_detect(sample_barcode, "NB|NM|M1"))

blocklist <- dbGetQuery(con, "SELECT * FROM analysis.blocklist")

surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
  filter(idh_codel_subtype != "IDHwt") %>%
  mutate(tumor_type = case_when(
    idh_codel_subtype == "IDHmut-codel"   ~ "Oligo.",
    idh_codel_subtype == "IDHmut-noncodel" ~ "Astro.",
    TRUE ~ NA_character_
  )) %>%
  select(case_barcode, sample_barcode, surgery_number, surgical_interval_mo, tumor_type)

survival.data <- dbGetQuery(con, "SELECT * FROM clinical.cases") %>%
  select(case_barcode, case_overall_survival_mo) %>%
  semi_join(surgeries, by = "case_barcode")

silver.set <- dbGetQuery(con, "SELECT * FROM analysis.silver_set") %>%
  semi_join(surgeries, by = "case_barcode") %>%
  select(tumor_barcode_a, tumor_barcode_b) %>%
  pivot_longer(cols = everything(), values_to = "aliquot_barcode") %>%
  select(aliquot_barcode) %>%
  mutate(silver = "yes")

gold.set <- dbGetQuery(con, "SELECT * FROM analysis.gold_set") %>%
  semi_join(surgeries, by = "case_barcode") %>%
  select(tumor_barcode_a, tumor_barcode_b) %>%
  pivot_longer(cols = everything(), values_to = "aliquot_barcode") %>%
  select(aliquot_barcode) %>%
  mutate(gold = "yes")

dbDisconnect(con)

# merge full data
full.data <- surgeries %>%
  left_join(survival.data, by = "case_barcode") %>%
  left_join(aliquots, by = "sample_barcode") %>%
  left_join(blocklist, by = "aliquot_barcode") %>%
  left_join(silver.set, by = "aliquot_barcode") %>%
  left_join(gold.set, by = "aliquot_barcode") %>%
  mutate(unique_identifier = paste0(case_barcode, "-", surgery_number))

# generate sankey_per_sample: one label per unique_identifier
sankey_per_sample_step1 <- full.data %>%
  arrange(silver, gold) %>%
  arrange(tumor_type, unique_identifier) %>%
  group_by(tumor_type, unique_identifier) %>% 
  mutate(
    preQC = if_else(row_number() == 1, "preQC", NA_character_),
    QC1   = if_else(any(str_detect(aliquot_analysis_type,"W")), "QC1", NA_character_),
    QC2   = if_else(any(!str_detect(aliquot_analysis_type,"RNA") & str_detect(fingerprint_exclusion,"allow")), "QC2", NA_character_),
    QC3   = if_else(any(!str_detect(aliquot_analysis_type,"RNA") & str_detect(coverage_exclusion,"allow") & str_detect(fingerprint_exclusion,"allow")), "QC3", NA_character_),
    QC4   = if_else(any(!str_detect(aliquot_analysis_type,"RNA") & str_detect(clinical_exclusion,"allow") & str_detect(coverage_exclusion,"allow") & str_detect(fingerprint_exclusion,"allow")), "QC4", NA_character_),
    QC5   = if_else(any(!str_detect(aliquot_analysis_type,"RNA") & str_detect(silver,"yes") & str_detect(clinical_exclusion,"allow") & str_detect(coverage_exclusion,"allow") & str_detect(fingerprint_exclusion,"allow")), "QC5", NA_character_),
    QC6   = if_else(any(!str_detect(aliquot_analysis_type,"RNA") & str_detect(silver,"yes") & str_detect(gold,"yes") & str_detect(cnv_exclusion,"allow|review") & str_detect(clinical_exclusion,"allow") & str_detect(coverage_exclusion,"allow") & str_detect(fingerprint_exclusion,"allow")), "QC6", NA_character_)
  ) %>%
  mutate(across(starts_with("QC"), ~ if_else(row_number() == 1, .x, NA_character_))) %>%
  ungroup()

full_sankey_data <- sankey_per_sample_step1 %>%
  group_by(tumor_type) %>%
  summarize(
    preQC_samples = paste0(tumor_type[1], " \n", n_distinct(case_barcode[!is.na(preQC)]),  " patients\n", sum(!is.na(preQC)),  " surgeries"),
    QC1_samples   = paste0(tumor_type[1], " \n", n_distinct(case_barcode[!is.na(QC1)]),  " patients\n", sum(!is.na(QC1)),  " surgeries"),
    QC2_samples   = paste0(tumor_type[1], " \n", n_distinct(case_barcode[!is.na(QC2)]),  " patients\n", sum(!is.na(QC2)),  " surgeries"),
    QC3_samples   = paste0(tumor_type[1], " \n", n_distinct(case_barcode[!is.na(QC3)]),  " patients\n", sum(!is.na(QC3)),  " surgeries"),
    QC4_samples   = paste0(tumor_type[1], " \n", n_distinct(case_barcode[!is.na(QC4)]),  " patients\n", sum(!is.na(QC4)),  " surgeries"),
    QC5_samples   = paste0(tumor_type[1], " \n", n_distinct(case_barcode[!is.na(QC5)]),  " patients\n", sum(!is.na(QC5)),  " surgeries"),
    QC6_samples   = paste0(tumor_type[1], " \n", n_distinct(case_barcode[!is.na(QC6)]),  " patients\n", sum(!is.na(QC6)),  " surgeries"),
  ) %>%
  ungroup()

sankey_per_sample <- sankey_per_sample_step1 %>%
  left_join(full_sankey_data, by = "tumor_type") %>%
  mutate(preQC = if_else(!is.na(preQC), preQC_samples, NA_character_),
         QC1 = if_else(!is.na(QC1), QC1_samples, NA_character_),
         QC2 = if_else(!is.na(QC2), QC2_samples, NA_character_),
         QC3 = if_else(!is.na(QC3), QC3_samples, NA_character_),
         QC4 = if_else(!is.na(QC4), QC4_samples, NA_character_),
         QC5 = if_else(!is.na(QC5), QC5_samples, NA_character_),
         QC6 = if_else(!is.na(QC6), QC6_samples, NA_character_),
  ) 

# generate sankey data and plot
sankey_data <- sankey_per_sample %>%
  make_long(preQC, QC1, QC2, QC3, QC4, QC5, QC6) %>%
  mutate(tumor_type = if_else(str_detect(node, "Astro"), "Astro.", "Oligo.")) %>%
  filter(!is.na(node))

ggplot(sankey_data, aes(x = x, next_x = next_x, node = node, next_node = next_node,
                        fill = factor(tumor_type), color = factor(tumor_type), label = node)) +
  geom_sankey(flow.alpha = 0.35) +
  geom_sankey_label(size = 3.5, color = 1, fill = "white", alpha = 0.9) +
  theme_sankey(base_size = 18) +
  scale_fill_manual(values = c("Oligo." = "#298C8C", "Astro." = "#800074")) +
  scale_color_manual(values = c("Oligo." = "#298C8C", "Astro." = "#800074")) +
  labs(x = "", fill = "Tumor type", color = "Tumor type") +
  theme(legend.position = "none")

table(sankey_per_sample$QC6, sankey_per_sample$aliquot_analysis_type, useNA= "ifany")