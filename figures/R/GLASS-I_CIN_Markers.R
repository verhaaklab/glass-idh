# GLASS-I: Chromosomal Instability Markers Figure
# Author: E. Kocakavuk
# Date: 2026-06-18

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(tidyverse)
library(survival)
library(survminer)
library(magrittr)
library(purrr)
library(stringr)
library(fs)
library(cowplot)
library(EnvStats)
library(httpgd)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggridges)
library(broom)
library(ggh4x)
library(cowplot)
library(ggpubr)
library(ggalluvial)
library(RPostgres)

hgd()
hgd_browse()

# Read in clinical and annotation data
con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")                  ## data available on Synapse (Synapse ID: syn17038081)

silver_set <- dbGetQuery(con, "SELECT * FROM analysis.silver_set")                     
gold_set <- dbGetQuery(con, "SELECT * FROM analysis.gold_set")                  
mut_freq <- dbGetQuery(con, "SELECT * FROM analysis_mut_freq.csv")
pairs <- dbGetQuery(con, "SELECT * FROM analysis_pairs.csv")
tumor_mut_comparison_anno <- dbGetQuery(con, "SELECT * FROM analysis_tumor_mut_comparison_anno.csv")
aneuploidy <- dbGetQuery(con, "SELECT * FROM analysis_gatk_aneuploidy.csv")
clinical_cases <- dbGetQuery(con, "SELECT * FROM clinical_cases.csv")
clinical_surgeries <- dbGetQuery(con, "SELECT * FROM clinical_surgeries.csv")

dbDisconnect(con)

# Read in copy number signatures & ascat data
cnsig_query <- read_tsv("~/glass5_cnsig_20260206.tsv") %>%                        ## created with cnsig_expl.R + sigpprofilertoolki.smk + sigp_matgen_cnv_prepare.R
                rename(aliquot_barcode = Samples)
ascat_query <- read_tsv("~/glass5_ascat_qc_20260206.tsv") %>%                     ## created with cnsig_expl.R + ascat.smk + ascat_merge.R + ascat_analysis.R
              rename(pair_barcode = aliquot_barcode)
ampsuite_query <- read_tsv("~/glass5_ampsuite_20260317.tsv")                      ## created with GLASS-I_Amplicon_Architect.R
seqz_wgd_query <- read_csv("~/WGD_sequenza_all.csv")                              ## created with WGD_sequenza_script.R
shatterseek <- read_csv("~/chromothripsis_shatterseek_samples.csv")               ## created with GLASS-I_Chromothripsis.R
shatterseek_count <- read_csv("~/chromothripsis_shatterseek_all_calls.csv")       ## created with GLASS-I_Chromothripsis.R
treatment_query <- read_csv("~/prior_treatment_gold_set.csv")%>%                  ## created with GLASS-I_Prior_Treatment_Table.R
driver_changes_query <- read_tsv("~/driver_changes.txt")                          ## created with GLASS-I_Drivers_Table.R
kataegis_query <- read_csv("~/glass5_kataegis_count_all_wgs.csv")                 ## created with kataegis_pipeline.sh
seqz_fga_query <- read_csv("~/fga_seqz_gold_IDHmut_update.csv")                   ## created with script_fga_calculation.R

cnsig <- cnsig_query %>% mutate(case_barcode = substr(aliquot_barcode,1,12)) %>% 
    filter(aliquot_barcode %in% c(gold_set$tumor_barcode_a, gold_set$tumor_barcode_b))
ascat <- ascat_query %>% left_join(pairs) %>% select(-pair_barcode, -normal_barcode) %>% rename(aliquot_barcode = tumor_barcode) %>% 
    filter(aliquot_barcode %in% c(gold_set$tumor_barcode_a, gold_set$tumor_barcode_b))

hm <- mut_freq %>% mutate(hm = ifelse(coverage_adj_mut_freq >=10, "HM", "Non-HM")) %>% 
    filter(aliquot_barcode %in% c(gold_set$tumor_barcode_a, gold_set$tumor_barcode_b))

subtype <- tumor_mut_comparison_anno %>% select(case_barcode, idh_codel_subtype)

# --> EK: We can ignore this part as we are not using the aa_query data
#aa <- aa_query %>% rename(aliquot_barcode = tumor_barcode, ecDNA = sample_classification) %>% 
#  filter(aliquot_barcode %in% c(gold_set$tumor_barcode_a, gold_set$tumor_barcode_b))

ampsuite <- ampsuite_query %>% 
  filter(aliquot_barcode %in% c(gold_set$tumor_barcode_a, gold_set$tumor_barcode_b))

kataegis <- kataegis_query %>% mutate(kataegis_event = ifelse(kataegis_event_count == 0, "no", "yes")) %>% 
  filter(aliquot_barcode %in% c(gold_set$tumor_barcode_a, gold_set$tumor_barcode_b))

seqz_fga <- seqz_fga_query %>% filter(aliquot_barcode %in% c(gold_set$tumor_barcode_a, gold_set$tumor_barcode_b))

post_rec <- gold_set %>% mutate(sample_barcode = substr(tumor_barcode_b, 1, 15)) %>% select(aliquot_barcode = tumor_barcode_b, sample_barcode) %>% 
  inner_join(clinical_surgeries) %>% inner_join(clinical_cases) %>% mutate(post_rec_time = as.numeric(case_overall_survival_mo) - as.numeric(surgical_interval_mo))

os_initial <- gold_set %>% mutate(sample_barcode = substr(tumor_barcode_a, 1, 15)) %>% select(aliquot_barcode = tumor_barcode_a, sample_barcode) %>% 
  inner_join(clinical_surgeries) %>% inner_join(clinical_cases) %>% filter(surgery_number == 1) %>% mutate(os_initial_time = as.numeric(case_overall_survival_mo))

hm_anno <- driver_changes_query %>% select(case_barcode, HM_group) %>% distinct()

chromothripsis_count <- shatterseek_count %>% filter(aliquot_barcode %in% c(gold_set$tumor_barcode_a, gold_set$tumor_barcode_b)) %>% 
  group_by(aliquot_barcode) %>% mutate(chromothripsis.any.conf = high.confidence.test1 | high.confidence.test2 | low.confidence.test) %>% 
  summarise(chromothripsis_count = sum(chromothripsis.any.conf), .groups = "drop")

sv_count <- shatterseek_count %>% filter(aliquot_barcode %in% c(gold_set$tumor_barcode_a, gold_set$tumor_barcode_b)) %>% 
  select(aliquot_barcode, sv_count = number_SVs_sample) %>% distinct()

# 1. WGD 
## Use correct way with ASCAT + Sequenza WGD calls
ascat_seqz_joined <-
  ascat %>%
  left_join(seqz_wgd_query, by = "aliquot_barcode") %>%
  mutate(case_barcode = substr(aliquot_barcode, 1, 12)) %>%
  filter(!is.na(WGD), idh_codel_subtype != "IDHwt") %>%
  mutate(type = case_when(aliquot_barcode %in% gold_set$tumor_barcode_a ~ "Initial", aliquot_barcode %in% gold_set$tumor_barcode_b ~ "Recurrence", TRUE ~ "Other"),
    WGD = ifelse(WGD == "0", "Absent", "Present"), WGD_final = ifelse(frac_major_ge2 < 0.5, "Absent", WGD)) %>%
  filter(type != "Other")

one_b_sankey <- ascat_seqz_joined %>%
  filter(!is.na(WGD_final)) %>%
  mutate(case_barcode = substr(aliquot_barcode, 1, 12),
    type = case_when(aliquot_barcode %in% gold_set$tumor_barcode_a ~ "Initial", aliquot_barcode %in% gold_set$tumor_barcode_b ~ "Recurrence", TRUE ~ "Other")) %>%
  filter(idh_codel_subtype != "IDHwt", type != "Other") %>%
  select(case_barcode, type, WGD_final, idh_codel_subtype) %>%
  pivot_wider(names_from = type, values_from = WGD_final, id_cols = c(case_barcode, idh_codel_subtype)) %>%
  filter(!is.na(Initial) & !is.na(Recurrence)) %>%
  count(idh_codel_subtype, Initial, Recurrence, name = "Freq")

totals <- one_b_sankey %>% group_by(idh_codel_subtype) %>% summarise(n = sum(Freq), .groups = "drop")

label_oligo <- paste0("Oligo. (n=", totals$n[totals$idh_codel_subtype == "IDHmut-codel"], ")")
label_astro <- paste0("Astro. (n=", totals$n[totals$idh_codel_subtype == "IDHmut-noncodel"], ")")

one_b_sankey <- one_b_sankey %>%
  left_join(totals, by = "idh_codel_subtype") %>%
  mutate(idh_codel_subtype = case_when(idh_codel_subtype == "IDHmut-noncodel" ~ label_astro, idh_codel_subtype == "IDHmut-codel" ~ label_oligo, TRUE ~ idh_codel_subtype),
    prop = Freq / n) %>% select(-n) %>%
  mutate(idh_codel_subtype = factor(idh_codel_subtype, levels = c(label_oligo, label_astro)))

one_b_mcnemar <- one_b_sankey %>%
  group_by(idh_codel_subtype) %>%
  summarise(Absent_Absent = coalesce(first(Freq[Initial == "Absent" & Recurrence == "Absent"]), 0),
    Absent_Present = coalesce(first(Freq[Initial == "Absent" & Recurrence == "Present"]), 0),
    Present_Absent = coalesce(first(Freq[Initial == "Present" & Recurrence == "Absent"]), 0),
    Present_Present = coalesce(first(Freq[Initial == "Present" & Recurrence == "Present"]), 0),
    .groups = "drop") %>%
  rowwise() %>%
  mutate(mcnemar_p = tryCatch(mcnemar.test(matrix(c(Absent_Absent, Absent_Present, Present_Absent, Present_Present),
        nrow = 2, byrow = TRUE))$p.value, error = function(e) NA_real_), label = paste0("McNemar, p = ", sprintf("%2.1e", mcnemar_p)),
        x = 1.5, y = 0.95) %>% ungroup()

panel_a <-  one_b_sankey %>%
  mutate(ever_present = ifelse(Initial == "Present" | Recurrence == "Present", "Present", "Absent")) %>%
  mutate(WGD_fill = case_when(ever_present == "Absent" ~ "Absent", 
    ever_present == "Present" & idh_codel_subtype == label_oligo ~ "Present - Oligo", 
    ever_present == "Present" & idh_codel_subtype == label_astro ~ "Present - Astro")) %>%
  ggplot(aes(axis1 = Initial, axis2 = Recurrence, y = prop, fill = WGD_fill)) +
  geom_alluvium(width = 0.25, knot.pos = 0.3, alpha = 0.8) +
  geom_stratum(width = 0.25, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  geom_text(data = one_b_mcnemar, aes(x = x, y = 1.025, label = label), inherit.aes = FALSE, size = 4) +
  scale_x_discrete(limits = c("Initial", "Recurrence"), expand = c(0.15, 0.05)) +
  scale_fill_manual(values = c( "Absent" = "#999999", "Present - Oligo" = "#298C8C", "Present - Astro" = "#800074"),
  breaks = c("Absent", "Present - Oligo", "Present - Astro"), labels = c("Absent", "Present, oligo.", "Present, astro."), name = "WGD ever") +
  ggh4x::facet_grid2(~ idh_codel_subtype, strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = c("#298C8C66", "#80007466")),
      text_x = ggh4x::elem_list_text(face = c("bold", "bold")))) +
  theme_bw() + labs(x = "WGD", y = "Proportion") + theme(panel.grid = element_blank(), legend.position = "none")

panel_a

# Chromothripsis higher in astro vs oligo at initial and recurrence 
## 2a. Chromothripsis enriched in Astros
pval_chrom_initial <- shatterseek %>% filter(!is.na(chromothripsis.any.conf)) %>% 
  mutate(case_barcode = substr(aliquot_barcode, 1, 12), 
        type = ifelse(aliquot_barcode %in% gold_set$tumor_barcode_a, "Initial", 
                              ifelse(aliquot_barcode %in% gold_set$tumor_barcode_b, "Recurrence", "Other"))) %>% 
  left_join(subtype) %>% filter(idh_codel_subtype != "IDHwt") %>% 
  filter(type == "Initial") %>% mutate(Chromothripsis = ifelse(chromothripsis.any.conf == "0", "No", "Yes")) %$% 
  table(idh_codel_subtype, Chromothripsis) %>% fisher.test() %>% .$p.value

pval_chrom_recurrence <- shatterseek %>% filter(!is.na(chromothripsis.any.conf)) %>% 
  mutate(case_barcode = substr(aliquot_barcode, 1, 12), 
        type = ifelse(aliquot_barcode %in% gold_set$tumor_barcode_a, "Initial", 
                              ifelse(aliquot_barcode %in% gold_set$tumor_barcode_b, "Recurrence", "Other"))) %>% 
  left_join(subtype) %>% filter(idh_codel_subtype != "IDHwt") %>% 
  filter(type == "Recurrence") %>% mutate(Chromothripsis = ifelse(chromothripsis.any.conf == "0", "No", "Yes")) %$% 
  table(idh_codel_subtype, Chromothripsis) %>% fisher.test() %>% .$p.value

pval_chrom_df <- tibble(type = factor(c("Initial", "Recurrence"), levels = c("Initial", "Recurrence")),
    Chromothripsis = c("Yes", "Yes"), idh_codel_subtype = c("Oligo", "Oligo"), x = 1.5, y = 0.55,
    label = c(paste0("Fisher, p = ", sprintf("%2.1e", pval_chrom_initial)), paste0("Fisher, p = ", sprintf("%2.1e", pval_chrom_recurrence))))

panel_b <- shatterseek %>% filter(!is.na(chromothripsis.any.conf)) %>% 
  mutate(case_barcode = substr(aliquot_barcode, 1, 12), 
    type = ifelse(aliquot_barcode %in% gold_set$tumor_barcode_a, "Initial", 
            ifelse(aliquot_barcode %in% gold_set$tumor_barcode_b, "Recurrence", "Other"))) %>% 
  left_join(subtype) %>% filter(idh_codel_subtype != "IDHwt") %>% 
  mutate(idh_codel_subtype = ifelse(idh_codel_subtype == "IDHmut-noncodel", "Astro", "Oligo"),
    idh_codel_subtype = factor(idh_codel_subtype, levels = c("Oligo", "Astro")),
    Chromothripsis = ifelse(chromothripsis.any.conf == "0", "No", "Yes")) %$% 
  table(type, Chromothripsis, idh_codel_subtype) %>% as.data.frame() %>%  group_by(idh_codel_subtype, type) %>%
  mutate(prop = 100 * Freq / sum(Freq), chrom_fill = case_when(Chromothripsis == "No" ~ "No",
      Chromothripsis == "Yes" & idh_codel_subtype == "Oligo" ~ "Yes - Oligo", Chromothripsis == "Yes" & idh_codel_subtype == "Astro" ~ "Yes - Astro")) %>%
  ungroup() %>% ggplot(aes(x = idh_codel_subtype, y = Freq, fill = chrom_fill)) +
  geom_bar(stat = "identity", position = "fill") + theme_bw() + 
  scale_fill_manual(values = c("No" = "#999999", "Yes - Oligo" = "#298C8C", "Yes - Astro" = "#800074")) + 
  geom_text(aes(label = paste0(Freq, " ", "(", round(prop), "%)")), position = position_fill(vjust = 0.5), color = "white") +  
  geom_text(data = pval_chrom_df, aes(x = x, y = 1.025, label = label), inherit.aes = FALSE,
    vjust = 0, hjust = 0.5, size = 4, parse = FALSE) + xlab("Chromothripsis") + ylab("Proportion") +
  labs(fill = "Chromothripsis") + facet_grid(~ type) +
  theme(legend.position = "none", strip.background = element_rect(fill = "white"))

panel_b

## 2b. Chromothripsis enriched in RT-treated Astros
tab_chr_rt <- shatterseek %>% 
  filter(aliquot_barcode %in% gold_set$tumor_barcode_b) %>% rename(Chromothripsis = chromothripsis.any.conf) %>% 
  inner_join(treatment_query) %>%  
  mutate(Chromothripsis = ifelse(Chromothripsis == 1, "Present", "Absent"), 
        RT = ifelse(RT == 1, "Yes", "No")) %$% table(RT, Chromothripsis, idh_codel_subtype)
  
# Separate p-values for each subtype
pval_chr_rt_oligo <- shatterseek %>% filter(aliquot_barcode %in% gold_set$tumor_barcode_b) %>% 
  filter(idh_codel_subtype == "IDHmut-codel") %>% rename(Chromothripsis = chromothripsis.any.conf) %>% 
  inner_join(treatment_query) %>% 
  mutate(Chromothripsis = ifelse(Chromothripsis == 1, "Present", "Absent"), 
        RT = ifelse(RT == 1, "Yes", "No")) %$% table(RT, Chromothripsis) %>% fisher.test() %>% .$p.value

pval_chr_rt_astro <- shatterseek %>% 
  filter(aliquot_barcode %in% gold_set$tumor_barcode_b) %>% filter(idh_codel_subtype == "IDHmut-noncodel") %>% 
  rename(Chromothripsis = chromothripsis.any.conf) %>% inner_join(treatment_query) %>% 
  mutate(Chromothripsis = ifelse(Chromothripsis == 1, "Present", "Absent"), 
        RT = ifelse(RT == 1, "Yes", "No")) %$% table(RT, Chromothripsis) %>% fisher.test() %>% .$p.value

pval_chr_rt_df <- tibble(idh_codel_subtype = factor(c("Oligo. (n=51)", "Astro. (n=38)"), levels = c("Oligo. (n=51)", "Astro. (n=38)")),
  Chromothripsis = c("Present", "Present"), RT = c("No", "No"), x = 1.5, y = 0.55,
  label = c(paste0("Fisher, p = ", sprintf("%2.1e", pval_chr_rt_oligo)), paste0("Fisher, p = ", sprintf("%2.1e", pval_chr_rt_astro))))

panel_c <- tab_chr_rt %>% as.data.frame() %>% 
  mutate(idh_codel_subtype = ifelse(idh_codel_subtype == "IDHmut-noncodel", "Astro. (n=38)", "Oligo. (n=51)"), 
        idh_codel_subtype = factor(idh_codel_subtype, levels=c("Oligo. (n=51)", "Astro. (n=38)"))) %>% 
  ggplot(aes(Chromothripsis, RT, fill = Freq)) + geom_tile(color = "black") +
  geom_text(aes(label = Freq), size = 6) +  geom_text(data = pval_chr_rt_df, aes(x = x, y = y+1.95, label = label), inherit.aes = FALSE, vjust = 0, hjust = 0.5, size = 4, parse = FALSE) +
  scale_fill_gradient(low = "white", high = "grey30") + theme_bw() + labs(fill = "Count", y = "Prior Radiotherapy") + 
  ggh4x::facet_grid2(~ idh_codel_subtype, strip = ggh4x::strip_themed(
  background_x = ggh4x::elem_list_rect(fill = c("#298C8C66", "#80007466")), text_x = ggh4x::elem_list_text(face = c("bold", "bold")))) + 
  theme(legend.position = "none")

panel_c

### 3. Kataegis higher in astro vs oligo at initial and recurrence
kataegis_plot_df <- kataegis %>% mutate(case_barcode = substr(aliquot_barcode, 1, 12)) %>% 
  left_join(subtype, by = "case_barcode") %>% filter(idh_codel_subtype != "IDHwt") %>% 
  left_join(hm_anno) %>% filter(HM_group == "NHM") %>%
  mutate(type = ifelse(aliquot_barcode %in% gold_set$tumor_barcode_a, "Initial", 
      ifelse(aliquot_barcode %in% gold_set$tumor_barcode_b, "Recurrence", "Other"))) %$%
  table(kataegis_event, idh_codel_subtype, type) %>% as.data.frame() %>% 
  mutate(kataegis_event = ifelse(kataegis_event == "yes", "Kataegis", "Non-Kataegis"),
    kataegis_event = factor(kataegis_event, levels = c("Non-Kataegis", "Kataegis")),
    type = factor(type, levels = c("Initial", "Recurrence", "Other"))) %>%
  group_by(idh_codel_subtype, type) %>% mutate(prop = 100 * Freq / sum(Freq)) %>% ungroup()

# Fisher's exact test per type
fisher_kataegis_df <- kataegis_plot_df %>% filter(type != "Other") %>% group_by(type) %>%
  summarise(p_value = fisher.test(xtabs(Freq ~ kataegis_event + idh_codel_subtype, data = cur_data()))$p.value, .groups = "drop") %>%
  mutate(p_label = paste0("Fisher, p = ", sprintf("%2.1e", p_value)), x = 1.5, y = 1)

# Plot with Fisher p-values
panel_d <- kataegis_plot_df %>% filter(type != "Other") %>%
  mutate(idh_codel_subtype = factor(idh_codel_subtype, levels = c("IDHmut-codel", "IDHmut-noncodel"), labels = c("Oligo", "Astro")), 
        kataegis_fill = case_when(kataegis_event == "Non-Kataegis" ~ "Non-Kataegis", kataegis_event == "Kataegis" & idh_codel_subtype == "Oligo" ~ "Kataegis - Oligo", kataegis_event == "Kataegis" & idh_codel_subtype == "Astro" ~ "Kataegis - Astro"), 
        kataegis_fill = factor(kataegis_fill, levels = c("Non-Kataegis", "Kataegis - Oligo", "Kataegis - Astro"))) %>%
  ggplot(aes(x = idh_codel_subtype, y = Freq, fill = kataegis_fill)) + geom_bar(stat = "identity", position = "fill") + 
  geom_text(aes(label = paste0(Freq, " (", round(prop), "%)")), position = position_fill(vjust = 0.5), color = "white") +
  geom_text(data = fisher_kataegis_df, aes(x = x, y = 1.025, label = p_label), inherit.aes = FALSE, size = 4) +
  facet_grid(~type) + scale_fill_manual(values = c("Non-Kataegis" = "#999999", "Kataegis - Oligo" = "#298C8C", "Kataegis - Astro" = "#800074")) +
  theme_bw() +  theme(strip.background = element_rect(fill = "white")) + 
  theme(legend.position = "none") + xlab("Kataegis") + ylab("Proportion")

panel_d

## 4. Chromothripsis higher in kataegis vs non-kataegis at recurrence in non-HM cases
panel_e_plot_df <- kataegis %>% mutate(case_barcode = substr(aliquot_barcode, 1, 12)) %>% 
  left_join(subtype) %>% filter(idh_codel_subtype != "IDHwt", aliquot_barcode %in% gold_set$tumor_barcode_b) %>% 
  left_join(shatterseek) %>% left_join(hm_anno) %>% filter(HM_group == "NHM") %$%
  table(kataegis_event, chromothripsis.any.conf) %>% as.data.frame() %>% 
  mutate(kataegis_event = ifelse(kataegis_event == "yes", "Kataegis", "Non-Kataegis"),
    kataegis_event = factor(kataegis_event, levels = c("Non-Kataegis", "Kataegis")),
    Chromothripsis = ifelse(chromothripsis.any.conf == "1", "Present", "Absent"),
    Chromothripsis = factor(Chromothripsis, levels = c("Absent", "Present")),
    type = factor("Recurrence", levels = "Recurrence")) %>%
  group_by(kataegis_event) %>% mutate(prop = 100 * Freq / sum(Freq)) %>% ungroup()

fisher_panel_e_df <- panel_e_plot_df %>%
  summarise(p_value = fisher.test(xtabs(Freq ~ kataegis_event + Chromothripsis, data = cur_data()))$p.value, .groups = "drop") %>%
  mutate(p_label = paste0("Fisher, p = ", sprintf("%2.1e", p_value)), x = 1.5, y = 1)

panel_e <- panel_e_plot_df %>%
  ggplot(aes(x = kataegis_event, y = Freq, fill = Chromothripsis)) + geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = paste0(Freq, " (", round(prop), "%)")), position = position_fill(vjust = 0.5), color = "white") +
  geom_text(data = fisher_panel_e_df, aes(x = x, y = 1.025, label = p_label), inherit.aes = FALSE, size = 4) +
  facet_grid(~ type) + scale_fill_manual(values = c("Absent" = "#999999", "Present" = "#377EB8")) +
  theme_bw() + theme(strip.background = element_rect(fill = "white"), legend.position = "none") +
  xlab("Chromothripsis") + ylab("Proportion")

panel_e

#5. Hypermutation is assciated with chromothripsis-signature CN7
panel_f <- cnsig %>% 
  tidyr::pivot_longer(cols = -c(aliquot_barcode, case_barcode), names_to = "CN", values_to = "value") %>%
  group_by(aliquot_barcode) %>% mutate(prop = value / sum(value, na.rm = TRUE)) %>%
  ungroup() %>% inner_join(subtype) %>% filter(idh_codel_subtype != "IDHwt") %>%
  left_join(hm_anno) %>% filter(aliquot_barcode %in% gold_set$tumor_barcode_b) %>%
  filter(!is.na(HM_group)) %>% group_by(CN) %>%
  summarise(mean_HM  = mean(prop[HM_group == "HM"], na.rm = TRUE), mean_NHM = mean(prop[HM_group == "NHM"], na.rm = TRUE),
    p_value = wilcox.test(prop ~ HM_group)$p.value, .groups = "drop") %>%
  mutate(sig = ifelse(p_value <= 0.05, "yes", "no"), number = CN, type = factor("Recurrence")) %>% filter(mean_HM!=0, mean_NHM!=0) %>% 
  ggplot(aes(x=log2((mean_HM)/(mean_NHM)), y=-log10(p_value))) + geom_point(aes(color = sig), size = 4, alpha=0.5) + 
  ggrepel::geom_text_repel(aes(label = number), size = 3.5, nudge_y = 0.1, max.overlaps = Inf, box.padding = 0.05, point.padding = 0.05, min.segment.length = 0, segment.color = "grey50", segment.size = 0.3, direction = "both", force_pull=5, seed = 1) +
  geom_hline(yintercept = -log10(0.05), color = "grey30") + geom_vline(xintercept = 0, color = "grey50", linetype="longdash") +
  annotate("text", x = 0.5, y = 1.5, label = "Enriched in\nHypermutant", color = "grey50") +
  annotate("text", x = -0.6, y = 1.5, label = "Enriched in\nNon-Hypermutant", color = "grey50") + facet_grid(~ type) +
  theme_classic() + scale_color_manual(values = c("black", "red")) + theme(legend.position="none", strip.background = element_rect(fill = "white")) + 
  ylab(expression(-log[10](p-value))) + xlab(expression(log[2](FC)))

panel_f

## Arrange panels into main figure
pdf("~/Desktop/GLASS_20260206/plots/plot_paper_20260517.pdf", width = 12, height = 12, useDingbats = FALSE)
  ggarrange(panel_a, panel_b, panel_c, panel_d, panel_e, panel_f, ncol = 2, nrow = 3, labels = "AUTO", heights = c(1, 1, 1))
dev.off()

## Extended data figures
#1. WGD enriched for CN2/CN20

plot_cn2_df <- cnsig %>% 
  tidyr::pivot_longer(cols = -c(aliquot_barcode, case_barcode), names_to = "CN", values_to = "value") %>%
  group_by(aliquot_barcode) %>% mutate(prop = value / sum(value, na.rm = TRUE)) %>% ungroup() %>%
  inner_join(subtype) %>% filter(idh_codel_subtype != "IDHwt") %>% inner_join(ascat_seqz_joined) %>%
  filter(!is.na(WGD_final)) %>% group_by(WGD_final, CN) %>%
  summarise(frac_with_CN = mean(prop > 0, na.rm = TRUE), median_prop = median(prop[prop > 0], na.rm = TRUE), .groups = "drop")

edf_panel_a <- plot_cn2_df %>% filter(CN %in% c("CN2", "CN20")) %>% 
  ggplot(aes(x = WGD_final, y = CN)) +
  geom_point(aes(size = frac_with_CN, fill = median_prop), shape = 21, color = "black") +
  #geom_text(data = pval_cn_df, aes(x = x, y = y, label = p_label), inherit.aes = FALSE, size = 3) +
  scale_size(range = c(0, 30), labels = scales::percent) +
  scale_fill_viridis_c(na.value = "grey90") +  theme_bw() +
  labs(x = "WGD", y = NULL, size = "Fraction with any contribution", fill = "Median contribution") +
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 7), legend.key.size = unit(0.35, "cm"))

edf_panel_a

#2. Chromothripsis enriched for CN4-6 and CN7-8
plot_cn48_df <- cnsig %>% 
  tidyr::pivot_longer(cols = -c(aliquot_barcode, case_barcode), names_to = "CN", values_to = "value") %>%
  group_by(aliquot_barcode) %>% mutate(prop = value / sum(value, na.rm = TRUE)) %>%
  ungroup() %>% filter(CN %in% c("CN4", "CN5", "CN6", "CN7", "CN8")) %>%
  mutate(CN_group = case_when(CN %in% c("CN4", "CN5", "CN6") ~ "CN4-6", CN %in% c("CN7", "CN8") ~ "CN7-8")) %>%
  group_by(aliquot_barcode, case_barcode, CN_group) %>% summarise(prop = sum(prop, na.rm = TRUE), .groups = "drop") %>%
  inner_join(subtype) %>% filter(idh_codel_subtype != "IDHwt") %>% inner_join(shatterseek) %>% filter(!is.na(chromothripsis.any.conf)) %>%
  mutate(Chromothripsis = ifelse(chromothripsis.any.conf == 0, "Absent", "Present")) %>% group_by(Chromothripsis, CN_group) %>%
  summarise(frac_with_CN = mean(prop > 0, na.rm = TRUE), median_prop = median(prop[prop > 0], na.rm = TRUE),.groups = "drop")

edf_panel_b <- plot_cn48_df  %>% filter(CN_group %in% c("CN4-6", "CN7-8")) %>% 
  ggplot(aes(x = Chromothripsis, y = CN_group)) +
  geom_point(aes(size = frac_with_CN, fill = median_prop), shape = 21, color = "black") +
  scale_size(range = c(1, 20), labels = scales::percent) + scale_fill_viridis_c(na.value = "grey90") +
  theme_bw() + labs(x = "Chromothripsis", y = NULL, size = "Fraction with any contribution", fill = "Median contribution") + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 7), legend.key.size = unit(0.35, "cm"))

edf_panel_b

##### Does WGD correlate with FGA? 
edf_panel_c <- ascat_seqz_joined %>% 
  left_join(seqz_fga) %>%  mutate(case_barcode = substr(aliquot_barcode, 1, 12)) %>% 
  left_join(subtype) %>% filter(idh_codel_subtype != "IDHwt") %>% mutate(WGD = ifelse(WGD == "0", "Absent", "Present")) %>% 
  mutate(WGD_final = ifelse(frac_major_ge2 <0.5, "Absent", WGD)) %>% mutate(WGD_final = factor(WGD_final, levels = c("Present", "Absent"))) %>%
  mutate(timepoint = ifelse(aliquot_barcode %in% gold_set$tumor_barcode_a, "Initial", "Recurrence")) %>% 
  ggplot(aes(x=WGD_final, y=FGA)) + geom_boxplot(aes(fill=WGD_final), linewidth = 0.7) +
  theme_bw() + stat_compare_means(label.x = 1.1, label.y = 1.05,
                     aes(label = sprintf("Wilcoxon, p = %2.1e", as.numeric(..p.format..)))) +
  labs(y="FGA", x="WGD status") + scale_fill_manual(values = c("Present" = "#377EB8", "Absent" = "#999999")) + 
  stat_summary(fun.data = function(x) {data.frame(y = -0.025,label = paste0("italic(n) == ", length(x)))}, geom = "text",parse = TRUE) +
  facet_grid(~timepoint) + theme(legend.position="none", strip.background = element_rect(fill = "white")) 

edf_panel_c

##### Does chromothripsis correlate with FGA? 
edf_panel_d <- 
  shatterseek %>% filter(!is.na(chromothripsis.any.conf)) %>% 
  left_join(seqz_fga) %>%  mutate(case_barcode = substr(aliquot_barcode, 1, 12)) %>% 
  left_join(subtype) %>% filter(idh_codel_subtype != "IDHwt") %>% 
  mutate(chromothripsis = ifelse(chromothripsis.any.conf == 0, "Absent", "Present"), chromothripsis = factor(chromothripsis, levels = c("Present", "Absent")), 
         timepoint = ifelse(aliquot_barcode %in% gold_set$tumor_barcode_a, "Initial", "Recurrence")) %>% 
  ggplot(aes(x=chromothripsis, y=FGA)) + geom_boxplot(aes(fill=chromothripsis), linewidth = 0.7) +
  theme_bw() + stat_compare_means(label.x = 1.1, label.y = 1.05,
                     aes(label = sprintf("Wilcoxon, p = %2.1e", as.numeric(..p.format..)))) +
  labs(y="FGA", x="Chromothripsis") + scale_fill_manual(values = c("Present" = "#377EB8", "Absent" = "#999999")) + 
  stat_summary(fun.data = function(x) {data.frame(y = -0.025,label = paste0("italic(n) == ", length(x)))}, geom = "text",parse = TRUE) +
  facet_grid(~timepoint) + theme(legend.position="none", strip.background = element_rect(fill = "white")) 

edf_panel_d

##### Does Hypermutation correlate with FGA? 
edf_panel_e <- 
seqz_fga %>% mutate(case_barcode = substr(aliquot_barcode, 1, 12)) %>% 
  left_join(hm_anno) %>% left_join(subtype) %>% filter(idh_codel_subtype != "IDHwt") %>% filter(!is.na(HM_group)) %>%
  mutate(timepoint = ifelse(aliquot_barcode %in% gold_set$tumor_barcode_a, "Initial", "Recurrence"), 
          hypermutation = ifelse(HM_group == "HM", "Present", "Absent"), hypermutation = factor(hypermutation, levels = c("Present", "Absent"))) %>% 
  ggplot(aes(x=hypermutation, y=FGA)) + geom_boxplot(aes(fill=hypermutation), linewidth = 0.7) +
  theme_bw() + stat_compare_means(label.x = 1.1, label.y = 1.05,
                     aes(label = sprintf("Wilcoxon, p = %2.1e", as.numeric(..p.format..)))) +
  labs(y="FGA", x="Hypermutation") + scale_fill_manual(values = c("Present" = "#377EB8", "Absent" = "#999999")) +
  stat_summary(fun.data = function(x) {data.frame(y = -0.025,label = paste0("italic(n) == ", length(x)))}, geom = "text",parse = TRUE) +
  facet_grid(~timepoint) + theme(legend.position="none", strip.background = element_rect(fill = "white")) 

edf_panel_e

pdf("~/plot_suppl_20260517.pdf", width = 14, height = 10, useDingbats = FALSE)
ggarrange(
  ggarrange(edf_panel_a, edf_panel_b, nrow=1, labels = c("A", "B")), 
  ggarrange(edf_panel_c, edf_panel_d, edf_panel_e, nrow=1, labels = c("C", "D", "E")), 
  nrow = 2, heights = c(1, 0.6))
dev.off()
