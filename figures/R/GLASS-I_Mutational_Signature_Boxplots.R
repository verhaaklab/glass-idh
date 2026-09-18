# GLASS-I: Treatment-Associated Mutational Signatures Box Plots
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(tidyverse)
library(rstatix)

# load data
mutational.signature.cases <- read_csv("~/summary_molecular_cohort.csv")%>%     ## created with GLASS-I_Summarizing_Table_Molecular_Cohort.R
  filter(sequencing_platform == "WGS")%>%
  filter(hypermutation_status == "hypermutant" | (hypermutation_status == "non-hypermutant" & initial_TMB < 10))%>%
  select(case_barcode, glioma_type, hypermutation_status, radiotherapy_prior_to_recurrence, temozolomide_prior_to_recurrence)

skeleton <- mutational.signature.cases %>%
  crossing(
    mut_sig = c("SBS11", "SBS119", "ID8"),
    variant_state = c("P", "R", "S")
  )

mutational.signature.data <- read.csv("~/mut_sig_prop.csv") %>%                 ## created with mutational_signature_pipeline.sh
  filter(mut_sig == "SBS11" | mut_sig == "SBS119" | mut_sig == "ID8") %>%
  right_join(skeleton, by = c("case_barcode", "mut_sig", "variant_state")) %>%
  mutate(sig_prop = replace_na(sig_prop, 0))

sbs119.nhm.tmz.treated <- mutational.signature.data%>%
  filter(mut_sig == "SBS119" & hypermutation_status == "non-hypermutant" & temozolomide_prior_to_recurrence == "treated")

kw_stats_sbs119 <- sbs119.nhm.tmz.treated %>%
  group_by(glioma_type) %>%
  kruskal_test(sig_prop ~ variant_state) %>%
  mutate(label = paste0("Kruskal-Wallis, p = ", signif(p, 2)))

ggplot(data = sbs119.nhm.tmz.treated)+
  geom_boxplot(aes(x=variant_state, y=sig_prop, fill = variant_state))+
  geom_text(
    data = kw_stats_sbs119,
    aes(x = 2, y = 1, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_classic()+
  facet_grid(~glioma_type)+
  scale_fill_manual(values = c("P" = "#CA2F66","S" = "#CA932F", "R" = "#2fb3ca"))+
  ylab("SBS119 proportion (%)")+
  xlab("Signature timing")+
  ylim(c(0,1))
  
sbs11.hm.tmz.treated <- mutational.signature.data%>%
  filter(mut_sig == "SBS11" & hypermutation_status == "hypermutant" & temozolomide_prior_to_recurrence == "treated")

kw_stats_sbs11 <- sbs11.hm.tmz.treated %>%
  group_by(glioma_type) %>%
  kruskal_test(sig_prop ~ variant_state) %>%
  mutate(label = paste0("Kruskal-Wallis, p = ", signif(p, 2)))

ggplot(data = sbs11.hm.tmz.treated)+
  geom_boxplot(aes(x=variant_state, y=sig_prop, fill = variant_state))+
  geom_text(
    data = kw_stats_sbs11,
    aes(x = 2, y = 1, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_classic()+
  facet_grid(~glioma_type)+
  scale_fill_manual(values = c("P" = "#CA2F66","S" = "#CA932F", "R" = "#2fb3ca"))+
  ylab("SBS11 proportion (%)")+
  xlab("Signature timing")+
  ylim(c(0,1))

id8.nhm.rt.treated <- mutational.signature.data%>%
  filter(mut_sig == "ID8" & hypermutation_status == "non-hypermutant" & radiotherapy_prior_to_recurrence == "treated")

kw_stats_id8 <- id8.nhm.rt.treated %>%
  group_by(glioma_type) %>%
  kruskal_test(sig_prop ~ variant_state) %>%
  mutate(label = paste0("Kruskal-Wallis, p = ", signif(p, 2)))

ggplot(data = id8.nhm.rt.treated)+
  geom_boxplot(aes(x=variant_state, y=sig_prop, fill = variant_state))+
  geom_text(
    data = kw_stats_id8,
    aes(x = 2, y = 1, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_classic()+
  facet_grid(~glioma_type)+
  scale_fill_manual(values = c("P" = "#CA2F66","S" = "#CA932F", "R" = "#2fb3ca"))+
  ylab("ID8 proportion (%)")+
  xlab("Signature timing")+
  ylim(c(0,1))

# Wilcoxon test for SBS119 comparing recurrence-fraction between oligodendrogliomas and astrocytomas
sbs119.nhm.tmz.treated %>%
filter(variant_state == "R") %>%
wilcox_test(sig_prop ~ glioma_type, paired = FALSE)