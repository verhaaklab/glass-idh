# GLASS-I: McNemar's Tests CNAs
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# Load required libraries
library(tidyverse)

# Load and preprocess data
real.primaries <- read.delim("~/fga_calc_annotated.txt")%>%                   ## created with script_fga_calculation.R
  select(case_barcode)%>%
  distinct(.)

cnv.data <- read.delim("~/driver_changes.txt")%>%                             ## created with GLASS-I_Drivers_Table.R
  filter(gene_symbol %in% c("CCND2", "CDK4", "CDK6", "EGFR", "MET", "MYC", "MYCN", "PDGFRA", "PTPRZ1", "Any Amplification",
                            "ATRX", "CDKN2A", "CDKN2B", "PTEN", "Any HD"))%>%
  semi_join(real.primaries, by = "case_barcode")
           
# Run McNemar's test for CNVs in all hypermutant oligos
cnv.data.oligo.hm <- cnv.data %>%
  filter(idh_codel_subtype == "IDHmut-codel") %>%
  filter(HM_group == "HM")


results.oligo.hm <- tibble(Gene = character(), pvalue = numeric())

for (gene in unique(cnv.data.oligo.hm$gene_symbol)) {
  test_data <- cnv.data.oligo.hm %>% filter(gene_symbol == gene) 
  gene <- unique(test_data$gene_symbol)
  mcnemar.matrix <-
    matrix(c(nrow(test_data %>% filter(driver_change == "S")), 
             nrow(test_data %>% filter(driver_change == "R")),     
             nrow(test_data %>% filter(driver_change == "P")), 
             nrow(test_data %>% filter(is.na(driver_change )))
    ),
    nrow = 2,
    dimnames = list("Primary" = c("AMP/HD", "No_amp/hd"),
                    "Recurrence" = c("AMP/HD", "No_amp/hd")))
  pvalue <- mcnemar.test(mcnemar.matrix)$p.value
  results.oligo.hm <- results.oligo.hm %>% add_row(Gene = as.character(gene), pvalue = pvalue)
}

oligo.hm.single.amp <- results.oligo.hm%>%
  filter(Gene == "CCND2" | Gene == "CDK4" | Gene == "CDK6" | Gene == "EGFR" | Gene == "MET" |Gene == "MYC" |Gene == "MYCN" |Gene == "PDGFRA" |Gene == "PTPRZ1")%>%
  mutate(p.adjust = p.adjust(as.numeric(pvalue), method = "fdr"))

print(oligo.hm.single.amp)

oligo.hm.single.hd <- results.oligo.hm%>%
  filter(Gene == "ATRX" | Gene == "CDKN2A" | Gene == "CDKN2B" | Gene == "PTEN")%>%
  mutate(p.adjust = p.adjust(as.numeric(pvalue), method = "fdr"))

print(oligo.hm.single.hd)

oligo.hm.panel <- results.oligo.hm%>%
  filter(Gene == "Any Amplification" | Gene == "Any HD")%>%
  mutate(p.adjust = p.adjust(as.numeric(pvalue), method = "fdr"))

print(oligo.hm.panel)

# Run McNemar's test for CNVs in all non-hypermutant oligos
cnv.data.oligo.nhm <- cnv.data %>%
  filter(idh_codel_subtype == "IDHmut-codel") %>%
  filter(HM_group == "NHM")


results.oligo.nhm <- tibble(Gene = character(), pvalue = numeric())

for (gene in unique(cnv.data.oligo.nhm$gene_symbol)) {
  test_data <- cnv.data.oligo.nhm %>% filter(gene_symbol == gene) 
  gene <- unique(test_data$gene_symbol)
  mcnemar.matrix <-
    matrix(c(nrow(test_data %>% filter(driver_change == "S")), 
             nrow(test_data %>% filter(driver_change == "R")),     
             nrow(test_data %>% filter(driver_change == "P")), 
             nrow(test_data %>% filter(is.na(driver_change )))
    ),
    nrow = 2,
    dimnames = list("Primary" = c("AMP/HD", "No_amp/hd"),
                    "Recurrence" = c("AMP/HD", "No_amp/hd")))
  pvalue <- mcnemar.test(mcnemar.matrix)$p.value
  results.oligo.nhm <- results.oligo.nhm %>% add_row(Gene = as.character(gene), pvalue = pvalue)
}

oligo.nhm.single.amp <- results.oligo.nhm%>%
  filter(Gene == "CCND2" | Gene == "CDK4" | Gene == "CDK6" | Gene == "EGFR" | Gene == "MET" |Gene == "MYC" |Gene == "MYCN" |Gene == "PDGFRA" |Gene == "PTPRZ1")%>%
  mutate(p.adjust = p.adjust(as.numeric(pvalue), method = "fdr"))

print(oligo.nhm.single.amp)

oligo.nhm.single.hd <- results.oligo.nhm%>%
  filter(Gene == "ATRX" | Gene == "CDKN2A" | Gene == "CDKN2B" | Gene == "PTEN")%>%
  mutate(p.adjust = p.adjust(as.numeric(pvalue), method = "fdr"))

print(oligo.nhm.single.hd)

oligo.nhm.panel <- results.oligo.nhm%>%
  filter(Gene == "Any Amplification" | Gene == "Any HD")%>%
  mutate(p.adjust = p.adjust(as.numeric(pvalue), method = "fdr"))

print(oligo.nhm.panel)

# Run McNemar's test for CNVs in all hypermutant astros
cnv.data.astro.hm <- cnv.data %>%
  filter(idh_codel_subtype == "IDHmut-noncodel")%>%
  filter(HM_group == "HM")

results.astro.hm <- tibble(Gene = character(), pvalue = numeric())

for (gene in unique(cnv.data.astro.hm$gene_symbol)) {
  test_data <- cnv.data.astro.hm %>% filter(gene_symbol == gene) 
  gene <- unique(test_data$gene_symbol)
  mcnemar.matrix <-
    matrix(c(nrow(test_data %>% filter(driver_change == "S")), 
             nrow(test_data %>% filter(driver_change == "R")),     
             nrow(test_data %>% filter(driver_change == "P")), 
             nrow(test_data %>% filter(is.na(driver_change )))
    ),
    nrow = 2,
    dimnames = list("Primary" = c("AMP/HD", "No_amp/hd"),
                    "Recurrence" = c("AMP/HD", "No_amp/hd")))
  pvalue <- mcnemar.test(mcnemar.matrix)$p.value
  results.astro.hm <- results.astro.hm %>% add_row(Gene = as.character(gene), pvalue = pvalue)
}

astro.hm.single.amp <- results.astro.hm%>%
  filter(Gene == "CCND2" | Gene == "CDK4" | Gene == "CDK6" | Gene == "EGFR" | Gene == "MET" |Gene == "MYC" |Gene == "MYCN" |Gene == "PDGFRA" |Gene == "PTPRZ1")%>%
  mutate(p.adjust = p.adjust(as.numeric(pvalue), method = "fdr"))

print(astro.hm.single.amp)

astro.hm.single.hd <- results.astro.hm%>%
  filter(Gene == "ATRX" | Gene == "CDKN2A" | Gene == "CDKN2B" | Gene == "PTEN")%>%
  mutate(p.adjust = p.adjust(as.numeric(pvalue), method = "fdr"))

print(astro.hm.single.hd)

astro.hm.panel <- results.astro.hm%>%
  filter(Gene == "Any Amplification" | Gene == "Any HD")%>%
  mutate(p.adjust = p.adjust(as.numeric(pvalue), method = "fdr"))

print(astro.hm.panel)

# Run McNemar's test for CNVs in all non-hypermutant astros
cnv.data.astro.nhm <- cnv.data %>%
  filter(idh_codel_subtype == "IDHmut-noncodel")%>%
  filter(HM_group == "NHM")

results.astro.nhm <- tibble(Gene = character(), pvalue = numeric())

for (gene in unique(cnv.data.astro.nhm$gene_symbol)) {
  test_data <- cnv.data.astro.nhm %>% filter(gene_symbol == gene) 
  gene <- unique(test_data$gene_symbol)
  mcnemar.matrix <-
    matrix(c(nrow(test_data %>% filter(driver_change == "S")), 
             nrow(test_data %>% filter(driver_change == "R")),     
             nrow(test_data %>% filter(driver_change == "P")), 
             nrow(test_data %>% filter(is.na(driver_change )))
    ),
    nrow = 2,
    dimnames = list("Primary" = c("AMP/HD", "No_amp/hd"),
                    "Recurrence" = c("AMP/HD", "No_amp/hd")))
  pvalue <- mcnemar.test(mcnemar.matrix)$p.value
  results.astro.nhm <- results.astro.nhm %>% add_row(Gene = as.character(gene), pvalue = pvalue)
}

astro.nhm.single.amp <- results.astro.nhm%>%
  filter(Gene == "CCND2" | Gene == "CDK4" | Gene == "CDK6" | Gene == "EGFR" | Gene == "MET" |Gene == "MYC" |Gene == "MYCN" |Gene == "PDGFRA" |Gene == "PTPRZ1")%>%
  mutate(p.adjust = p.adjust(as.numeric(pvalue), method = "fdr"))

print(astro.nhm.single.amp)

astro.nhm.single.hd <- results.astro.nhm%>%
  filter(Gene == "ATRX" | Gene == "CDKN2A" | Gene == "CDKN2B" | Gene == "PTEN")%>%
  mutate(p.adjust = p.adjust(as.numeric(pvalue), method = "fdr"))

print(astro.nhm.single.hd)

astro.nhm.panel <- results.astro.nhm%>%
  filter(Gene == "Any Amplification" | Gene == "Any HD")%>%
  mutate(p.adjust = p.adjust(as.numeric(pvalue), method = "fdr"))

print(astro.nhm.panel)