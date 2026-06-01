# GLASS-I: Mutual Exclusivity Driver sSNVs
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(tidyverse)
library(RPostgres)
library(corrplot)
library(discover)

# load data
molecular.data <- read.delim("~/driver_changes.txt")%>%                       ## created with GLASS-I_Drivers_Table.R  
  mutate(timing.recurrent = substr(tumor_pair_barcode, 20,21),
         timing.primary = substr(tumor_pair_barcode, 14,15) )%>%
  select(-idh_codel_subtype) 

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")             ## data available on Synapse (Synapse ID: syn17038081)

surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
  mutate(timing.surgery = substr(sample_barcode, 14,15)) %>%
  filter(idh_codel_subtype != "IDHwt")

dbDisconnect(con)

snv.at.primary <- surgeries%>%
  left_join(molecular.data, by ="case_barcode")%>%
  filter(surgery_number == 1) %>% 
  filter(timing.primary == timing.surgery) %>%  
  group_by(case_barcode)%>%
  summarise(HM_group =         HM_group[1],
            tumor_type =       if_else(idh_codel_subtype[1] != "IDHmut-codel", "IDH-mutant Astrocytomas", "IDH-mutant Oligodendrogliomas"),
            NOTCH1 =           if_else(any(gene_symbol == "NOTCH1" & driver_status == "mutant" & driver_change == "P"), 1, 
                               if_else(any(gene_symbol == "NOTCH1" & driver_status == "mutant" & driver_change == "S"), 1, 0)),
            NOTCH2 =           if_else(any(gene_symbol == "NOTCH2" & driver_status == "mutant" & driver_change == "P"), 1,
                               if_else(any(gene_symbol == "NOTCH2" & driver_status == "mutant" & driver_change == "S"), 1, 0)),
            NOTCH3 =           if_else(any(gene_symbol == "NOTCH3" & driver_status == "mutant" & driver_change == "P"), 1, 
                               if_else(any(gene_symbol == "NOTCH3" & driver_status == "mutant" & driver_change == "S"), 1, 0)),
            PIK3CA =           if_else(any(gene_symbol == "PIK3CA" & driver_status == "mutant" & driver_change == "P"), 1, 
                               if_else(any(gene_symbol == "PIK3CA" & driver_status == "mutant" & driver_change == "S"), 1, 0)),
            PIK3R1 =           if_else(any(gene_symbol == "PIK3R1" & driver_status == "mutant" & driver_change == "P"), 1, 
                               if_else(any(gene_symbol == "PIK3R1" & driver_status == "mutant" & driver_change == "S"), 1, 0)),
            PI3K =             if_else(any(gene_symbol == "PI3K" & driver_status == "mutant" & driver_change == "P"), 1, 
                               if_else(any(gene_symbol == "PI3K" & driver_status == "mutant" & driver_change == "S"), 1, 0))
  )%>%
  ungroup()%>%
  mutate(  NOTCH =             if_else(NOTCH1 ==1 | NOTCH2 == 1 | NOTCH3 == 1, 1, 0))

snv.at.recurrence <- surgeries%>%
  left_join(molecular.data, by ="case_barcode")%>%
  filter(timing.recurrent == timing.surgery) %>%  
  group_by(case_barcode)%>%
  summarise(HM_group =         HM_group[1],
            tumor_type =       if_else(idh_codel_subtype[1] != "IDHmut-codel", "IDH-mutant Astrocytomas", "IDH-mutant Oligodendrogliomas"),
            NOTCH1 =           if_else(any(gene_symbol == "NOTCH1" & driver_status == "mutant" & driver_change == "R"), 1, 
                               if_else(any(gene_symbol == "NOTCH1" & driver_status == "mutant" & driver_change == "S"), 1, 0)),
            NOTCH2 =           if_else(any(gene_symbol == "NOTCH2" & driver_status == "mutant" & driver_change == "R"), 1, 
                               if_else(any(gene_symbol == "NOTCH2" & driver_status == "mutant" & driver_change == "S"), 1, 0)),
            NOTCH3 =           if_else(any(gene_symbol == "NOTCH3" & driver_status == "mutant" & driver_change == "R"), 1, 
                               if_else(any(gene_symbol == "NOTCH3" & driver_status == "mutant" & driver_change == "S"), 1, 0)),
            PIK3CA =           if_else(any(gene_symbol == "PIK3CA" & driver_status == "mutant" & driver_change == "R"), 1, 
                               if_else(any(gene_symbol == "PIK3CA" & driver_status == "mutant" & driver_change == "S"), 1, 0)),
            PIK3R1 =           if_else(any(gene_symbol == "PIK3R1" & driver_status == "mutant" & driver_change == "R"), 1, 
                               if_else(any(gene_symbol == "PIK3R1" & driver_status == "mutant" & driver_change == "S"), 1, 0)),
            PI3K =             if_else(any(gene_symbol == "PI3K" & driver_status == "mutant" & driver_change == "R"), 1, 
                               if_else(any(gene_symbol == "PI3K" & driver_status == "mutant" & driver_change == "S"), 1, 0))
  )%>%
  ungroup()%>%
  mutate(  NOTCH =             if_else(NOTCH1 ==1 | NOTCH2 == 1 | NOTCH3 == 1, 1, 0))

# Primary
snv.at.primary.oligo <- snv.at.primary %>%
  filter(tumor_type == "IDH-mutant Oligodendrogliomas")%>%
  select(-case_barcode, -tumor_type, -HM_group)%>%
  drop_na()

events <- discover.matrix(t(snv.at.primary.oligo))
result.mutex <- pairwise.discover.test(events[, ])
df <- as.data.frame(result.mutex, q.threshold=1)[, c("gene1", "gene2", "q.value")]
genes <- sort(unique(c(df$gene1, df$gene2)))
gene_matrix <- matrix(NA, nrow = length(genes), ncol = length(genes), dimnames = list(genes, genes))

for (i in 1:nrow(df)) {
  gene_matrix[df$gene1[i], df$gene2[i]] <- df$q.value[i]
  gene_matrix[df$gene2[i], df$gene1[i]] <- df$q.value[i]
}

p.mat.oligo.primary <- apply(gene_matrix, 2, as.numeric)
rownames(p.mat.oligo.primary) <- genes

corrplot(cor(snv.at.primary.oligo, method = "spearman"),
         order = "alphabet",
         title = "Correlations NOTCH and PI3K mutations \nat primary surgery",
         type="upper", 
         col=c("black", "white"),
         diag=F,
         tl.pos = "n",
         method = "shade",
         p.mat = p.mat.oligo.primary,
         sig.level = c(0.001, 0.01, 0.05), 
         mar = c(10,0,3,0),
         pch.cex = 4,
         insig = 'label_sig', 
         na.label =  "NA",
         pch.col = '#DAEEF8',
         addgrid.col = '#298C8C',
         na.label.col = "#298C8C",
         bg="#DAEEF8")

snv.at.primary.astro <- snv.at.primary %>%
  filter(tumor_type == "IDH-mutant Astrocytomas")%>%
  select(-case_barcode, -tumor_type, -HM_group)%>%
  drop_na()

events <- discover.matrix(t(snv.at.primary.astro))
result.mutex <- pairwise.discover.test(events[, ])
df <- as.data.frame(result.mutex, q.threshold=1)[, c("gene1", "gene2", "q.value")]
genes <- sort(unique(c(df$gene1, df$gene2)))
gene_matrix <- matrix(NA, nrow = length(genes), ncol = length(genes), dimnames = list(genes, genes))

for (i in 1:nrow(df)) {
  gene_matrix[df$gene1[i], df$gene2[i]] <- df$q.value[i]
  gene_matrix[df$gene2[i], df$gene1[i]] <- df$q.value[i]
}

p.mat.astro.primary <- apply(gene_matrix, 2, as.numeric)
rownames(p.mat.astro.primary) <- genes

corrplot(cor(snv.at.primary.astro, method = "spearman"), 
         order = 'alphabet',
         type="lower", 
         col=c("black", "white"),
         diag=F,
         cl.pos = "n",
         tl.pos = "ld",
         tl.col="black", 
         tl.srt=90,
         tl.cex =1.25,
         method = "shade",
         p.mat = p.mat.astro.primary,
         sig.level = c(0.001, 0.01, 0.05), 
         pch.cex = 4,
         insig = 'label_sig', 
         pch.col = '#F5E2EE',
         addgrid.col = '#800074',
         na.label =  "NA",
         na.label.col = '#800074',
         add=T,
         bg="#F5E2EE")

cor.primary.plot <- recordPlot()

# HM at recurrence
snv.at.recurrence.hm <- snv.at.recurrence %>%
  filter(HM_group == "HM")

snv.at.recurrence.hm.oligo <- snv.at.recurrence.hm %>%
  filter(tumor_type == "IDH-mutant Oligodendrogliomas")%>%
  select(-case_barcode, -tumor_type, -HM_group)%>%
  drop_na()

events <- discover.matrix(t(snv.at.recurrence.hm.oligo))
result.mutex <- pairwise.discover.test(events[, ])
df <- as.data.frame(result.mutex, q.threshold=1)[, c("gene1", "gene2", "q.value")]
genes <- sort(unique(c(df$gene1, df$gene2)))
gene_matrix <- matrix(NA, nrow = length(genes), ncol = length(genes), dimnames = list(genes, genes))

for (i in 1:nrow(df)) {
  gene_matrix[df$gene1[i], df$gene2[i]] <- df$q.value[i]
  gene_matrix[df$gene2[i], df$gene1[i]] <- df$q.value[i]
}

p.mat.oligo.recurrence.hm <- apply(gene_matrix, 2, as.numeric)
rownames(p.mat.oligo.recurrence.hm) <- genes

corrplot(cor(snv.at.recurrence.hm.oligo, method = "spearman"), 
         order = "alphabet",
         title = "Correlations NOTCH and PI3K mutations \nat recurrence in hypermutant gliomas",
         type="upper", 
         col=c("black", "white"),
         diag=F,
         tl.pos = "n",
         method = "shade",
         p.mat = p.mat.oligo.recurrence.hm,
         sig.level = c(0.001, 0.01, 0.05), 
         mar = c(10,0,3,0),
         pch.cex = 4,
         insig = 'label_sig', 
         na.label =  "NA",
         pch.col = '#DAEEF8',
         addgrid.col = '#298C8C',
         na.label.col = "#298C8C",
         bg="#DAEEF8")

snv.at.recurrence.hm.astro <- snv.at.recurrence.hm %>%
  filter(tumor_type == "IDH-mutant Astrocytomas")%>%
  select(-case_barcode, -tumor_type, -HM_group)%>%
  drop_na()

events <- discover.matrix(t(snv.at.recurrence.hm.astro))
result.mutex <- pairwise.discover.test(events[, ])
df <- as.data.frame(result.mutex, q.threshold=1)[, c("gene1", "gene2", "q.value")]
genes <- sort(unique(c(df$gene1, df$gene2)))
gene_matrix <- matrix(NA, nrow = length(genes), ncol = length(genes), dimnames = list(genes, genes))

for (i in 1:nrow(df)) {
  gene_matrix[df$gene1[i], df$gene2[i]] <- df$q.value[i]
  gene_matrix[df$gene2[i], df$gene1[i]] <- df$q.value[i]
}

p.mat.astro.recurrence.hm <- apply(gene_matrix, 2, as.numeric)
rownames(p.mat.astro.recurrence.hm) <- genes

corrplot(cor(snv.at.recurrence.hm.astro, method = "spearman"), 
         order = 'alphabet',
         type="lower", 
         col=c("black", "white"),
         diag=F,
         cl.pos = "n",
         tl.pos = "ld",
         tl.col="black", 
         tl.srt=90,
         tl.cex =1.25,
         method = "shade",
         p.mat = p.mat.astro.recurrence.hm,
         sig.level = c(0.001, 0.01, 0.05), 
         pch.cex = 4,
         insig = 'label_sig', 
         pch.col = '#F5E2EE',
         addgrid.col = '#800074',
         na.label =  "NA",
         na.label.col = '#800074',
         add=T,
         bg="#F5E2EE")

cor.recurrence.hm.plot <- recordPlot()

# NHM at recurrence
snv.at.recurrence.nhm <- snv.at.recurrence %>%
  filter(HM_group == "NHM")

snv.at.recurrence.nhm.oligo <- snv.at.recurrence.nhm %>%
  filter(tumor_type == "IDH-mutant Oligodendrogliomas")%>%
  select(-case_barcode, -tumor_type, -HM_group)%>%
  drop_na()

events <- discover.matrix(t(snv.at.recurrence.nhm.oligo))
result.mutex <- pairwise.discover.test(events[, ])
df <- as.data.frame(result.mutex, q.threshold=1)[, c("gene1", "gene2", "q.value")]
genes <- sort(unique(c(df$gene1, df$gene2)))
gene_matrix <- matrix(NA, nrow = length(genes), ncol = length(genes), dimnames = list(genes, genes))

for (i in 1:nrow(df)) {
  gene_matrix[df$gene1[i], df$gene2[i]] <- df$q.value[i]
  gene_matrix[df$gene2[i], df$gene1[i]] <- df$q.value[i]
}

p.mat.oligo.recurrence.nhm <- apply(gene_matrix, 2, as.numeric)
rownames(p.mat.oligo.recurrence.nhm) <- genes

corrplot(cor(snv.at.recurrence.nhm.oligo, method = "spearman"), 
         order = "alphabet",
         title = "Correlations NOTCH and PI3K mutations \nat recurrence in non-hypermutant gliomas",
         type="upper", 
         col=c("black", "white"),
         diag=F,
         tl.pos = "n",
         method = "shade",
         p.mat = p.mat.oligo.recurrence.nhm,
         sig.level = c(0.001, 0.01, 0.05), 
         mar = c(10,0,3,0),
         pch.cex = 4,
         insig = 'label_sig', 
         na.label =  "NA",
         pch.col = '#DAEEF8',
         addgrid.col = '#298C8C',
         na.label.col = "#298C8C",
         bg="#DAEEF8")

snv.at.recurrence.nhm.astro <- snv.at.recurrence.nhm %>%
  filter(tumor_type == "IDH-mutant Astrocytomas")%>%
  select(-case_barcode, -tumor_type, -HM_group)%>%
  drop_na()

events <- discover.matrix(t(snv.at.recurrence.nhm.astro))
result.mutex <- pairwise.discover.test(events[, ])
df <- as.data.frame(result.mutex, q.threshold=1)[, c("gene1", "gene2", "q.value")]
genes <- sort(unique(c(df$gene1, df$gene2)))
gene_matrix <- matrix(NA, nrow = length(genes), ncol = length(genes), dimnames = list(genes, genes))

for (i in 1:nrow(df)) {
  gene_matrix[df$gene1[i], df$gene2[i]] <- df$q.value[i]
  gene_matrix[df$gene2[i], df$gene1[i]] <- df$q.value[i]
}

p.mat.astro.recurrence.nhm <- apply(gene_matrix, 2, as.numeric)
rownames(p.mat.astro.recurrence.nhm) <- genes

corrplot(cor(snv.at.recurrence.nhm.astro, method = "spearman"), 
         order = 'alphabet',
         type="lower", 
         col=c("black", "white"),
         diag=F,
         cl.pos = "n",
         tl.pos = "ld",
         tl.col="black", 
         tl.srt=90,
         tl.cex =1.25,
         method = "shade",
         p.mat = p.mat.astro.recurrence.nhm,
         sig.level = c(0.001, 0.01, 0.05), 
         pch.cex = 4,
         insig = 'label_sig', 
         pch.col = '#F5E2EE',
         addgrid.col = '#800074',
         na.label =  "NA",
         na.label.col = '#800074',
         add=T,
         bg="#F5E2EE")

cor.recurrence.nhm.plot <- recordPlot()

cor.primary.plot
cor.recurrence.hm.plot
cor.recurrence.nhm.plot