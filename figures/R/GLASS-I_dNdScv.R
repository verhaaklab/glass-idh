# GLASS-I: dNdScv Analysis
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

library(dndscv)
library(tidyverse)
library(RPostgres)
library(gridExtra)
library(ggplot2)
library(patchwork)

# load data
data("cancergenes_cgc81", package="dndscv")

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")                   ## data available on Synapse (Synapse ID: syn17038081)

res_fraction_nhm <- dbGetQuery(con, read_file("~/GLASS-I_dNdScv_NHM_input.sql"))    ## see SQL script of same name
res_fraction_hm <- dbGetQuery(con, read_file("~/GLASS-I_dNdScv_HM_input.sql"))      ## see SQL script of same name

dbDisconnect(con)

###################
## Run dNdSCV seperately for private/shared variants and IDH status in non-hypermutants
###################

result_list_nhm <- lapply(na.omit(unique(res_fraction_nhm$subtype)), function(st) {
  result_list_nhm <- lapply(na.omit(unique(res_fraction_nhm$fraction)), function(fr) {
    message("Computing dNdS for ", fr, " and ", st)
    qres_subset = res_fraction_nhm %>% filter(fraction == fr, subtype == st) %>% select(case_barcode,chrom,pos,ref,mut) %>% distinct()
    dnds_subset = dndscv(qres_subset, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = NULL)
    message(".. Ran dNdSCV for ", n_distinct(qres_subset$case_barcode) - length(dnds_subset$exclsamples))
    globaldnds <- cbind(fraction = fr, subtype = st, dnds_subset$globaldnds)
    sel_cv <- cbind(fraction = fr, subtype = st, dnds_subset$sel_cv)
    gene_ci <- cbind(fraction = fr, subtype = st, geneci(dnds_subset, gene_list = known_cancergenes))
    list(globaldnds, sel_cv, gene_ci)
  })
  res1 <- data.table::rbindlist(lapply(result_list_nhm,'[[',1))
  res2 <- data.table::rbindlist(lapply(result_list_nhm,'[[',2))
  res3 <- data.table::rbindlist(lapply(result_list_nhm,'[[',3))
  list(res1,res2,res3)
})

dnds_fraction_sel_cv_nhm <- data.table::rbindlist(lapply(result_list_nhm,'[[',2))

###################
## Plot dNdSCV genes by subtype - fraction in non-hypermutants
###################

dat_nhm <- dnds_fraction_sel_cv_nhm %>%
  mutate(fraction = factor(fraction, levels = c("P", "S", "R"))) %>%
  complete(gene_name,fraction,subtype, fill = list(qglobal_cv=1)) %>%
  group_by(fraction,subtype) %>%
  arrange(qglobal_cv, pglobal_cv) %>%
  ungroup() %>%
  arrange(fraction)

myplots_nhm <- lapply(split(dat_nhm, paste(dat_nhm$subtype, dat_nhm$fraction))[c(1,3,2,4,6,5,7,9,8)], function(df){
  df$gene_name = factor(df$gene_name, levels = rev(unique(df$gene_name)))
  p <- ggplot(df[1:5,], aes(x=gene_name, y=-log10(qglobal_cv), fill=fraction), width=0.75) +
    geom_bar(stat="identity", position = position_dodge()) + 
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) +
    labs(y="-log10(p), FDR-adjusted", x="") +
    guides(fill = F) +
    facet_wrap(~ subtype, scales = "free_y") + 
    coord_flip(ylim = c(0,10)) +
    scale_fill_manual(values = c("P" = "#CA2F66","S" = "#CA932F", "R" = "#2fb3ca"), drop=F) +
    theme_linedraw(base_size = 16) +
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.ticks = element_line(size = 0))
  
  return(p)
})

pdf('~/GLASS-I_dNdSCV_NHM.pdf', height = 8, width = 14)
myplots_nhm$`IDHmut-codel P` + myplots_nhm$`IDHmut-codel S` + myplots_nhm$`IDHmut-codel R`+
myplots_nhm$`IDHmut-noncodel P` + myplots_nhm$`IDHmut-noncodel S` + myplots_nhm$`IDHmut-noncodel R`
dev.off()

###################
## Run dNdSCV seperately for private/shared variants and IDH status in hypermutants
###################

result_list_hm <- lapply(na.omit(unique(res_fraction_hm$subtype)), function(st) {
  result_list_hm <- lapply(na.omit(unique(res_fraction_hm$fraction)), function(fr) {
    message("Computing dNdS for ", fr, " and ", st)
    qres_subset = res_fraction_hm %>% filter(fraction == fr, subtype == st) %>% select(case_barcode,chrom,pos,ref,mut) %>% distinct()
    dnds_subset = dndscv(qres_subset, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = NULL)
    message(".. Ran dNdSCV for ", n_distinct(qres_subset$case_barcode) - length(dnds_subset$exclsamples))
    globaldnds <- cbind(fraction = fr, subtype = st, dnds_subset$globaldnds)
    sel_cv <- cbind(fraction = fr, subtype = st, dnds_subset$sel_cv)
    gene_ci <- cbind(fraction = fr, subtype = st, geneci(dnds_subset, gene_list = known_cancergenes))
    list(globaldnds, sel_cv, gene_ci)
  })
  res1 <- data.table::rbindlist(lapply(result_list_hm,'[[',1))
  res2 <- data.table::rbindlist(lapply(result_list_hm,'[[',2))
  res3 <- data.table::rbindlist(lapply(result_list_hm,'[[',3))
  list(res1,res2,res3)
})

dnds_fraction_sel_cv_hm <- data.table::rbindlist(lapply(result_list_hm,'[[',2))

###################
## Plot dNdSCV genes by subtype - fraction in hypermutants
###################

dat_hm <- dnds_fraction_sel_cv_hm %>%
  mutate(fraction = factor(fraction, levels = c("P", "S", "R"))) %>%
  complete(gene_name,fraction,subtype, fill = list(qglobal_cv=1)) %>%
  group_by(fraction,subtype) %>%
  arrange(qglobal_cv, pglobal_cv) %>%
  ungroup() %>%
  arrange(fraction)

myplots_hm <- lapply(split(dat_hm, paste(dat_hm$subtype, dat_hm$fraction))[c(1,3,2,4,6,5,7,9,8)], function(df){
  df$gene_name = factor(df$gene_name, levels = rev(unique(df$gene_name)))
  p <- ggplot(df[1:5,], aes(x=gene_name, y=-log10(qglobal_cv), fill=fraction), width=0.75) +
    geom_bar(stat="identity", position = position_dodge()) + 
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) +
    labs(y="-log10(p), FDR-adjusted", x="") +
    guides(fill = F) +
    facet_wrap(~ subtype, scales = "free_y") + 
    coord_flip(ylim = c(0,10)) +
    scale_fill_manual(values = c("P" = "#CA2F66","S" = "#CA932F", "R" = "#2fb3ca"), drop=F) +
    theme_linedraw(base_size = 16) +
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.ticks = element_line(size = 0))
  
  return(p)
})

pdf('~/GLASS-I_dNdSCV_HM.pdf', height = 8, width = 14)
myplots_hm$`IDHmut-codel P` + myplots_hm$`IDHmut-codel S` + myplots_hm$`IDHmut-codel R`+
myplots_hm$`IDHmut-noncodel P` + myplots_hm$`IDHmut-noncodel S` + myplots_hm$`IDHmut-noncodel R`
dev.off()

pdf('~/GLASS-I_dNdSCV_All.pdf', height = 9, width = 14)
myplots_nhm$`IDHmut-codel P` + myplots_nhm$`IDHmut-codel S` + myplots_nhm$`IDHmut-codel R`+
myplots_hm$`IDHmut-codel P` + myplots_hm$`IDHmut-codel S` + myplots_hm$`IDHmut-codel R`+
myplots_nhm$`IDHmut-noncodel P` + myplots_nhm$`IDHmut-noncodel S` + myplots_nhm$`IDHmut-noncodel R`+
myplots_hm$`IDHmut-noncodel P` + myplots_hm$`IDHmut-noncodel S` + myplots_hm$`IDHmut-noncodel R`+ plot_layout(nrow = 4)
dev.off()