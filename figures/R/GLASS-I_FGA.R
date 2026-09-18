# GLASS-I: FGA
# Author: C.M.S. Tesileanu
# Date: 2026-03-04

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggh4x)

molecular.data <- read.delim("~/driver_changes.txt")%>%                       ## created with GLASS-I_Drivers_Table.R
  group_by(case_barcode)%>%
  summarise(driver = if_else(any(HM_group == "HM"), "HM", "NHM"))%>%
  ungroup()

cnv_summary <- read.delim("~/fga_calc_annotated.txt")%>%           ## created with script_fga_calculation.R
  left_join(molecular.data, by = "case_barcode")

cnv_summary_all <- cnv_summary %>%
  mutate(HM_group = "All")

cnv_summary_hm_groups <- cnv_summary %>%
  filter(HM_group == "HM"| HM_group == "NHM")

cnv_summary_combined <- bind_rows(
  cnv_summary_all,
  cnv_summary_hm_groups
)

# Set factor order so All appears first, then HM, then NHM
cnv_summary_combined$HM_group <- factor(
  cnv_summary_combined$HM_group,
  levels = c("All", "HM", "NHM")
)

gg_fga_norm_combined <- ggplot(
  cnv_summary_combined,
  aes(x = sample_type, y = fga_norm, fill = sample_type)
) +
  geom_violin(trim = FALSE) +
  geom_line(aes(group = case_barcode, color = platform)) +
  geom_point(position = position_jitter(width = 0.05), size = 2) +
  stat_compare_means(
    method = "wilcox.test",
    paired = TRUE
  ) +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = round(..y.., 2)),
    position = position_nudge(x = 0.5),
    size = 5,
    color = "black"
  ) +
  facet_grid2(
    idh_codel_subtype ~  HM_group,
    scales = "free",
    strip = strip_nested(
      bleed = FALSE,
      background_x = element_rect(color = "black"),
      background_y = element_rect(color = "black")
    )
  ) +
  scale_fill_manual(values = c("#CA2F66", "#2FB3CA")) +
  scale_color_manual(values = c("#98c127", "#ffcd8e")) +
  theme_classic(base_size = 20) +
  labs(title = "", x = "", y = "Fraction of genome altered") +
  theme(axis.text.x = element_text(angle = 90)) +
  ylim(0, 1.15)

gg_fga_norm_combined

oligo <- cnv_summary_combined%>%
  filter(idh_codel_subtype == "Oligodendroglioma")
nlevels(as.factor(oligo$case_barcode))

non.hypermut.oligo <- oligo%>%
  filter(HM_group == "NHM")
nlevels(as.factor(non.hypermut.oligo$case_barcode))

hypermut.oligo <- oligo%>%
  filter(HM_group == "HM")
nlevels(as.factor(hypermut.oligo$case_barcode))

astro <- cnv_summary_combined%>%
  filter(idh_codel_subtype == "Astrocytoma")
nlevels(as.factor(astro$case_barcode))

non.hypermut.astro <- astro%>%
  filter(HM_group == "NHM")
nlevels(as.factor(non.hypermut.astro$case_barcode))

hypermut.astro <- astro%>%
  filter(HM_group == "HM")
nlevels(as.factor(hypermut.astro$case_barcode))

paired_df <- cnv_summary_combined %>%
  select(case_barcode, platform, idh_codel_subtype, HM_group, sample_type, fga_norm) %>%
  pivot_wider(
    names_from = sample_type,
    values_from = fga_norm
  ) %>%
  drop_na(Primary, Recurrent)

oligo.wgs <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WGS")
wilcox.test(oligo.wgs$Primary, oligo.wgs$Recurrent, paired = TRUE)

oligo.wgs.hm <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WGS" & HM_group == "HM")
wilcox.test(oligo.wgs.hm$Primary, oligo.wgs.hm$Recurrent, paired = TRUE)

oligo.wgs.nhm <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WGS" & HM_group == "NHM")
wilcox.test(oligo.wgs.nhm$Primary, oligo.wgs.nhm$Recurrent, paired = TRUE)

oligo.wxs <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WXS")
wilcox.test(oligo.wxs$Primary, oligo.wxs$Recurrent, paired = TRUE)

oligo.wxs.hm <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WXS" & HM_group == "HM")
wilcox.test(oligo.wxs.hm$Primary, oligo.wxs.hm$Recurrent, paired = TRUE)

oligo.wxs.nhm <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WXS" & HM_group == "NHM")
wilcox.test(oligo.wxs.nhm$Primary, oligo.wxs.nhm$Recurrent, paired = TRUE)

astro.wgs <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WGS")
wilcox.test(astro.wgs$Primary, astro.wgs$Recurrent, paired = TRUE)

astro.wgs.hm <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WGS" & HM_group == "HM")
wilcox.test(astro.wgs.hm$Primary, astro.wgs.hm$Recurrent, paired = TRUE)

astro.wgs.nhm <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WGS" & HM_group == "NHM")
wilcox.test(astro.wgs.nhm$Primary, astro.wgs.nhm$Recurrent, paired = TRUE)

astro.wxs <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WXS")
wilcox.test(astro.wxs$Primary, astro.wxs$Recurrent, paired = TRUE)

astro.wxs.hm <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WXS" & HM_group == "HM")
wilcox.test(astro.wxs.hm$Primary, astro.wxs.hm$Recurrent, paired = TRUE)

astro.wxs.nhm <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WXS" & HM_group == "NHM")
wilcox.test(astro.wxs.nhm$Primary, astro.wxs.nhm$Recurrent, paired = TRUE)

ggsave("fga_norm_combined.pdf", plot = gg_fga_norm_combined, width = 20, height = 15, units = "in", dpi = 300)