# GLASS-I: FGA Non-Hypermutants
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggh4x)

molecular.data <- read.delim("~/driver_changes.txt")%>%                       ## created with GLASS-I_Drivers_Table.R
  mutate(driver = case_when(idh_codel_subtype == "IDHmut-noncodel" & gene_symbol == "Any Amplification" & driver_status == "HLAMP" & driver_change != "P" ~ "Driver+",
                            idh_codel_subtype == "IDHmut-noncodel" & gene_symbol == "CDKN2A/CDKN2B" & driver_status == "HLDEL" & driver_change != "P" ~ "Driver+",
                            idh_codel_subtype == "IDHmut-codel" & gene_symbol == "PI3K" & driver_status == "mutant" & driver_change != "P" ~ "Driver+",
                            idh_codel_subtype == "IDHmut-codel" & gene_symbol == "NOTCH1" & driver_status == "mutant" & driver_change != "P" ~ "Driver+",
                            TRUE ~ "Driver-")
         )%>%
  group_by(case_barcode)%>%
  summarise(driver = if_else(any(driver == "Driver+"), "Driver+", "Driver-"))%>%
  ungroup()

cnv_summary_filter_primary <- read.delim("~/fga_calc_annotated.txt")%>%       ## created with script_fga_calculation.R
  filter(HM_group == "NHM")%>%
  left_join(molecular.data, by = "case_barcode")

gg_fga_norm_combine_all_violin <- ggplot(
  cnv_summary_filter_primary,
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
    HM_group ~ idh_codel_subtype,
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

gg_fga_norm_combine_all_violin

gg_fga_norm_combine_driver_violin <- ggplot(
  cnv_summary_filter_primary,
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
    HM_group ~ idh_codel_subtype + driver,
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

gg_fga_norm_combine_driver_violin

non.hypermut.astro <- cnv_summary_filter_primary%>%
  filter(idh_codel_subtype == "Astrocytoma")
nlevels(as.factor(non.hypermut.astro$case_barcode))

non.hypermut.astro.cnv.pos <- non.hypermut.astro%>%
  filter(driver == "Driver+")
nlevels(as.factor(non.hypermut.astro.cnv.pos$case_barcode))

non.hypermut.astro.cnv.min <- non.hypermut.astro%>%
  filter(driver == "Driver-")
nlevels(as.factor(non.hypermut.astro.cnv.min$case_barcode))

non.hypermut.oligo <- cnv_summary_filter_primary%>%
  filter(idh_codel_subtype == "Oligodendroglioma")
nlevels(as.factor(non.hypermut.oligo$case_barcode))

non.hypermut.oligo.snv.pos <- non.hypermut.oligo%>%
  filter(driver == "Driver+")
nlevels(as.factor(non.hypermut.oligo.snv.pos$case_barcode))

non.hypermut.oligo.snv.min <- non.hypermut.oligo%>%
  filter(driver == "Driver-")
nlevels(as.factor(non.hypermut.oligo.snv.min$case_barcode))

paired_df <- cnv_summary_filter_primary %>%
  select(case_barcode, platform, idh_codel_subtype, driver, sample_type, fga_norm) %>%
  pivot_wider(
    names_from = sample_type,
    values_from = fga_norm
  ) %>%
  drop_na(Primary, Recurrent)

astro.wgs <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WGS")
wilcox.test(astro.wgs$Primary, astro.wgs$Recurrent, paired = TRUE)

astro.wgs.driv.plus <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WGS" & driver == "Driver+")
wilcox.test(astro.wgs.driv.plus$Primary, astro.wgs.driv.plus$Recurrent, paired = TRUE)

astro.wgs.driv.min <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WGS" & driver == "Driver-")
wilcox.test(astro.wgs.driv.min$Primary, astro.wgs.driv.min$Recurrent, paired = TRUE)

astro.wxs <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WXS")
wilcox.test(astro.wxs$Primary, astro.wxs$Recurrent, paired = TRUE)

astro.wxs.driv.plus <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WXS" & driver == "Driver+")
wilcox.test(astro.wxs.driv.plus$Primary, astro.wxs.driv.plus$Recurrent, paired = TRUE)

astro.wxs.driv.min <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WXS" & driver == "Driver-")
wilcox.test(astro.wxs.driv.min$Primary, astro.wxs.driv.min$Recurrent, paired = TRUE)

oligo.wgs <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WGS")
wilcox.test(oligo.wgs$Primary, oligo.wgs$Recurrent, paired = TRUE)

oligo.wgs.driv.plus <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WGS" & driver == "Driver+")
wilcox.test(oligo.wgs.driv.plus$Primary, oligo.wgs.driv.plus$Recurrent, paired = TRUE)

oligo.wgs.driv.min <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WGS" & driver == "Driver-")
wilcox.test(oligo.wgs.driv.min$Primary, oligo.wgs.driv.min$Recurrent, paired = TRUE)

oligo.wxs <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WXS")
wilcox.test(oligo.wxs$Primary, oligo.wxs$Recurrent, paired = TRUE)

oligo.wxs.driv.plus <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WXS" & driver == "Driver+")
wilcox.test(oligo.wxs.driv.plus$Primary, oligo.wxs.driv.plus$Recurrent, paired = TRUE)

oligo.wxs.driv.min <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WXS" & driver == "Driver-")
wilcox.test(oligo.wxs.driv.min$Primary, oligo.wxs.driv.min$Recurrent, paired = TRUE)

ggsave("~/fga_norm_combine_all_violin.pdf", plot = gg_fga_norm_combine_all_violin, width = 14, height = 10, units = "in", dpi = 300)
ggsave("~/fga_norm_combine_driver_violin.pdf", plot = gg_fga_norm_combine_driver_violin, width = 20, height = 10, units = "in", dpi = 300)