# GLASS-I: FGA Grade-Specified Non-Hypermutants
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
library(RPostgres)

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

cdkn2ab.data <- read.delim("~/driver_changes.txt")%>%                         ## created with GLASS-I_Drivers_Table.R
  filter(idh_codel_subtype == "IDHmut-noncodel" & gene_symbol == "CDKN2A/CDKN2B")%>%
  mutate(CDKN2AB = if_else(is.na(driver_change), "No HD at primary",
                   if_else(driver_change=="P"|driver_change=="S", "HD at primary", "No HD at primary")))%>%
  dplyr::select(case_barcode, CDKN2AB)

real.primary.surgeries <- read.delim("~/fga_calc_annotated.txt")%>%           ## created with script_fga_calculation.R
  filter(sample_type == "Primary")%>%
  mutate(sample_barcode = substr(aliquot_barcode, 1,15))

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")             ## data available on Synapse (Synapse ID: syn17038081)

surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
  mutate(timing.surgery = substr(sample_barcode, 14,15)) %>%
  filter(idh_codel_subtype != "IDHwt")%>%
  semi_join(real.primary.surgeries, by = "sample_barcode")%>%
  dplyr::select(case_barcode, grade)

dbDisconnect(con)

cnv_summary_filter_primary <- read.delim("~/fga_calc_annotated.txt")%>%       ## created with script_fga_calculation.R
  filter(HM_group == "NHM")%>%
  left_join(molecular.data, by = "case_barcode")%>%
  left_join(cdkn2ab.data, by = "case_barcode")%>%
  left_join(surgeries, by = "case_barcode")%>%
  mutate(grade_primary =  case_when(CDKN2AB == "HD at primary" ~ "4",
                                    CDKN2AB == "No HD at primary" & grade == "II" ~ "2",
                                    CDKN2AB == "No HD at primary" & grade == "III" ~ "3",
                                    CDKN2AB == "No HD at primary" & grade == "IV" ~ "4",
                                    is.na(CDKN2AB) & grade == "II" ~ "2",
                                    is.na(CDKN2AB) & grade == "III" ~ "3"))

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
    HM_group + grade_primary~ idh_codel_subtype,
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
    HM_group + grade_primary~ idh_codel_subtype + driver,
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

non.hypermut.astro.grade2 <- cnv_summary_filter_primary%>%
  filter(idh_codel_subtype == "Astrocytoma" & grade_primary == "2")
nlevels(as.factor(non.hypermut.astro.grade2$case_barcode))

non.hypermut.astro.grade2.cnv.pos <- non.hypermut.astro.grade2%>%
  filter(driver == "Driver+")
nlevels(as.factor(non.hypermut.astro.grade2.cnv.pos$case_barcode))

non.hypermut.astro.grade2.cnv.min <- non.hypermut.astro.grade2%>%
  filter(driver == "Driver-")
nlevels(as.factor(non.hypermut.astro.grade2.cnv.min$case_barcode))

non.hypermut.astro.grade3 <- cnv_summary_filter_primary%>%
  filter(idh_codel_subtype == "Astrocytoma" & grade_primary == "3")
nlevels(as.factor(non.hypermut.astro.grade3$case_barcode))

non.hypermut.astro.grade3.cnv.pos <- non.hypermut.astro.grade3%>%
  filter(driver == "Driver+")
nlevels(as.factor(non.hypermut.astro.grade3.cnv.pos$case_barcode))

non.hypermut.astro.grade3.cnv.min <- non.hypermut.astro.grade3%>%
  filter(driver == "Driver-")
nlevels(as.factor(non.hypermut.astro.grade3.cnv.min$case_barcode))

non.hypermut.astro.grade4 <- cnv_summary_filter_primary%>%
  filter(idh_codel_subtype == "Astrocytoma" & grade_primary == "4")
nlevels(as.factor(non.hypermut.astro.grade4$case_barcode))

non.hypermut.astro.grade4.cnv.pos <- non.hypermut.astro.grade4%>%
  filter(driver == "Driver+")
nlevels(as.factor(non.hypermut.astro.grade4.cnv.pos$case_barcode))

non.hypermut.astro.grade4.cnv.min <- non.hypermut.astro.grade4%>%
  filter(driver == "Driver-")
nlevels(as.factor(non.hypermut.astro.grade4.cnv.min$case_barcode))

non.hypermut.oligo.grade2 <- cnv_summary_filter_primary%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & grade_primary == "2")
nlevels(as.factor(non.hypermut.oligo.grade2$case_barcode))

non.hypermut.oligo.grade2.snv.pos <- non.hypermut.oligo.grade2%>%
  filter(driver == "Driver+")
nlevels(as.factor(non.hypermut.oligo.grade2.snv.pos$case_barcode))

non.hypermut.oligo.grade2.snv.min <- non.hypermut.oligo.grade2%>%
  filter(driver == "Driver-")
nlevels(as.factor(non.hypermut.oligo.grade2.snv.min$case_barcode))

non.hypermut.oligo.grade3 <- cnv_summary_filter_primary%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & grade_primary == "3")
nlevels(as.factor(non.hypermut.oligo.grade3$case_barcode))

non.hypermut.oligo.grade3.snv.pos <- non.hypermut.oligo.grade3%>%
  filter(driver == "Driver+")
nlevels(as.factor(non.hypermut.oligo.grade3.snv.pos$case_barcode))

non.hypermut.oligo.grade3.snv.min <- non.hypermut.oligo.grade3%>%
  filter(driver == "Driver-")
nlevels(as.factor(non.hypermut.oligo.grade3.snv.min$case_barcode))

paired_df <- cnv_summary_filter_primary %>%
  select(case_barcode, platform, idh_codel_subtype, driver, sample_type, fga_norm, grade_primary) %>%
  pivot_wider(
    names_from = sample_type,
    values_from = fga_norm
  ) %>%
  drop_na(Primary, Recurrent)

astro.wgs.grade2 <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WGS" & grade_primary == "2")
wilcox.test(astro.wgs.grade2$Primary, astro.wgs.grade2$Recurrent, paired = TRUE)

astro.wgs.grade2.driv.plus <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WGS" & driver == "Driver+" & grade_primary == "2")
wilcox.test(astro.wgs.grade2.driv.plus$Primary, astro.wgs.grade2.driv.plus$Recurrent, paired = TRUE)

astro.wgs.grade2.driv.min <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WGS" & driver == "Driver-" & grade_primary == "2")
wilcox.test(astro.wgs.grade2.driv.min$Primary, astro.wgs.grade2.driv.min$Recurrent, paired = TRUE)

astro.wxs.grade2 <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WXS" & grade_primary == "2")
wilcox.test(astro.wxs.grade2$Primary, astro.wxs.grade2$Recurrent, paired = TRUE)

astro.wxs.grade2.driv.plus <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WXS" & driver == "Driver+" & grade_primary == "2")
wilcox.test(astro.wxs.grade2.driv.plus$Primary, astro.wxs.grade2.driv.plus$Recurrent, paired = TRUE)

astro.wxs.grade2.driv.min <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WXS" & driver == "Driver-" & grade_primary == "2")
wilcox.test(astro.wxs.grade2.driv.min$Primary, astro.wxs.grade2.driv.min$Recurrent, paired = TRUE)

astro.wgs.grade3 <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WGS" & grade_primary == "3")
wilcox.test(astro.wgs.grade3$Primary, astro.wgs.grade3$Recurrent, paired = TRUE)

astro.wgs.grade3.driv.plus <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WGS" & driver == "Driver+" & grade_primary == "3")
wilcox.test(astro.wgs.grade3.driv.plus$Primary, astro.wgs.grade3.driv.plus$Recurrent, paired = TRUE)

astro.wgs.grade3.driv.min <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WGS" & driver == "Driver-" & grade_primary == "3")
wilcox.test(astro.wgs.grade3.driv.min$Primary, astro.wgs.grade3.driv.min$Recurrent, paired = TRUE)

astro.wxs.grade3 <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WXS" & grade_primary == "3")
wilcox.test(astro.wxs.grade3$Primary, astro.wxs.grade3$Recurrent, paired = TRUE)

astro.wxs.grade3.driv.plus <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WXS" & driver == "Driver+" & grade_primary == "3")
wilcox.test(astro.wxs.grade3.driv.plus$Primary, astro.wxs.grade3.driv.plus$Recurrent, paired = TRUE)

astro.wxs.grade3.driv.min <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WXS" & driver == "Driver-" & grade_primary == "3")
wilcox.test(astro.wxs.grade3.driv.min$Primary, astro.wxs.grade3.driv.min$Recurrent, paired = TRUE)

astro.wgs.grade4 <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WGS" & grade_primary == "4")
wilcox.test(astro.wgs.grade4$Primary, astro.wgs.grade4$Recurrent, paired = TRUE)

astro.wgs.grade4.driv.plus <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WGS" & driver == "Driver+" & grade_primary == "4")
wilcox.test(astro.wgs.grade4.driv.plus$Primary, astro.wgs.grade4.driv.plus$Recurrent, paired = TRUE)

astro.wgs.grade4.driv.min <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WGS" & driver == "Driver-" & grade_primary == "4")
wilcox.test(astro.wgs.grade4.driv.min$Primary, astro.wgs.grade4.driv.min$Recurrent, paired = TRUE)

astro.wxs.grade4 <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WXS" & grade_primary == "4")
wilcox.test(astro.wxs.grade4$Primary, astro.wxs.grade4$Recurrent, paired = TRUE)

astro.wxs.grade4.driv.plus <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WXS" & driver == "Driver+" & grade_primary == "4")
wilcox.test(astro.wxs.grade4.driv.plus$Primary, astro.wxs.grade4.driv.plus$Recurrent, paired = TRUE)

astro.wxs.grade4.driv.min <- paired_df%>%
  filter(idh_codel_subtype == "Astrocytoma" & platform == "WXS" & driver == "Driver-" & grade_primary == "4")
wilcox.test(astro.wxs.grade4.driv.min$Primary, astro.wxs.grade4.driv.min$Recurrent, paired = TRUE)

oligo.wgs.grade2 <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WGS" & grade_primary == "2")
wilcox.test(oligo.wgs.grade2$Primary, oligo.wgs.grade2$Recurrent, paired = TRUE)

oligo.wgs.grade2.driv.plus <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WGS" & driver == "Driver+" & grade_primary == "2")
wilcox.test(oligo.wgs.grade2.driv.plus$Primary, oligo.wgs.grade2.driv.plus$Recurrent, paired = TRUE)

oligo.wgs.grade2.driv.min <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WGS" & driver == "Driver-" & grade_primary == "2")
wilcox.test(oligo.wgs.grade2.driv.min$Primary, oligo.wgs.grade2.driv.min$Recurrent, paired = TRUE)

oligo.wxs.grade2 <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WXS" & grade_primary == "2")
wilcox.test(oligo.wxs.grade2$Primary, oligo.wxs.grade2$Recurrent, paired = TRUE)

oligo.wxs.grade2.driv.plus <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WXS" & driver == "Driver+" & grade_primary == "2")
wilcox.test(oligo.wxs.grade2.driv.plus$Primary, oligo.wxs.grade2.driv.plus$Recurrent, paired = TRUE)

oligo.wxs.grade2.driv.min <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WXS" & driver == "Driver-" & grade_primary == "2")
wilcox.test(oligo.wxs.grade2.driv.min$Primary, oligo.wxs.grade2.driv.min$Recurrent, paired = TRUE)

oligo.wgs.grade3 <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WGS" & grade_primary == "3")
wilcox.test(oligo.wgs.grade3$Primary, oligo.wgs.grade3$Recurrent, paired = TRUE)

oligo.wgs.grade3.driv.plus <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WGS" & driver == "Driver+" & grade_primary == "3")
wilcox.test(oligo.wgs.grade3.driv.plus$Primary, oligo.wgs.grade3.driv.plus$Recurrent, paired = TRUE)

oligo.wgs.grade3.driv.min <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WGS" & driver == "Driver-" & grade_primary == "3")
wilcox.test(oligo.wgs.grade3.driv.min$Primary, oligo.wgs.grade3.driv.min$Recurrent, paired = TRUE)

oligo.wxs.grade3 <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WXS" & grade_primary == "3")
wilcox.test(oligo.wxs.grade3$Primary, oligo.wxs.grade3$Recurrent, paired = TRUE)

oligo.wxs.grade3.driv.plus <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WXS" & driver == "Driver+" & grade_primary == "3")
wilcox.test(oligo.wxs.grade3.driv.plus$Primary, oligo.wxs.grade3.driv.plus$Recurrent, paired = TRUE)

oligo.wxs.grade3.driv.min <- paired_df%>%
  filter(idh_codel_subtype == "Oligodendroglioma" & platform == "WXS" & driver == "Driver-" & grade_primary == "3")
wilcox.test(oligo.wxs.grade3.driv.min$Primary, oligo.wxs.grade3.driv.min$Recurrent, paired = TRUE)

ggsave("~/fga_norm_combine_all_violin_grade.pdf", plot = gg_fga_norm_combine_all_violin, width = 14, height = 25, units = "in", dpi = 300)
ggsave("~/fga_norm_combine_driver_violin_grade.pdf", plot = gg_fga_norm_combine_driver_violin, width = 20, height = 25, units = "in", dpi = 300)