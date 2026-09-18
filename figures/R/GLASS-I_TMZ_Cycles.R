# GLASS-I: TMZ Cycles
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(RPostgres)
library(tidyverse)
library(patchwork)

# load data
hypermut.data <- read.delim("~/driver_changes.txt")%>%                        ## created with GLASS-I_Drivers_Table.R
  select(case_barcode, HM_group)%>%
  distinct(.)

prior.treatment <- read_csv("~/prior_treatment_gold_set.csv")%>%              ## created with GLASS-I_Prior_Treatment_Table.R
  select(-case_barcode, -surgery_number)

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")             ## data available on Synapse (Synapse ID: syn17038081)

surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
  mutate(timing.surgery = substr(sample_barcode, 14,15)) %>%
  filter(idh_codel_subtype != "IDHwt")%>%
  left_join(prior.treatment, by ="sample_barcode")

gold.set <- dbGetQuery(con, "SELECT * FROM analysis.gold_set")%>%
  mutate(timing.recurrent = substr(tumor_pair_barcode, 20,21),
         timing.primary = substr(tumor_pair_barcode, 14,15) )%>%
  semi_join(surgeries, by = "case_barcode")

initial.tmb.data <-  dbGetQuery(con, "SELECT * FROM analysis.tumor_clinical_comparison_2024") %>%
  filter(str_detect(tumor_pair_barcode, "-WGS"))%>%
  filter((sbs11_r_hypermutation == F & mutation_burden_a <10) |sbs11_r_hypermutation == T)%>%
  semi_join(gold.set, by = "case_barcode")

dbDisconnect(con)

signature.data <-  read_csv("~/sbs11_sbs119_id8_Rprop.csv")%>%                  ## created with mutational_signature_pipeline.sh
  semi_join(initial.tmb.data, by = "case_barcode")

sbs11.data <- signature.data %>%
  mutate(sig_prop_sbs11 = SBS11_R)%>%
  select(case_barcode, sig_prop_sbs11)

sbs119.data <- signature.data %>%
  mutate(sig_prop_sbs119 = SBS119_R)%>%
  select(case_barcode, sig_prop_sbs119)

tmz.dose.data.mut.sig <-   gold.set%>%
  left_join(surgeries, by ="case_barcode")%>%
  left_join(sbs119.data, by ="case_barcode")%>%
  left_join(sbs11.data, by ="case_barcode")%>%
  left_join(hypermut.data,by ="case_barcode")%>%
  filter(str_detect(tumor_pair_barcode, "-WGS"))%>%
  filter(timing.surgery== timing.recurrent )%>%
  mutate(tumor_type = case_when(idh_codel_subtype == "IDHmut-noncodel" ~ "Astrocytomas",
                                idh_codel_subtype == "IDHmut-codel" ~ "Oligodendrogliomas"),
         sig_prop_sbs11 = if_else(!is.na(sig_prop_sbs11) & HM_group == "HM", sig_prop_sbs11*100, NA),
         sig_prop_sbs119 = if_else(!is.na(sig_prop_sbs119) &HM_group == "NHM", sig_prop_sbs119*100, NA),
  )%>%
  filter(!is.na(TMZ_cycles))

tmz.dose.data.hypermut <-   gold.set%>%
  left_join(surgeries, by ="case_barcode")%>%
  left_join(hypermut.data,by ="case_barcode")%>%
  filter(timing.surgery== timing.recurrent)%>%
  mutate(tumor_type = case_when(idh_codel_subtype == "IDHmut-noncodel" ~ "Astrocytomas",
                                idh_codel_subtype == "IDHmut-codel" ~ "Oligodendrogliomas"),
         HM_group =  case_when(HM_group=="HM" ~ 1, 
                               HM_group=="NHM" ~0),
  )

astro.dose.hypermut <-   tmz.dose.data.hypermut%>%
  filter(tumor_type == "Astrocytomas")

oligo.dose.hypermut <-   tmz.dose.data.hypermut%>%
  filter(tumor_type == "Oligodendrogliomas")

# Create correlation plot with Spearman test
p1 <- ggplot(tmz.dose.data.mut.sig,aes(x = TMZ_cycles, y = sig_prop_sbs119, color = tumor_type)) +
  geom_smooth(method = "lm", color = "gray50", se = TRUE, fill = "gray85", linetype = "dashed") +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Astrocytomas" = "#800074","Oligodendrogliomas"  = "#298c8c")
  ) +
  stat_cor(
    data = tmz.dose.data.mut.sig,
    inherit.aes = FALSE,
    aes(
      x = TMZ_cycles,
      y = sig_prop_sbs119,
      label = paste(..r.label.., ..p.label.., sep = "~`,`~")
    ),
    method = "spearman",
    label.x.npc = "left",
    label.y.npc = "top"
  ) +
  stat_cor(
    method = "spearman",
    show.legend = FALSE,
    label.x.npc = "middle",
    label.y.npc = "top",
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))
  ) +
  labs(
    x = "Temozolomide cycles",
    y = "SBS119 recurrence proportion (%)",
    color = "Tumor type"
  ) +
  ggtitle("Temozolomide-treated") + 
  ylim(0,100)+
  theme_linedraw(base_size = 14)

# Create correlation plot with Spearman test
p2 <- ggplot(tmz.dose.data.mut.sig,aes(x = TMZ_cycles, y = sig_prop_sbs11, color = tumor_type)) +
  geom_smooth(method = "lm", color = "gray50", se = TRUE, fill = "gray85", linetype = "dashed") +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Astrocytomas" = "#800074","Oligodendrogliomas"  = "#298c8c")
  ) +
  stat_cor(
    data = tmz.dose.data.mut.sig,
    inherit.aes = FALSE,
    aes(
      x = TMZ_cycles,
      y = sig_prop_sbs11,
      label = paste(..r.label.., ..p.label.., sep = "~`,`~")
    ),
    method = "spearman",
    label.x.npc = "left",
    label.y.npc = "bottom"
  ) +
  stat_cor(
    method = "spearman",
    show.legend = FALSE,
    label.x.npc = "middle",
    label.y.npc = "bottom",
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))
  ) +
  labs(
    x = "Temozolomide cycles",
    y = "SBS11 recurrence proportion (%)",
    color = "Tumor type"
  ) +
  ggtitle("Temozolomide-treated") + 
  ylim(0,100)+
  xlim(0,25)+
  theme_linedraw(base_size = 14)

p.astro <- summary(glm(HM_group ~ TMZ_cycles,
                       data = astro.dose.hypermut,
                       family = binomial))$coefficients["TMZ_cycles", "Pr(>|z|)"]

p_label.astro <- ifelse(p.astro < 0.001, "p < 0.001", paste0("p = ", signif(p.astro, 3)))

p3 <- ggplot(astro.dose.hypermut, aes(x = TMZ_cycles, y =HM_group)) +
  geom_point(width = 0.3, height = 0.05, alpha = 0.7, color = "#800074") +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),
              color = "gray50", se = TRUE, fill = "gray85") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Wald test: " ,p_label.astro),
           hjust = 1.1, vjust = 5.0,
           size = 4) +
  scale_y_continuous(breaks = c(0,1), labels = c("Non-Hypermutant", "Hypermutant")) +
  labs(
    x = "Temozolomide cycles",
    y = ""  ) +
  theme_linedraw()

p.oligo <- summary(glm(HM_group ~ TMZ_cycles,
                       data = oligo.dose.hypermut,
                       family = binomial))$coefficients["TMZ_cycles", "Pr(>|z|)"]

p_label.oligo <- ifelse(p.oligo < 0.001, "p < 0.001", paste0("p = ", signif(p.oligo, 3)))

p4 <- ggplot(oligo.dose.hypermut, aes(x = TMZ_cycles, y =HM_group)) +
  geom_point(width = 0.3, height = 0.05, alpha = 0.7, color = "#298C8C") +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),
              color = "gray50", se = TRUE, fill = "gray85") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Wald test: " ,p_label.oligo),
           hjust = 1.1, vjust = 5.0,
           size = 4) +
  scale_y_continuous(breaks = c(0,1), labels = c("Non-Hypermutant", "Hypermutant")) +
  labs(
    x = "Temozolomide cycles",
    y = ""  ) +
  theme_linedraw()

p2 + p1 + 
  plot_annotation() +
  plot_layout(widths = c(2, 2))&
  theme(plot.tag = element_text(face = 'bold'),
        plot.title = element_text(face = 'bold', size =12, hjust = 0.5))

p3 + p4 + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(widths = c(2))&
  theme(plot.tag = element_text(face = 'bold'),
        plot.title = element_text(face = 'bold', size =12, hjust = 0.5))