# GLASS-I: TCGA OS Comparison
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(tidyverse)
library(RPostgres)
library(survival)
library(survminer)
library(patchwork)
library(gtsummary)

# load data
tcga.data <- read.delim("~/data_clinical_patient.txt", skip = 4)%>%           ## data available on TCGA (lgg_tcga_pan_can_atlas_2018)
  filter(!str_detect(SUBTYPE, "wt"))%>%
  filter(!duplicated(PATIENT_ID))%>%
  filter(SUBTYPE !="")%>%
  mutate(case_barcode = PATIENT_ID,
         tumor_type = if_else(SUBTYPE == "LGG_IDHmut-non-codel", "IDH-mutant Astrocytoma",
                              if_else(SUBTYPE == "LGG_IDHmut-codel", "IDH-mutant Oligodendroglioma", NA)),
         age_at_diagnosis = AGE,
         os_from_diagnosis = OS_MONTHS/12,
         `Follow-up (years)` = os_from_diagnosis,
         os_status = if_else(OS_STATUS == "0:LIVING", 0, 1),
         os_from_birth = (AGE+(OS_MONTHS/12)),
         Cohort = "TCGA")%>%
  select(case_barcode,tumor_type, age_at_diagnosis,`Follow-up (years)`, os_from_birth, os_from_diagnosis,os_status, Cohort)%>%
  filter(!is.na(age_at_diagnosis))

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")             ## data available on Synapse (Synapse ID: syn17038081)

surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
  filter(idh_codel_subtype != "IDHwt")%>%
  filter(!is.na(surgical_interval_mo))

survival.data <- dbGetQuery(con, "SELECT * FROM clinical.cases")[,c(2,6:8)]%>%
  semi_join(surgeries, by="case_barcode")%>%
  filter(!is.na(case_overall_survival_mo))

dbDisconnect(con)

glass.data <- survival.data %>%
  left_join(surgeries, by = "case_barcode") %>%
  filter(!duplicated(case_barcode))%>%
  mutate(tumor_type = if_else(idh_codel_subtype == "IDHmut-noncodel", "IDH-mutant Astrocytoma",
                              if_else(idh_codel_subtype == "IDHmut-codel", "IDH-mutant Oligodendroglioma", NA)),
         age_at_diagnosis = case_age_diagnosis_years,
         os_from_diagnosis = case_overall_survival_mo/12,
         `Follow-up (years)` = os_from_diagnosis,
         os_status = if_else(case_vital_status == "alive", 0, 1),
         os_from_birth = case_age_diagnosis_years+(case_overall_survival_mo/12),
         Cohort = "GLASS")%>%
  select(case_barcode,tumor_type, age_at_diagnosis,`Follow-up (years)`, os_from_birth, os_from_diagnosis,os_status, Cohort)

# comparison table
full.data <- rbind(tcga.data,glass.data)

medians <- full.data %>%
  group_by(Cohort, tumor_type) %>%
  summarise(median_age = median(age_at_diagnosis, na.rm = TRUE))%>%
  ungroup()

# Define a function to perform the Wilcoxon test and return the p-value
get_wilcox_pvalue <- function(df, tumor) {
  df.subset <- df %>% filter(tumor_type == tumor)
  wilcox.test(age_at_diagnosis ~ Cohort, data = df.subset)$p.value
}

# Calculate and annotate p-values for each cohort and each tumor type
p_values <- full.data %>%
  group_by(tumor_type) %>%
  summarise(p_value = get_wilcox_pvalue(full.data, tumor_type))%>%
  ungroup()%>%
  mutate(p_value_formatted = if_else(p_value < 0.001, "p < 0.001", paste0("p = ", round(p_value, 3))))

p1<- ggplot(full.data)+
  geom_histogram(aes(x=age_at_diagnosis, color = Cohort, fill = Cohort), size =1, bins = 40,  alpha = 0.2, position="identity")+
  theme_linedraw()+
  scale_color_manual(values=c("GLASS"="#f1a226","TCGA"= "#298c8c"))+
  scale_fill_manual(values=c("GLASS"="#f1a226","TCGA"= "#298c8c"))+
  facet_wrap(~tumor_type)+
  xlab("Age at diagnosis (years)")+
  ylab("Patients (n)") + 
  labs(fill='IDH-mutant Glioma Patient Cohort', color='IDH-mutant Glioma Patient Cohort')+
  theme(legend.position = "top",
        legend.title = element_text(size =13),
        legend.text = element_text(size = 13),
        strip.text = element_text(face = "bold", size =13)) + 
  scale_x_continuous(breaks=seq(0, 80, 10), limits = c(0,80))+
  geom_vline(data = medians, aes(xintercept = median_age, color = Cohort),
             linetype = "dashed", size = 1) +
  geom_text(data = p_values, aes(x = 60, y = 21, label = p_value_formatted))

survival.data.os <- full.data%>%
  mutate(survival_type = "a_os_from_diagnosis")

survival.data.birth <- full.data%>%
  mutate(survival_type = "b_os_from_birth",
         os_from_diagnosis = os_from_birth)

survival.data.full <- rbind(survival.data.os, survival.data.birth)

fit <- survfit( Surv(os_from_diagnosis, os_status) ~ Cohort, data = survival.data.full )
p2 <- ggsurvplot(fit, 
           survival.data.full, 
           facet.by = c( "survival_type", "tumor_type") ,
           palette = c("#f1a226", "#298c8c"), 
           ggtheme = theme_linedraw(),
           short.panel.labs = T,
           legend.title = "IDH-mutant Glioma Patient Cohort",
           legend.labs = c("GLASS", "TCGA"),
           pval.coord = c(62.5, 0.95),
           surv.median.line = "v",
           break.x.by = 10,
           font.legend = c(13, "plain", "black"),
           censor.shape = "|",
           censor.size = 4,
           panel.labs = list(survival_type = c("From Diagnosis", "From Birth")),
           panel.labs.font = list(face = "bold", size =13),
           xlab ="Overall survival (years)",
           pval = T)

p3 <- ggplot() + theme_void()

p1/p2/p3 +
plot_annotation(tag_levels = 'A')+
plot_layout(heights = c(1, 2,2)) &
theme(plot.tag = element_text(face = 'bold', size =16))

oligo.data <- full.data%>%
  filter(tumor_type == "IDH-mutant Oligodendroglioma")

oligo.table <- oligo.data %>%
  tbl_summary(include = c(
    `Follow-up (years)`,
    Cohort
  ),
  by = Cohort,
  missing = "always",
  type = all_continuous() ~ "continuous2",
  statistic = list(
    all_continuous() ~ c(
      "{median} [{p25}, {p75}]"
    ),
    all_categorical() ~ "{p}% ({n}/{N})"
  ))%>% 
  modify_header(label ~ "") %>%
  add_p()%>%
  italicize_labels()%>%
  modify_caption("")

survfit(Surv(os_from_diagnosis, os_status) ~ Cohort, data = oligo.data)
summary(coxph(Surv(os_from_diagnosis, os_status) ~ Cohort, data = oligo.data))

survfit(Surv(os_from_birth, os_status) ~ Cohort, data = oligo.data)
summary(coxph(Surv(os_from_birth, os_status) ~ Cohort, data = oligo.data))

astro.data <- full.data%>%
  filter(tumor_type == "IDH-mutant Astrocytoma")

astro.table <- astro.data %>%
  tbl_summary(include = c(
    `Follow-up (years)`,
    Cohort
  ),
  by = Cohort,
  missing = "always",
  type = all_continuous() ~ "continuous2",
  statistic = list(
    all_continuous() ~ c(
      "{median} [{p25}, {p75}]"
    ),
    all_categorical() ~ "{p}% ({n}/{N})"
  ))%>% 
  modify_header(label ~ "") %>%
  add_p()%>%
  italicize_labels()%>%
  modify_caption("")

survfit(Surv(os_from_diagnosis, os_status) ~ Cohort, data = astro.data)
summary(coxph(Surv(os_from_diagnosis, os_status) ~ Cohort, data = astro.data))

survfit(Surv(os_from_birth, os_status) ~ Cohort, data = astro.data)
summary(coxph(Surv(os_from_birth, os_status) ~ Cohort, data = astro.data))

tbl_merge(
  tbls = list(oligo.table, astro.table),
  tab_spanner = c("**IDH-mutant Oligodendroglioma**", "**IDH-mutant Astrocytoma**")
)