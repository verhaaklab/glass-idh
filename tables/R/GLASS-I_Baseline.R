# GLASS-I: Baseline Table
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(RPostgres)
library(tidyverse)
library(gtsummary)
library(gt)

# load data
molecular.data <- read.delim("~/driver_changes.txt")%>%                       ## created with GLASS-I_Drivers_Table.R
  filter(gene_symbol == "CDKN2A/CDKN2B") %>%
  mutate(CDKN2AB =          if_else(!str_detect(tumor_pair_barcode, "-TP-"), NA,
                            if_else(driver_status == "No_AMP/HD", "absent", 
                            if_else(driver_change == "R", "absent", "present"))),
         timing.recurrent = substr(tumor_pair_barcode, 20,21),
         timing.primary =   substr(tumor_pair_barcode, 14,15) )%>%
  select(-idh_codel_subtype) 

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")             ## data available on Synapse (Synapse ID: syn17038081)

surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
  filter(idh_codel_subtype != "IDHwt")%>%
  filter(!is.na(surgical_interval_mo))

survival.data <- dbGetQuery(con, "SELECT * FROM clinical.cases")[,c(2,5:8)]%>%
  semi_join(surgeries, by="case_barcode")%>%
  filter(!is.na(case_overall_survival_mo))

dbDisconnect(con)

patient.lvl.data.demographics <- survival.data%>%
  left_join(surgeries, by ="case_barcode")%>%
  left_join(molecular.data, by ="case_barcode")%>%
  mutate(`Age at diagnosis (years)` =     case_age_diagnosis_years,
         Sex =                            if_else(case_sex == "female", "Female", "Male"),
         `Follow-up (years)` =            case_overall_survival_mo/12,
         glioma.type =                    if_else(idh_codel_subtype != "IDHmut-codel", "IDH-mutant<br>Astrocytomas", "IDH-mutant<br>Oligodendrogliomas"))%>%
  select(case_barcode,`Age at diagnosis (years)`, Sex,`Follow-up (years)`, glioma.type)%>%
  distinct()

patient.lvl.data.primary <- survival.data%>%
  left_join(surgeries, by ="case_barcode")%>%
  left_join(molecular.data, by ="case_barcode")%>%
  filter(surgery_number ==1)%>%
  mutate(Surgery =                        case_when(surgery_type == "Biopsy" & (surgery_extent_of_resection == "Biopsy"|is.na(surgery_extent_of_resection)) ~"Biopsy",
                                                    surgery_type == "Craniotomy" & (surgery_extent_of_resection == "Subtotal"|surgery_extent_of_resection == "Total"|is.na(surgery_extent_of_resection)) ~ "Resection",
                                                    surgery_type == "Craniotomy" & surgery_extent_of_resection == "Biopsy"  ~ "Biopsy",
                                                    is.na(surgery_type) & surgery_extent_of_resection == "Biopsy" ~ "Biopsy",
                                                    is.na(surgery_type) & (surgery_extent_of_resection == "Subtotal"|surgery_extent_of_resection == "Total") ~ "Resection"),
         Radiotherapy =                   case_when(treatment_radiotherapy == T ~ "Treated",
                                                    treatment_radiotherapy == F ~"Non-treated"),
         `Alkylating chemotherapy`=       case_when(treatment_alkylating_agent == T  &  (treatment_concurrent_tmz == T|treatment_concurrent_tmz == F|is.na(treatment_concurrent_tmz))~ "Treated",
                                                    treatment_concurrent_tmz == T  &  (treatment_alkylating_agent == T |treatment_alkylating_agent == F|is.na(treatment_alkylating_agent))~ "Treated",
                                                    treatment_alkylating_agent == F  &  (treatment_concurrent_tmz == F |  is.na(treatment_concurrent_tmz))~ "Non-treated",
                                                    is.na(treatment_alkylating_agent)  &  treatment_concurrent_tmz == F~ "Non-treated"),
         `Platinum-based chemotherapy`=   if_else(str_detect(treatment_chemotherapy_other, "platin"), "Treated", "Non-treated"),
         `Other therapies` =              if_else(str_detect(treatment_chemotherapy_other, "guanine"), "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "Accutane"), "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "AG881"), "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "umab"), "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "tecan"), "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "tretino"), "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "Other"), "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "Tyrosine"), "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "Hydroxyurea"),  "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "Etoposide"),  "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "parib"),  "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "Mab"), "Treated",                        
                                          if_else(str_detect(treatment_chemotherapy_other, "olimus"), "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "inib"), "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "Tamoxifen"), "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "idomide"), "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "Celecoxib"), "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "ICLC"), "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "enib"), "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "Yes"), "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "Capecitabine"), "Treated",
                                          if_else(str_detect(treatment_chemotherapy_other, "RESIST"), "Treated",
                                          "Non-treated")))))))))))))))))))))),
         Grade =                          case_when(CDKN2AB == "present"& idh_codel_subtype == "IDHmut-noncodel" & (grade == "II"|grade == "III"|grade == "IV"|is.na(grade)) ~ 4,
                                                    (CDKN2AB == "absent"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-noncodel" & grade == "II" ~ 2,
                                                    (CDKN2AB == "absent"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-noncodel" & grade == "III" ~ 3,
                                                    (CDKN2AB == "absent"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-noncodel" & grade == "IV" ~ 4,
                                                    (CDKN2AB == "absent"|CDKN2AB == "present"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-codel" & grade == "II" ~ 2,
                                                    (CDKN2AB == "absent"|CDKN2AB == "present"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-codel" & grade == "III" ~ 3)
        )%>%
  select(case_barcode, Surgery, Radiotherapy, `Alkylating chemotherapy`, `Platinum-based chemotherapy`, `Other therapies`, Grade)

patient.lvl.data <- patient.lvl.data.demographics%>%
left_join(patient.lvl.data.primary, by = "case_barcode")

patient.lvl.table <- patient.lvl.data %>%
  tbl_summary(include = c(
                          `Age at diagnosis (years)`,
                          Sex,
                          `Follow-up (years)`,
                          Surgery,
                          Radiotherapy,
                          `Alkylating chemotherapy`,
                          `Platinum-based chemotherapy`,
                          `Other therapies`,
                          Grade,
                          glioma.type
  ),
  by = glioma.type,
  missing = "always",
  type = all_continuous() ~ "continuous2",
  statistic = list(
    all_continuous() ~ c(
      "{median} [{min}, {max}]"
    ),
    all_categorical() ~ "{p}% ({n}/{N})"
  ))%>% 
  modify_header(label ~ "") %>%
  add_p() %>%
  italicize_labels()%>%
  modify_caption("")%>%
  as_gt()%>%
  fmt_markdown(columns = c(2,3))

# Table 1: Baseline Table
Table.1 <- patient.lvl.table %>%
  tab_row_group(
    label = md("**Primary Pathology**"),
    rows = 31:35
  ) %>%
  tab_row_group(
    label = md("**Primary Treatment**"),
    rows = 11:30
  ) %>%
  tab_row_group(
    label = md("**Demographics**"),
    rows = 1:10
  ) %>%
  tab_style(
    style = cell_fill(color = "gray97"),
    locations = cells_row_groups()
  ) %>%
  sub_values(values = c("0% (0/73)"), replacement = "")%>%
  tab_options(  heading.padding = 3,  data_row.padding.horizontal = 30, row_group.padding = 5, table.font.size = 15, data_row.padding = 3, footnotes.padding = 3)

Table.1

table(patient.lvl.data.primary$Radiotherapy, patient.lvl.data.primary$`Alkylating chemotherapy`)

astrocytomas <- patient.lvl.data%>%
  filter(glioma.type == "IDH-mutant<br>Astrocytomas")%>%
  mutate(Treatment = if_else(Radiotherapy == "Treated" & `Alkylating chemotherapy` == "Treated", "Chemoradiation",
                             if_else(Radiotherapy == "Treated" & `Alkylating chemotherapy` == "Non-treated", "Radiotherapy only",
                                     if_else(Radiotherapy == "Non-treated" & `Alkylating chemotherapy` == "Treated", "Chemotherapy only",
                                             if_else(Radiotherapy == "Non-treated" & `Alkylating chemotherapy` == "Non-treated", "Non-treated",
                                                     NA)))))

oligodendrogliomas <- patient.lvl.data%>%
  filter(glioma.type == "IDH-mutant<br>Oligodendrogliomas")%>%
  mutate(Treatment = if_else(Radiotherapy == "Treated" & `Alkylating chemotherapy` == "Treated", "Chemoradiation",
                             if_else(Radiotherapy == "Treated" & `Alkylating chemotherapy` == "Non-treated", "Radiotherapy only",
                                     if_else(Radiotherapy == "Non-treated" & `Alkylating chemotherapy` == "Treated", "Chemotherapy only",
                                             if_else(Radiotherapy == "Non-treated" & `Alkylating chemotherapy` == "Non-treated", "Non-treated",
                                                     NA)))))%>%
  filter(!is.na(Grade))

# Primary treatment vs grade at primary surgery: astrocytomas
astrocytomas %>%
  tbl_summary(include = c(
    Treatment,
    Grade
  ),
  by = Grade,
  missing = "always",
  type = all_continuous() ~ "continuous2",
  statistic = list(
    all_categorical() ~ "{p}% ({n}/{N})"
  ))%>% 
  add_p()

# Primary treatment vs grade at primary surgery: oligodendrogliomas
oligodendrogliomas %>%
  tbl_summary(include = c(
    Treatment,
    Grade
  ),
  by = Grade,
  missing = "always",
  type = all_continuous() ~ "continuous2",
  statistic = list(
    all_categorical() ~ "{p}% ({n}/{N})"
  ))%>% 
  add_p()