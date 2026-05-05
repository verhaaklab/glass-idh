#################################################################################################################
# Content: Baseline table of patients with IDH-mutant gliomas
# Author: Mircea Tesileanu
# Date: 2025.07.03
# R-version: 4.3.2
# Tables: Table 1
#################################################################################################################

# Clean start and load libraries
rm(list = ls())
gc()

library(tidyverse)    # v2.0.0
library(RPostgres)    # v1.4.7
library(gtsummary)    # v2.0.4
library(gt)           # v0.11.1

#################################################################################################################

# Load relevant data
molecular.data <-  read.delim("/path/to/driver_changes_14042025.txt")%>%
  filter(gene_symbol == "CDKN2A/CDKN2B") %>%
  mutate(CDKN2AB =          if_else(!str_detect(tumor_pair_barcode, "-TP-"), NA,
                            if_else(driver_status == "No_AMP/HD", "absent", 
                            if_else(driver_change == "R", "absent", "present"))),
         timing.recurrent = substr(tumor_pair_barcode, 20,21),
         timing.primary =   substr(tumor_pair_barcode, 14,15) )%>%
  select(-idh_codel_subtype) 

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")

surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
  mutate(sample_barcode =   ifelse(surgery_number == 1, paste0(case_barcode, "-TP"),
                            ifelse(surgery_number == 2, paste0(case_barcode, "-R1"),
                            ifelse(surgery_number == 3, paste0(case_barcode, "-R2"),
                            ifelse(surgery_number == 4, paste0(case_barcode, "-R3"),
                            ifelse(surgery_number == 5, paste0(case_barcode, "-R4"),
                            ifelse(surgery_number == 6, paste0(case_barcode, "-R5"),
                            ifelse(surgery_number == 7, paste0(case_barcode, "-R6"),
                            ifelse(surgery_number == 8, paste0(case_barcode, "-R7"),
                            NA)))))))),
         timing.surgery = substr(sample_barcode, 14,15)) %>%
  filter(idh_codel_subtype != "IDHwt")%>%  
  filter(case_barcode != "GLSS-DK-0015")%>% # missing clinical data
  filter(case_barcode != "OLIG-PT-0030")    # missing survival data

silver.set <- dbGetQuery(con, "SELECT * FROM analysis.silver_set")%>%
  semi_join(surgeries, by = "case_barcode")

survival.data <- dbGetQuery(con, "SELECT * FROM clinical.cases")[,c(2,5:8)]%>%
  semi_join(silver.set, by="case_barcode")

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
  filter(timing.surgery =="TP")%>%
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
         `Histological characteristics` = histology,
         Grade =                          case_when(CDKN2AB == "present"& idh_codel_subtype == "IDHmut-noncodel" & (grade == "II"|grade == "III"|grade == "IV"|is.na(grade)) ~ 4,
                                                    (CDKN2AB == "absent"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-noncodel" & grade == "II" ~ 2,
                                                    (CDKN2AB == "absent"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-noncodel" & grade == "III" ~ 3,
                                                    (CDKN2AB == "absent"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-noncodel" & grade == "IV" ~ 4,
                                                    (CDKN2AB == "absent"|CDKN2AB == "present"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-codel" & grade == "II" ~ 2,
                                                    (CDKN2AB == "absent"|CDKN2AB == "present"|is.na(CDKN2AB))& idh_codel_subtype == "IDHmut-codel" & grade == "III" ~ 3)
        )%>%
  select(case_barcode, Surgery, Radiotherapy, `Alkylating chemotherapy`, `Platinum-based chemotherapy`, `Other therapies`, `Histological characteristics`, Grade)


#################################################################################################################

# Create baseline table and include statistical tests to compare patients with IDH-mutant oligodendrogliomas with those with IDH-mutant astrocytomas
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
                          `Histological characteristics`,
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
  add_p()%>%
  italicize_labels()%>%
  modify_caption("")%>%
  as_gt()%>%
  fmt_markdown(columns = c(2,3))


baseline.table <- patient.lvl.table %>%
  tab_row_group(
    label = md("**Primary Pathology**"),
    rows = 31:41
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

#################################################################################################################