# GLASS-I: Cohort Description
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(RPostgres)
library(tidyverse)

# load data
molecular.data <- read.delim("~/driver_changes.txt")%>%                       ## created with GLASS-I_Drivers_Table.R
  select(case_barcode, HM_group) 

con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")             ## data available on Synapse (Synapse ID: syn17038081)

surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
  filter(idh_codel_subtype != "IDHwt")%>%
  filter(!is.na(surgical_interval_mo))

full.surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
  filter(idh_codel_subtype != "IDHwt")

survival.data <- dbGetQuery(con, "SELECT * FROM clinical.cases")[,c(2,6:8)]%>%
  semi_join(surgeries, by="case_barcode")%>%
  filter(!is.na(case_overall_survival_mo))

full.survival.data <- dbGetQuery(con, "SELECT * FROM clinical.cases")[,c(-1)]

dbDisconnect(con)

surgeries <- surgeries%>%
  semi_join(survival.data, by="case_barcode")

full.data <-  full.survival.data%>%
  right_join(full.surgeries, by ="case_barcode")%>%
  mutate(case_source =         if_else(case_source == "19"|case_source == "FG", "CW", 
                               if_else(case_source == "06"|case_source == "DU", "HF",case_source)),
         institution_country = if_else(case_source == "14"| case_source == "CU"| case_source == "CW" | case_source == "DF"| case_source == "DH"|
                                       case_source == "HF"| case_source == "JX"| case_source == "MA" | case_source == "MD"| case_source == "MG"| 
                                       case_source == "NW"| case_source == "PT"| case_source == "SF" | case_source == "SJ", "USA",
                               if_else(case_source == "AT", "Austria",
                               if_else(case_source == "DK", "Germany",
                               if_else(case_source == "ER", "Netherlands",
                               if_else(case_source == "HK", "China",
                               if_else(case_source == "LU"|case_source == "TM", "UK",
                               if_else(case_source == "LX", "Luxembourg",
                               if_else(case_source == "NS"|case_source == "ON", "Australia",
                               if_else(case_source == "PS", "France",
                               if_else(case_source == "SN" |case_source == "SM", "South-Korea",
                               if_else(case_source == "SU"|case_source == "TK", "Japan",
                               if_else(case_source == "TO", "Canada",
                               if_else(case_source == "TQ", "Brazil",
                               NA))))))))))))),
         tumor_type =          if_else(idh_codel_subtype == "IDHmut-codel", "IDH-mutant Oligodendroglioma", "IDH-mutant Astrocytoma"),
         institution_city =    if_else(case_source == "CU", 'New York',
                               if_else(case_source == "CW", 'Cleveland',
                               if_else(case_source == "DH", 'Gainesville',
                               if_else(case_source == "DK", 'Heidelberg',
                               if_else(case_source == "HF", 'Detroit',
                               if_else(case_source == "HK", 'Hong Kong',
                               if_else(case_source == "JX", 'Farmington',
                               if_else(case_source == "LX", 'Luxembourg',
                               if_else(case_source == "MA", 'Rochester', 
                               if_else(case_source == "MD", 'Houston',
                               if_else(case_source == "NW", 'Chicago',
                               if_else(case_source == "ON", 'Melbourne',
                               if_else(case_source == "PS", 'Paris',
                               if_else(case_source == "PT", 'Durham',
                               if_else(case_source == "SF", 'San Francisco',
                               if_else(case_source == "SJ", 'Phoenix', 
                               if_else(case_source == "SM" | case_source == "SN", 'Seoul',
                               if_else(case_source == "SU", 'Tokyo',
                               if_else(case_source == "TM", 'Cardiff',
                               if_else(case_source == "TQ", 'Sao Paolo',
                               NA))))))))))))))))))))
  )%>%
  distinct(case_barcode, surgery_number, .keep_all = TRUE)

patient.lvl.data <- survival.data%>%
  left_join(surgeries, by ="case_barcode")%>%
  left_join(molecular.data, by ="case_barcode")%>%
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
         TMZ=                             case_when(treatment_tmz == T  &  (treatment_concurrent_tmz == T|treatment_concurrent_tmz == F|is.na(treatment_concurrent_tmz))~ "Treated",
                                                    treatment_concurrent_tmz == T  &  (treatment_tmz == T |treatment_tmz == F|is.na(treatment_tmz))~ "Treated",
                                                    treatment_tmz == F  &  (treatment_concurrent_tmz == F |  is.na(treatment_concurrent_tmz))~ "Non-treated",
                                                    is.na(treatment_tmz)  &  treatment_concurrent_tmz == F~ "Non-treated"),
         therapy_other =                  if_else(str_detect(treatment_chemotherapy_other, "guanine"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "Accutane"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "AG881"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "umab"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "tecan"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "tretino"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "Other"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "Hydroxyurea"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "Tyrosine"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "Etoposide"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "olimus"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "inib"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "parib"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "Tamoxifen"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "idomide"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "Mab"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "ICLC"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "RESIST"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "Celecoxib"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "enib"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "Capecitabine"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "Yes"),1,
                               0)))))))))))))))))))))),
         retinol_based =                  if_else(str_detect(treatment_chemotherapy_other, "Accutane"),1,
                                          if_else(str_detect(treatment_chemotherapy_other, "Isotretinoin"),1,
                               0)),
         anti_vegf =                      if_else(str_detect(treatment_chemotherapy_other, "Bevacizumab"),1,0)
  )%>%
  distinct(case_barcode, surgery_number, .keep_all = TRUE)

### data for text manuscript###
# number of institutions participating in manuscript
length(unique(full.data$case_source))

# number of countries participating in manuscript
length(unique(full.data$institution_country))
unique(full.data$institution_country)

# number of clinical and/or molecular patients and surgeries regardless of QC
nlevels(as.factor(full.data$case_barcode))
table(full.data$idh_codel_subtype)

# number of clinical patients and surgeries
nlevels(as.factor(patient.lvl.data$case_barcode))
table(patient.lvl.data$idh_codel_subtype)

# percentage of resections
as.data.frame(table(patient.lvl.data$Surgery))[2,2] / (as.data.frame(table(patient.lvl.data$Surgery))[1,2] +as.data.frame(table(patient.lvl.data$Surgery))[2,2]) *100

#extra data for text
text.data <- patient.lvl.data%>%
  group_by(case_barcode)%>%
  summarise(idh_codel_subtype = idh_codel_subtype[1],
            RT = if_else(any(Radiotherapy == "Treated"), 1, 0),
            TMZ = if_else(any(TMZ == "Treated"), 1, 0),
            ALK = if_else(any(`Alkylating chemotherapy` == "Treated"), 1, 0),
            PLAT = if_else(any(str_detect(treatment_chemotherapy_other, "platin")), 1, 0),
            OTHER = if_else(any(therapy_other ==1),1,0),
            RETINOL = if_else(any(retinol_based ==1),1,0),
            ANTIVEGF = if_else(any(anti_vegf ==1),1,0)
            )%>%
  ungroup()%>%  
  replace(is.na(.), 0) %>%
  mutate(no_treatment = if_else(RT ==0 & ALK ==0 & PLAT == 0 & OTHER ==0, 1, 0))

text.data.astros <- text.data%>%
  filter(idh_codel_subtype == "IDHmut-noncodel")

text.data.oligos <- text.data%>%
  filter(idh_codel_subtype == "IDHmut-codel")

# percentage of patients which received alkylating chemotherapy
sum(text.data$ALK)/nrow(text.data)*100 

# percentage of patients which received temozolomide
sum(text.data.astros$TMZ)/nrow(text.data.astros)*100 
sum(text.data.oligos$TMZ)/nrow(text.data.oligos)*100 

# percentage of patients which received platinum-based chemotherapy
sum(text.data.astros$PLAT)/nrow(text.data.astros)*100 
sum(text.data.oligos$PLAT)/nrow(text.data.oligos)*100 

# numbers of patients which received platinum-based chemotherapy compared to alkylating chemotherapy
table(text.data$PLAT, text.data$ALK)

# percentage of patients which received radiotherapy
sum(text.data.astros$RT)/nrow(text.data.astros)*100 
sum(text.data.oligos$RT)/nrow(text.data.oligos)*100 

# percentage of patients which received other therapies
sum(text.data$OTHER)/nrow(text.data)*100 

# percentage of patients which received anti-VEGF therapy
sum(text.data$ANTIVEGF)/nrow(text.data)*100 

# percentage of patients which received retinol-based drugs
sum(text.data$RETINOL)/nrow(text.data)*100 

# percentage of patients which were not reported to receive therapy
sum(text.data.astros$no_treatment)/nrow(text.data.astros)*100 
sum(text.data.oligos$no_treatment)/nrow(text.data.oligos)*100 

# Map data
map.data <- full.data%>%
  filter(!duplicated(case_barcode))

table(map.data$institution_city, map.data$tumor_type)

# Summary Figure data
rt.tmz.data <- patient.lvl.data%>%
  group_by(case_barcode)%>%
  summarise(idh_codel_subtype = idh_codel_subtype[1],
            HM = HM_group[1],
            RT = if_else(any(Radiotherapy == "Treated"), 1, 0),
            TMZ = if_else(any(TMZ == "Treated"), 1, 0)
  )%>%
  ungroup()%>%  
  mutate(rt_tmz = if_else(RT == 1 & TMZ ==1, 1, 0),
         no_rt_or_tmz = if_else(RT == 0 & TMZ ==0, 1, 0),
         only_rt = if_else(RT == 1 & TMZ ==0, 1, 0), 
         only_tmz = if_else(RT == 0 & TMZ ==1, 1, 0))

rt.tmz.data.astros <- rt.tmz.data%>%
  filter(idh_codel_subtype == "IDHmut-noncodel")

table(rt.tmz.data.astros$RT)
table(rt.tmz.data.astros$TMZ)
table(rt.tmz.data.astros$rt_tmz)
table(rt.tmz.data.astros$no_rt_or_tmz)
table(rt.tmz.data.astros$only_rt)
table(rt.tmz.data.astros$only_tmz)
table(rt.tmz.data.astros$TMZ , rt.tmz.data.astros$HM)

rt.tmz.data.oligos <- rt.tmz.data%>%
  filter(idh_codel_subtype == "IDHmut-codel")

table(rt.tmz.data.oligos$RT)
table(rt.tmz.data.oligos$TMZ)
table(rt.tmz.data.oligos$rt_tmz)
table(rt.tmz.data.oligos$no_rt_or_tmz)
table(rt.tmz.data.oligos$only_rt)
table(rt.tmz.data.oligos$only_tmz)
table(rt.tmz.data.oligos$TMZ , rt.tmz.data.oligos$HM)