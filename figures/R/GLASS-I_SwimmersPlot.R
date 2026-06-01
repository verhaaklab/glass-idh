# GLASS-I: Swimmers Plot
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(RPostgres)
library(tidyverse)
library(patchwork)
library(survival)
library(forcats)

# load data
con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")             ## data available on Synapse (Synapse ID: syn17038081)

surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
  filter(idh_codel_subtype != "IDHwt")%>%
  filter(!is.na(surgical_interval_mo))

survival.data <- dbGetQuery(con, "SELECT * FROM clinical.cases")[,c(2,5:8)]%>%
  semi_join(surgeries, by="case_barcode")%>%
  filter(!is.na(case_overall_survival_mo))

dbDisconnect(con)

surgeries <- surgeries%>%
  semi_join(survival.data, by="case_barcode")

long.survival <- survival.data %>%
  mutate(alive_at_fu = if_else(case_vital_status == "alive", 1, 0),
         time = case_overall_survival_mo/12,
         follow_up_time = case_overall_survival_mo/12) %>%
  select(-case_age_diagnosis_years, -case_vital_status, -case_overall_survival_mo)

long.age <- survival.data %>%
  mutate(age_under_25 = if_else(case_age_diagnosis_years <25, 1, 0),
         age_25_to_39 = if_else(case_age_diagnosis_years >24 & case_age_diagnosis_years <40, 1, 0),
         age_40_to_55 = if_else(case_age_diagnosis_years >39 & case_age_diagnosis_years <56, 1, 0),
         age_above_55 = if_else(case_age_diagnosis_years >55, 1, 0),
         time = -0.5)%>%
  select(-case_age_diagnosis_years, -case_vital_status, -case_overall_survival_mo)

long.surgeries <- surgeries%>%
  mutate(surgery_any = 1,
         time = surgical_interval_mo/12) %>%
  select(case_barcode, time, surgery_any)

missing.first.surgeries <- long.surgeries %>%
  group_by(case_barcode) %>%
  filter(!any(time == 0 & surgery_any == 1)) %>%
  distinct(case_barcode)%>%
  ungroup() %>%
  mutate(time = 0, surgery_any = 1)

long.surgeries <- long.surgeries %>%
  bind_rows(missing.first.surgeries) 

long.radiotherapy <- surgeries%>%
  mutate(radiotherapy = if_else(treatment_radiotherapy == T, 1,0),
         time = (surgical_interval_mo/12)+0.25) %>%
  select(case_barcode, time, radiotherapy)

short.surgeries <- surgeries%>%
  mutate(tumor_type = if_else(idh_codel_subtype == "IDHmut-codel", "IDH-mutant Oligodendroglioma", "IDH-mutant Astrocytoma")) %>%
  select(case_barcode, tumor_type)

long.chemotherapy.alk <- surgeries%>%
  mutate(tmz = if_else(treatment_tmz==T, 1, 0),
         alkylating_other = case_when(treatment_tmz == T & (treatment_concurrent_tmz == T| treatment_concurrent_tmz == F|is.na(treatment_concurrent_tmz)) & (treatment_alkylating_agent ==T|treatment_alkylating_agent ==F|is.na(treatment_alkylating_agent)) ~ 0,
                                      treatment_concurrent_tmz == T & (treatment_tmz == T| treatment_tmz == F|is.na(treatment_tmz)) & (treatment_alkylating_agent ==T|treatment_alkylating_agent ==F|is.na(treatment_alkylating_agent)) ~ 0,
                                      treatment_alkylating_agent == T & (treatment_concurrent_tmz == F|is.na(treatment_concurrent_tmz)) &  (treatment_tmz == F|is.na(treatment_tmz)) ~ 1
                                      ),
         time = (surgical_interval_mo/12)+0.6) %>%
  select(case_barcode, time, tmz, alkylating_other)

long.tmz.concurrent <- surgeries%>%
  mutate(tmz_concurrent = if_else(treatment_concurrent_tmz==T, 1, 0),
         time = (surgical_interval_mo/12)+0.25) %>%
  select(case_barcode, time, tmz_concurrent)

long.chemotherapy.platinum <- surgeries%>%
  mutate(platinum_based = if_else(str_detect(treatment_chemotherapy_other, "platin"),1, 0 ),
         time = (surgical_interval_mo/12)+1) %>%
  select(case_barcode, time, platinum_based)

long.therapy.other <- surgeries%>%
  mutate(therapy_other = if_else(str_detect(treatment_chemotherapy_other, "guanine"),1,
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
                         if_else(str_detect(treatment_chemotherapy_other, "parib"),1,
                         if_else(str_detect(treatment_chemotherapy_other, "inib"),1,
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
         time = (surgical_interval_mo/12)+1.3) %>%
  select(case_barcode, time, therapy_other)

combined_table <- long.age%>%
  full_join(long.survival, by = c("case_barcode", "time"))%>%
  full_join(long.surgeries, by = c("case_barcode", "time"))%>%
  full_join(long.radiotherapy, by = c("case_barcode", "time"))%>%
  full_join(long.chemotherapy.alk, by = c("case_barcode", "time"))%>%
  full_join(long.tmz.concurrent, by = c("case_barcode", "time"))%>%
  full_join(long.chemotherapy.platinum, by = c("case_barcode", "time"))%>%
  full_join(long.therapy.other, by = c("case_barcode", "time"))%>%
  full_join(short.surgeries, by = c("case_barcode"))%>%
  replace(is.na(.), 0) %>%
  distinct()%>%
  mutate(surgery_any_this_year = case_when(surgery_any == 1 ~ time),
         radiotherapy_this_year = case_when(radiotherapy == 1 ~ time),
         alive_at_fu_this_year = case_when(alive_at_fu == 1 ~ time),
         tmz_this_year = case_when(tmz == 1 ~ time),
         tmz_concurrent_this_year = case_when(tmz_concurrent == 1 ~ time),
         alkylating_other_this_year = case_when(alkylating_other == 1 ~ time),
         platinum_based_this_year = case_when(platinum_based == 1 ~ time),
         therapy_other_this_year = case_when(therapy_other == 1 ~ time),
         age_under_25_this_year = case_when(age_under_25 == 1 ~ time),
         age_25_to_39_this_year = case_when(age_25_to_39 == 1 ~ time),
         age_40_to_55_this_year = case_when(age_40_to_55 == 1 ~ time),
         age_above_55_this_year = case_when(age_above_55 == 1 ~ time))%>%
  group_by(case_barcode)%>%
  mutate(max_time = max(time),
         radiotherapy_this_year = if_else(radiotherapy_this_year < max(follow_up_time), radiotherapy_this_year, max(follow_up_time)),
         tmz_this_year = if_else(tmz_this_year < max(follow_up_time), tmz_this_year, max(follow_up_time)),
         tmz_concurrent_this_year = if_else(tmz_concurrent_this_year < max(follow_up_time), tmz_concurrent_this_year, max(follow_up_time)),
         alkylating_other_this_year = if_else(alkylating_other_this_year < max(follow_up_time), alkylating_other_this_year, max(follow_up_time)),
         platinum_based_this_year = if_else(platinum_based_this_year < max(follow_up_time), platinum_based_this_year, max(follow_up_time)),
         therapy_other_this_year = if_else(therapy_other_this_year < max(follow_up_time), therapy_other_this_year, max(follow_up_time))
         )%>%
  ungroup()%>%
  mutate(case_barcode = fct_reorder(factor(case_barcode), max_time),
         ` ` = 'Alive at last follow-up' )

astrocytomas <- combined_table %>%
  filter(tumor_type == "IDH-mutant Astrocytoma")

astro <- astrocytomas%>%
  filter(follow_up_time!=0)%>%
  mutate(os = if_else(alive_at_fu ==1, 0, 1))

median.os.astro <- unname(unlist(median(survfit(Surv(astro$follow_up_time, astro$os) ~ tumor_type, data = astro))))

oligodendrogliomas <- combined_table %>%
  filter(tumor_type == "IDH-mutant Oligodendroglioma")


oligo <- oligodendrogliomas%>%
  filter(follow_up_time!=0)%>%
  mutate(os = if_else(alive_at_fu ==1, 0, 1))

median.os.oligo <- unname(unlist(median(survfit(Surv(oligo$follow_up_time, oligo$os) ~ tumor_type, data = oligo))))

cols <- c("<25 yo"="gray75",
          "25-39 yo"= "gray50",
          "40-55 yo"= "gray25",
          ">55 yo"= "gray0",
          "Temozolomide (Adjuvant)" = "#ae282c",
          "Radiotherapy + Temozolomide (Concurrent)" = "#ae282c",
          "Alkylating Chemotherapy (Non-Temozolomide)" = "#ae282c",
          "Platinum-based Chemotherapy" = "#edac2d",
          "Other therapy" = "#007832",
          "Radiotherapy" = "#2066a8",
          "Surgery"= "black"
)

fills <- c("<25 yo"="gray95",
           "25-39 yo"= "gray70",
           "40-55 yo"= "gray45",
           ">55 yo"= "gray20",
           "Temozolomide (Adjuvant)" = "#f6d6c2",
           "Radiotherapy + Temozolomide (Concurrent)" = "#f6d6c2",
           "Alkylating Chemotherapy (Non-Temozolomide)" = "#f6d6c2",
           "Platinum-based Chemotherapy" = "#f0c571",
           "Other therapy" = "#007832",
           "Radiotherapy" = "#cde1ec",
           "Surgery"= "black"
)

p1 <-  astrocytomas %>%
  group_by(case_barcode)%>%
  mutate(max_follow_up_time = max(follow_up_time))%>%
  ungroup()%>%
  mutate(case_barcode = fct_reorder(case_barcode, max_follow_up_time, .desc = T)) %>%
  ggplot(aes(x = case_barcode, y = follow_up_time)) +
  geom_col(col = "#800074", fill = "#a8009c", alpha = 0.05,
           width = 0.6, linewidth = 0.4) +
  labs(y = "", x = "")+
  geom_hline(yintercept = median.os.astro, linetype = "twodash", color = "#a8009c", size = 1) + 
  geom_segment(aes(y=alive_at_fu_this_year, yend=alive_at_fu_this_year+1.5, xend=case_barcode, linetype = ` `),
                                       arrow = arrow(length=unit(.3, 'cm')),
                                       col = "gray10",
                                       lwd=0.5)+
  geom_segment(aes(y=alive_at_fu_this_year, yend=alive_at_fu_this_year+1.5, xend=case_barcode), 
               arrow = arrow(length=unit(.3, 'cm')),
               col = "#800074",
               lwd=0.6)+  
  geom_point(aes(y=radiotherapy_this_year, col = "Radiotherapy", fill ="Radiotherapy"),
             size=2,
             stroke =1,
             shape = 23) +
  geom_point(aes(y=tmz_concurrent_this_year, col = "Radiotherapy + Temozolomide (Concurrent)" , fill ="Radiotherapy + Temozolomide (Concurrent)" ),
             size=2,
             stroke =1,
             shape = 23) +
  geom_point(aes(y=tmz_this_year, col = "Temozolomide (Adjuvant)", fill ="Temozolomide (Adjuvant)"),
             size=2,
             stroke =1,
             shape = 24) +
  geom_point(aes(y=alkylating_other_this_year, col =  "Alkylating Chemotherapy (Non-Temozolomide)", fill = "Alkylating Chemotherapy (Non-Temozolomide)"),
             size=2,
             stroke =1,
             shape = 25) +
  geom_point(aes(y=platinum_based_this_year, col = "Platinum-based Chemotherapy", fill ="Platinum-based Chemotherapy"),
             size=1.5,
             stroke =1,
             shape = 22) +
  geom_point(aes(y=therapy_other_this_year, col = "Other therapy", fill ="Other therapy"),
             size=1.5,
             stroke =1,
             shape = 4) +
  geom_point(aes(y=age_under_25_this_year, col = "<25 yo", fill ="<25 yo"),
             size=2,
             stroke =1,
             shape = 21) +
  geom_point(aes(y=age_25_to_39_this_year, col = "25-39 yo", fill ="25-39 yo"),
             size=2,
             stroke =1,
             shape = 21) +
  geom_point(aes(y=age_40_to_55_this_year, col = "40-55 yo", fill ="40-55 yo"),
             size=2,
             stroke =1,
             shape = 21) +
  geom_point(aes(y=age_above_55_this_year, col = ">55 yo", fill =">55 yo"),
             size=2,
             stroke =1,
             shape = 21) +
  geom_point(aes(y=surgery_any_this_year, col = "Surgery", fill ="Surgery"),
             size=1.5,
             stroke =10,
             shape = 45) + 
  scale_y_reverse( limits = c(-1, 34), breaks = c(0,5,10,15,20,25,30, 35))+
  facet_wrap(~tumor_type)+
  theme_classic(base_size = 8)+
  theme(axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size =16),
        legend.title = element_blank(),
        legend.margin = unit(0, "mm")) +
  scale_color_manual(values = cols,
                     name="")+
  scale_fill_manual(values = fills,
                    name="")

p2 <- oligodendrogliomas %>%
  group_by(case_barcode)%>%
  mutate(max_follow_up_time = max(follow_up_time))%>%
  ungroup()%>%
  mutate(case_barcode = fct_reorder(case_barcode, max_follow_up_time, .desc = F)) %>%
  ggplot(aes(x = case_barcode, y = follow_up_time)) +
  geom_col(col = "#298c8c", fill= "#29b4b4",alpha=0.05, width = 0.6, linewidth =0.4)+
  labs(y = "Overall survival (years)", x = "")+
  geom_hline(yintercept = median.os.oligo, linetype = "twodash", color = "#29b4b4", size = 1) + 
  geom_segment(aes(y=alive_at_fu_this_year, yend=alive_at_fu_this_year+1.5, xend=case_barcode, linetype = ` `),
               arrow = arrow(length=unit(.3, 'cm')),
               col = "gray10",
               lwd=0.5)+
  geom_segment(aes(y=alive_at_fu_this_year, yend=alive_at_fu_this_year+1.5, xend=case_barcode), 
               arrow = arrow(length=unit(.3, 'cm')),
               col = "#298c8c",
               lwd=0.6)+
  geom_point(aes(y=radiotherapy_this_year, col = "Radiotherapy", fill ="Radiotherapy"),
             size=2,
             stroke =1,
             shape = 23) +
  geom_point(aes(y=tmz_concurrent_this_year, col = "Radiotherapy + Temozolomide (Concurrent)" , fill ="Radiotherapy + Temozolomide (Concurrent)" ),
             size=2,
             stroke =1,
             shape = 23) +
  geom_point(aes(y=tmz_this_year, col = "Temozolomide (Adjuvant)", fill ="Temozolomide (Adjuvant)"),
             size=2,
             stroke =1,
             shape = 24) +
  geom_point(aes(y=alkylating_other_this_year, col =  "Alkylating Chemotherapy (Non-Temozolomide)", fill = "Alkylating Chemotherapy (Non-Temozolomide)"),
             size=2,
             stroke =1,
             shape = 25) +
  geom_point(aes(y=platinum_based_this_year, col = "Platinum-based Chemotherapy", fill ="Platinum-based Chemotherapy"),
             size=1.5,
             stroke =1,
             shape = 22) +
  geom_point(aes(y=therapy_other_this_year, col = "Other therapy", fill ="Other therapy"),
             size=1.5,
             stroke =1,
             shape = 4) +
  geom_point(aes(y=age_under_25_this_year, col = "<25 yo", fill ="<25 yo"),
             size=2,
             stroke =1,
             shape = 21) +
  geom_point(aes(y=age_25_to_39_this_year, col = "25-39 yo", fill ="25-39 yo"),
             size=2,
             stroke =1,
             shape = 21) +
  geom_point(aes(y=age_40_to_55_this_year, col = "40-55 yo", fill ="40-55 yo"),
             size=2,
             stroke =1,
             shape = 21) +
  geom_point(aes(y=age_above_55_this_year, col = ">55 yo", fill =">55 yo"),
             size=2,
             stroke =1,
             shape = 21) +
  geom_point(aes(y=surgery_any_this_year, col = "Surgery", fill ="Surgery"),
             size=1.5,
             stroke =10,
             shape = 45) + 
  scale_y_reverse( limits = c(-1, 34), breaks = c(0,5,10,15,20,25,30, 35))+
  facet_wrap(~tumor_type)+
  theme_classic(base_size = 8)+
  theme(axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size =16),
        legend.title = element_blank(),
        legend.margin = unit(0, "mm")) +
  scale_color_manual(values = cols,
                     name="")+
  scale_fill_manual(values = fills,
                    name="") 

p2 + p1 + plot_layout(guides = "collect") &
  theme(legend.position='bottom')

survfit(Surv(oligo$follow_up_time, oligo$os) ~ tumor_type, data = oligo)
survfit(Surv(astro$follow_up_time, astro$os) ~ tumor_type, data = astro)