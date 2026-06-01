# GLASS-I: Prior Treatment Table
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

# Let's start fresh
rm(list = ls())
gc()

# load libraries
library(RPostgres)
library(tidyverse)

# Functions to change TMZ, RT, ALK_CHEMO and TMZ_cycles columns for prior treatment
modify_tmz <- function(df) {
  df %>%
    arrange(case_barcode, surgery_number) %>%
    group_by(case_barcode) %>%
    mutate(
      TMZ = {
        if (all(is.na(tmz))) {
          rep(NA_real_, n())
        } else {
          x <- cummax(if_else(is.na(tmz), -Inf, tmz))
          x <- lag(x, default = NA_real_)
          if_else(x == -Inf, NA_real_, x)
        }
      }
    ) %>%
    ungroup()
}

modify_rt <- function(df) {
  df %>%
    arrange(case_barcode, surgery_number) %>%
    group_by(case_barcode) %>%
    mutate(
      RT = {
        if (all(is.na(radiotherapy))) {
          rep(NA_real_, n())
        } else {
          x <- cummax(if_else(is.na(radiotherapy), -Inf, radiotherapy))
          x <- lag(x, default = NA_real_)
          if_else(x == -Inf, NA_real_, x)
        }
      }
    ) %>%
    ungroup()
}

modify_alk_chemo <- function(df) {
  df %>%
    arrange(case_barcode, surgery_number) %>%
    group_by(case_barcode) %>%
    mutate(
      ALK_CHEMO = {
        if (all(is.na(alk_chemo))) {
          rep(NA_real_, n())
        } else {
          x <- cummax(if_else(is.na(alk_chemo), -Inf, alk_chemo))
          x <- lag(x, default = NA_real_)
          if_else(x == -Inf, NA_real_, x)
        }
      }
    ) %>%
    ungroup()
}

modify_tmz_cycles <- function(df) {
  df %>%
    arrange(case_barcode, surgery_number) %>%
    group_by(case_barcode) %>%
    mutate(
      TMZ_cycles = {
        if (all(is.na(treatment_tmz_cycles))) {
          rep(NA_real_, n())
        } else {
          x <- cummax(if_else(is.na(treatment_tmz_cycles), -Inf, treatment_tmz_cycles))
          x <- lag(x, default = NA_real_)
          if_else(x == -Inf, NA_real_, x)
        }
      }
    ) %>%
    ungroup()
}

# load data
con <- dbConnect(RPostgres::Postgres(), service = "glass5reader")             ## data available on Synapse (Synapse ID: syn17038081)
  
surgeries <- dbGetQuery(con, "SELECT * FROM clinical.surgeries") %>%
    mutate(timing.surgery = substr(sample_barcode, 14,15),
           tmz = case_when(treatment_tmz == T  &  (treatment_concurrent_tmz == T|treatment_concurrent_tmz == F|is.na(treatment_concurrent_tmz))~ 1,
                           treatment_concurrent_tmz == T  &  (treatment_tmz == T |treatment_tmz == F|is.na(treatment_tmz))~ 1,
                           treatment_tmz == F  &  (treatment_concurrent_tmz == F |  is.na(treatment_concurrent_tmz))~ 0,
                           is.na(treatment_tmz)  &  treatment_concurrent_tmz == F~ 0),
           radiotherapy = if_else(treatment_radiotherapy == T, 1, 0),
           alk_chemo=       case_when(treatment_alkylating_agent == T  &  (treatment_concurrent_tmz == T|treatment_concurrent_tmz == F|is.na(treatment_concurrent_tmz))~ 1,
                                      treatment_concurrent_tmz == T  &  (treatment_alkylating_agent == T |treatment_alkylating_agent == F|is.na(treatment_alkylating_agent))~ 1,
                                      treatment_alkylating_agent == F  &  (treatment_concurrent_tmz == F |  is.na(treatment_concurrent_tmz))~ 0,
                                      is.na(treatment_alkylating_agent)  &  treatment_concurrent_tmz == F~ 0)) %>%
    filter(idh_codel_subtype != "IDHwt")

gold.set.samples<- dbGetQuery(con, "SELECT * FROM analysis.gold_set")%>%
    semi_join(surgeries, by = "case_barcode")%>%
    pivot_longer(
      cols = c(tumor_barcode_a, tumor_barcode_b),
      names_to = NULL,
      values_to = "aliquot_barcode"
    )%>%
    mutate(sample_barcode = substr(aliquot_barcode, 1, 15))%>%
  dplyr::select(sample_barcode)
  
dbDisconnect(con)

surgeries <- modify_alk_chemo(surgeries)
surgeries <- modify_tmz(surgeries)
surgeries <- modify_rt(surgeries)
surgeries <- modify_tmz_cycles(surgeries)

treatment.data <- gold.set.samples%>%
  left_join(surgeries, by = "sample_barcode")%>%
  dplyr::select(case_barcode, sample_barcode, surgery_number, TMZ, RT, ALK_CHEMO, TMZ_cycles)

write_csv(treatment.data, "~/prior_treatment_gold_set.csv")