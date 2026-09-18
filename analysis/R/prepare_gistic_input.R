##################################################
# Prepare input files for running GISTIC seperately for primaries and recurrences
# Ignores multi-sector samples (one sample per patient and timepoint)
# Updated: 2024.12.11
# Author: Floris B
# Editor: Tamrin C
##################################################

library(DBI)
library(tidyverse)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "glass5")

q <- read_file("gistic_prepare.sql")

qres <- dbGetQuery(con, q)
seg = qres %>% filter(complete.cases(start,end,num_snps))

seg_p = seg %>% filter(sample_type == "P") %>% select(-sample_type)
seg_r = seg %>% filter(sample_type == "R") %>% select(-sample_type)

write.table(seg_p, file = "input/primary.seg", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(seg_r, file = "input/recurrence.seg", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

markers = data.frame(id = 1:(2*nrow(seg)), chr = c(seg$chrom, seg$chrom), pos = c(seg$start, seg$end), stringsAsFactors = FALSE) %>% distinct()
write.table(markers, file = "input/markers.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## Input separated by subtypes

clinical_subtypes <- dbGetQuery(con, "SELECT * FROM clinical.subtypes")
seg$case_barcode <- sapply(strsplit(seg$aliquot_barcode, "-"), function(x) paste(x[1:3], collapse = "-"))
seg <- merge(seg, clinical_subtypes, by = "case_barcode", all.x = TRUE)

seg_wt <- seg %>% filter(idh_codel_subtype=="IDHwt")
seg_codel <- seg %>% filter(idh_codel_subtype=="IDHmut-codel")
seg_noncodel <- seg %>% filter(idh_codel_subtype=="IDHmut-noncodel")

seg_wt_p = seg_wt %>% filter(sample_type == "P") %>% select(-sample_type, -case_barcode, -idh_codel_subtype)
seg_wt_r = seg_wt %>% filter(sample_type == "R") %>% select(-sample_type, -case_barcode, -idh_codel_subtype)
write.table(seg_wt_p, file = "input/idhwt_primary.seg", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(seg_wt_r, file = "input/idhwt_recurrence.seg", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

seg_codel_p = seg_codel %>% filter(sample_type == "P") %>% select(-sample_type, -case_barcode, -idh_codel_subtype)
seg_codel_r = seg_codel %>% filter(sample_type == "R") %>% select(-sample_type, -case_barcode, -idh_codel_subtype)
write.table(seg_codel_p, file = "input/codel_primary.seg", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(seg_codel_r, file = "input/codel_recurrence.seg", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

seg_noncodel_p = seg_noncodel %>% filter(sample_type == "P") %>% select(-sample_type, -case_barcode, -idh_codel_subtype)
seg_noncodel_r = seg_noncodel %>% filter(sample_type == "R") %>% select(-sample_type, -case_barcode, -idh_codel_subtype)
write.table(seg_noncodel_p, file = "input/noncodel_primary.seg", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(seg_noncodel_r, file = "input/noncodel_recurrence.seg", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
