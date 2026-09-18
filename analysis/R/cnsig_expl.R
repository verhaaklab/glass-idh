### Load libraries
library(tidyverse)
library(odbc)
library(DBI)
library(stringr)
library(fs)
library(vroom)

### GLASS database connection
con <- DBI::dbConnect(odbc::odbc(), "GLASS5")

## Read in CNsig results
cnsig_files <- dir_ls(
 c("/vast/palmer/pi/verhaak/shared/glass_oligo_plus/results/sigprofiler/assignment/cnv/cosmic/",
  "/vast/palmer/pi/verhaakplus/glass/results/sigprofiler/assignment/cnv/cosmic/"
  ),
  recurse = TRUE,
  glob    = "*/*/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
)

cnsig_qc_query <- vroom::vroom(
  cnsig_files,
  delim = "\t",
  .name_repair = "minimal",
  progress = TRUE
)

ascat_files <- dir_ls(
 c(
  "/vast/palmer/pi/verhaak/shared/glass_oligo_plus/results/ascat/single/qc/",
  "/vast/palmer/pi/verhaakplus/glass/results/ascat/single/qc/"
  ),
  recurse = TRUE,
  glob    = "*.QC_metrics.txt"
)

ascat_qc_query <- vroom::vroom(
  ascat_files,
  delim = "\t",
  .name_repair = "minimal",
  progress = TRUE
)

#Export dataframes
write_tsv(cnsig_qc_query, "/vast/palmer/pi/verhaak/ek766/vscode_output/glass5_cnsig_20260206.tsv")
write_tsv(ascat_qc_query, "/vast/palmer/pi/verhaak/ek766/vscode_output/glass5_ascat_qc_20260206.tsv")


#Exploratory analyses
cnsig_qc_query %>% view()

cnsig_qc_query %>% gather(key, value, -Samples) %>% 
  ggplot(aes(x=key, y = value)) + geom_boxplot()

cnsig_qc_query %>% gather(key, value, -Samples) %>% 
  ggplot(aes(x=Samples, y = value)) + geom_boxplot() + facet_wrap(~key)
