library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE62408/"

expr_id   <- "syn11915424"
upload_id <- "syn11915389"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

rpkm_df <- expr_id %>% 
    create_df_from_synapse_id(unzip = T) %>% 
    dplyr::rename("Hugo" = `-`) %>% 
    group_by(Hugo) %>% 
    summarise_all(.funs = sum) %>% 
    ungroup %>% 
    .[,order(colnames(.))] %>% 
    select(Hugo, everything()) %>% 
    write_tsv("rpkm.tsv")


activity_obj <- Activity(
    name = "create",
    description = "process GEO file into usable table",
    used = list(expr_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE62408/create_processed_tables.R")
)

upload_file_to_synapse("rpkm.tsv", upload_id, activity_obj = activity_obj)


