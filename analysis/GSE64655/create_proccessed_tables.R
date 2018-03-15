library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE64655/"


hugo_id          <- "syn11536071"
GSE64655_expr_id <- "syn11969378"
GSE64655_anno_id <- "syn11969387"

upload_id  <- "syn11969377"


setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

hugo_df <-  create_df_from_synapse_id(hugo_id)

GSE64655_expr_df <- GSE64655_expr_id %>%
    create_df_from_synapse_id(unzip = T, skip = 3) %>% 
    select(-c(`Gene Type`, Description, `Gene Symbol`)) %>% 
    left_join(hugo_df, by = c("Gene ID" = "ensembl_gene_id")) %>% 
    rename("Ensembl" = `Gene ID`) %>% 
    rename("Hugo" = hgnc_symbol) %>% 
    .[,order(colnames(.))] %>% 
    select(Ensembl, Hugo, everything())

GSE64655_anno_df <- GSE64655_anno_id %>% 
    create_df_from_synapse_id(unzip = T, skip = 30, nrow = 7) %>%
    rename("title" = `!Sample_title`) %>% 
    filter(title == "!Sample_source_name_ch1") %>% 
    transpose_df("title", "sample") %>% 
    set_colnames(c("sample", "cell_type")) %>% 
    mutate(patient = str_sub(sample, end = 4)) %>% 
    mutate(days = str_sub(sample, start = -9, end = -9)) %>% 
    mutate(ABV1 = str_sub(sample, start = 6, end = -11)) %>% 
    mutate(ABV2 = ifelse(ABV1 == "DC", "mDC", 
                         ifelse(ABV1 == "MO", "Mono", 
                                ifelse(ABV1 == "NEU", "Neut", 
                                       ifelse(ABV1 == "PMBC", "PBMC", ABV1))))) %>% 
    mutate(sample = str_c(patient, "_", ABV2, "_d", days)) %>% 
    select(-c(ABV1, ABV2)) %>% 
    arrange(sample)

activity_obj <- Activity(
    name = "create",
    description = "process GEO files into usable tables",
    used = list(hugo_id, GSE64655_anno_id, GSE64655_expr_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE64655/create_processed_tables.R")
)

write_tsv(GSE64655_expr_df, "expression.tsv")
write_tsv(GSE64655_anno_df, "annotation.tsv")

upload_file_to_synapse("expression.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)

