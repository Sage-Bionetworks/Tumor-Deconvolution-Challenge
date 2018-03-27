library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(biomaRt)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE83115/"



expr_id  <- "syn11976809"


upload_id  <- "syn11976795"


setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()


ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

query_df <- getBM(attributes = 
                      c('ensembl_transcript_id',
                        'ensembl_gene_id', 
                        'hgnc_symbol'), 
                  mart = ensembl_mart)



expr_path <- download_from_synapse(expr_id)
system(str_c("gunzip ", expr_path))

expr_df <- expr_path %>% 
    str_sub(end = -4) %>% 
    read.table %>% 
    rownames_to_column("ensembl_transcript_id") %>% 
    as_data_frame %>% 
    left_join(query_df) %>% 
    dplyr::select(-ensembl_transcript_id) %>% 
    rename("Ensembl" = ensembl_gene_id) %>% 
    rename("Hugo" = hgnc_symbol) %>% 
    group_by(Ensembl, Hugo) %>% 
    .[complete.cases(.),] %>% 
    summarise_all(.funs = sum) %>% 
    .[,order(colnames(.))] %>% 
    dplyr::select(
        Ensembl,
        Hugo,
        B_cells_BGI, 
        B_cells_COH,
        CD4._T_cells_BGI,
        CD4._T_cells_COH,
        CD8._T_cells_BGI,
        CD8._T_cells_COH,
        Monocytes_BGI,
        Monocytes_COH)


activity_obj <- Activity(
    name = "create",
    description = "process GEO file into usable table",
    used = list(expr_id ),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE83115/create_processed_tables.R")
)

write_tsv(expr_df, "fpkm.tsv")

upload_file_to_synapse("fpkm.tsv", upload_id, activity_obj = activity_obj)

