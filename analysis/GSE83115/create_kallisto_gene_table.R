library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(biomaRt)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE83115/"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

upload_id <- "syn12180291"

registerDoMC(cores = detectCores() -1)


ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

query_df <- getBM(attributes = 
                      c('ensembl_transcript_id',
                        'ensembl_gene_id', 
                        'hgnc_symbol'), 
                  mart = ensembl_mart)



ids <- synQuery('select id from file where parentId=="syn12180291"') %>% 
    use_series(results) %>% 
    unlist %>% 
    unname 

paths <- map_chr(ids, download_from_synapse)

names <- paths %>% 
    basename %>% 
    str_remove("_kallisto_abundance.tsv")

paths %>% 
    llply(fread, .parallel = T) %>% 
    map(as_data_frame) %>% 
    map(dplyr::select, target_id, tpm) %>% 
    reduce(left_join, by = "target_id") %>% 
    set_colnames(c("transcript_id", names)) %>% 
    mutate(transcript_id = str_sub(transcript_id, end = -3)) %>% 
    left_join(query_df, by = c("transcript_id" = "ensembl_transcript_id")) %>% 
    dplyr::select(-c(transcript_id, ensembl_gene_id)) %>% 
    rename("Hugo" = hgnc_symbol) %>% 
    group_by(Hugo) %>%
    summarise_all(sum) %>% 
    filter(!Hugo == "") %>% 
    write_tsv("kallisto_gene_tpms.tsv")

activity_obj = Activity(
    name = "create and upload",
    description = "download kallisto quant files, merge, and collapse by hugo gene name",
    used = list(ids),
    executed = "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE83115/create_kallisto_gene_table.R"
)
    
upload_file_to_synapse("kallisto_gene_tpms.tsv", upload_id, activity_obj = activity_obj)    
    

