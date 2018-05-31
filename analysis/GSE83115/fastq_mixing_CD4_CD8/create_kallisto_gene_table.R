library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(biomaRt)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE83115/"
upload_id <- "syn12231588"
download_id <- upload_id

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()


registerDoMC(cores = detectCores() -1)


ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

query_df <- getBM(attributes = 
                      c('ensembl_transcript_id',
                        'ensembl_gene_id', 
                        'hgnc_symbol'), 
                  mart = ensembl_mart)


file_df <- str_c('select id, name from file where parentId=="', download_id, '"') %>% 
    synQuery %>%
    use_series("results") %>% 
    map(data.frame) %>% 
    bind_rows %>% 
    as_data_frame 

tpm_df <- file_df %>% 
    use_series(file.id) %>% 
    llply(create_df_from_synapse_id, .parallel = T) %>% 
    map(dplyr::select, target_id, tpm) %>% 
    reduce(left_join, by = "target_id") %>% 
    set_colnames(c("transcript_id", file_df$file.name)) %>% 
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
    used = list(file_df$file.id),
    executed = "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE83115/fastq_mixing_CD4_CD8/create_kallisto_gene_table.R"
)
    
upload_file_to_synapse("kallisto_gene_tpms.tsv", upload_id, activity_obj = activity_obj)    
    

