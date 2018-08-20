library(doMC)
library(plyr)
library(synapser)
library(tidyverse)
library(data.table)
library(magrittr)
library(biomaRt)


home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE64655/"

manifest_id <- "syn13865548"
file_dir_id <- "syn13841771"

n_cores <- detectCores() -1
registerDoMC(n_cores)
synLogin()
setwd(tmp_dir)

source(str_c(home_dir, "scripts/utils.R"))

ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

query_df <- getBM(attributes = 
                      c('ensembl_transcript_id',
                        'ensembl_gene_id', 
                        'hgnc_symbol'), 
                  mart = ensembl_mart)

manifest_df <- manifest_id %>% 
    create_df_from_synapse_id 

file_df <- file_dir_id %>% 
    get_file_df_from_synapse_dir_id %>% 
    inner_join(manifest_df, by = c("file.name" = "output_name"))

count_df <- file_df %>% 
    use_series(file.id) %>% 
    llply(create_df_from_synapse_id, .parallel = T) %>% 
    map(dplyr::select, target_id, est_counts) %>% 
    reduce(left_join, by = c("target_id")) %>% 
    set_colnames(c(
        "transcript_id", 
        str_c(file_df$patient, file_df$cell_type, sep = "_"))) %>% 
    mutate(transcript_id = str_sub(transcript_id, end = -3)) %>% 
    left_join(query_df, by = c("transcript_id" = "ensembl_transcript_id")) %>% 
    dplyr::select(-c(transcript_id, ensembl_gene_id)) %>% 
    rename("Hugo" = hgnc_symbol) %>% 
    group_by(Hugo) %>%
    summarise_all(sum) %>% 
    filter(!Hugo == "") 

write_tsv(count_df, "kallisto_gene_counts.tsv")



activity_obj = Activity(
    name = "create and upload",
    description = "download kallisto quant files, merge, and collapse by hugo gene name",
    used = list(file_df$file.id),
    executed = "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE64655/fastq_mixing_CD4_CD8_2/create_kallisto_gene_table.R"
)

upload_file_to_synapse("kallisto_gene_counts.tsv", file_dir_id , activity_obj = activity_obj)   

