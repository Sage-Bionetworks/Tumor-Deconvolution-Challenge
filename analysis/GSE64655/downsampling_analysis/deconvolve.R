library(plyr)
library(doMC)
library(synapseClient)
library(tidyverse)
library(data.table)
library(magrittr)
library(biomaRt)


tmp_dir     <- "/home/aelamb/tmp/tumor_deconvolution/GSE64655/"
download_id <- "syn12678224"
upload_id   <- "syn12678224"

create_df_from_synapse_id <- function(syn_id, location = NULL, unzip = F, ...){
    path <- download_from_synapse(syn_id, location)
    if(unzip) path <- stringr::str_c("zcat ", path)
    path %>% 
        data.table::fread(...) %>% 
        dplyr::as_data_frame() 
}

download_from_synapse <- function(syn_id, location = NULL){
    path <- synapseClient::synGet(syn_id, downloadLocation = location)@filePath
    return(path)
}

df_to_matrix <- function(df, id_column){
    df %>% 
        data.frame() %>% 
        tibble::column_to_rownames(id_column) %>% 
        as.matrix()
}


upload_file_to_synapse <- function(
    path, synapse_id, 
    annotation_list = NULL, 
    activity_obj = NULL, 
    ret = "entity"){
    
    entity <- synapseClient::File(
        path = path, 
        parentId = synapse_id, 
        annotations = annotation_list)
    entity <- synapseClient::synStore(entity, activity = activity_obj)
    if(ret == "entity") return(entity)
    if(ret == "syn_id") return(entity$properties$id)
}

get_file_df_from_synapse_dir_id <- function(syn_id){
    str_c('select id, name from file where parentId=="', download_id, '"') %>% 
        synQuery() %>%
        as_data_frame() 
}


setwd(tmp_dir)
synapseClient::synapseLogin()
registerDoMC(cores = detectCores() -1)

ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

query_df <- getBM(attributes = 
                      c('ensembl_transcript_id',
                        'ensembl_gene_id', 
                        'hgnc_symbol'), 
                  mart = ensembl_mart)



file_df <- download_id %>% 
    get_file_df_from_synapse_dir_id %>% 
    filter(!str_detect(file.name, "results"))

tpm_df <- file_df %>% 
    use_series(file.id) %>% 
    llply(create_df_from_synapse_id, .parallel = T) %>% 
    llply(dplyr::select, target_id, tpm, .parallel = T) %>% 
    reduce(left_join, by = c("target_id")) %>% 
    set_colnames(c("target_id", file_df$file.name)) %>% 
    mutate(ensembl_transcript_id = str_sub(target_id, end = -3)) %>% 
    dplyr::select(ensembl_transcript_id, everything()) %>% 
    left_join(query_df) %>% 
    dplyr::select(-c(target_id, ensembl_gene_id, ensembl_transcript_id)) %>% 
    filter(!hgnc_symbol == "") %>%  
    group_by(hgnc_symbol) %>% 
    summarise_all(sum) %>% 
    ungroup %>% 
    dplyr::rename(Hugo = hgnc_symbol)



write_tsv(tpm_df, "expr.tsv")

tpm_df %>% 
    df_to_matrix("Hugo") %>% 
    write.table("expr_matrix.tsv", sep = "\t", quote = F)

system(str_c(
    "cwltool /home/aelamb/repos/irwg/iatlas-tool-cibersort/Dockstore.cwl", 
    "--mixture_file expr.tsv", 
    "--sig_matrix_file /home/aelamb/repos/irwg/iatlas-tool-cibersort/sample.references.matrix.txt", 
    "--QN",
    "--output_file_string cibersort_results.tsv",
    sep = " "))

system(str_c(
    "cwltool /home/aelamb/repos/irwg/iatlas-tool-mcpcounter/Dockstore.cwl",
    "--input_expression_file expr_matrix.tsv",
    "--output_file_string mcpcounter_results.tsv",
    "--features_type HUGO_symbols",
    sep = " "))

activity_obj <- Activity(
    name = "create",
    description = "create and upload deconvolution results using cibersort and mcpcounter cwl files",
    used = list(download_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE64655/downsampling_analysis/deconvolve.R")
)

upload_file_to_synapse("cibersort_results.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("mcpcounter_results.tsv", upload_id, activity_obj = activity_obj)
