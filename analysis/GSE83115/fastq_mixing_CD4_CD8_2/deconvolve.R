library(synapseClient)
library(tidyverse)
library(data.table)
library(magrittr)


tmp_dir   <- "/home/aelamb/tmp/tumor_deconvolution/GSE83115/"
tpm_id    <- "syn12385635"
upload_id <- "syn12369867"

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


setwd(tmp_dir)
synapseClient::synapseLogin()



tpm_df <- tpm_id %>%
    create_df_from_synapse_id()

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
    used = list(tpm_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE83115/fastq_mixing_CD4_CD8_2/deconvolve.R")
)

upload_file_to_synapse("cibersort_results.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("mcpcounter_results.tsv", upload_id, activity_obj = activity_obj)
