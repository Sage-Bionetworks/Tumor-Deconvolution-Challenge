library(tidyverse)
library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)

hugo_id <- "syn11536071"
expr_id <- "syn11969378"
anno_id <- "syn11969387"

upload_id  <- "syn12667653"

script_url <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE64655/create_processed_tables.R"
dataset <- "GSE64655"
activity_name <- "process GEO files into usable tables"


source("../../scripts/utils.R")
synLogin()


anno_df <- anno_id %>% 
    create_df_from_synapse_id(unzip = T, skip = 30, nrow = 7) %>%
    dplyr::rename("title" = `!Sample_title`) %>% 
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
    arrange(sample) %>% 
    filter(days == 0) %>% 
    select(-days)



hugo_df <-  create_df_from_synapse_id(hugo_id)


tmm_df <- expr_id %>%
    create_df_from_synapse_id(unzip = T, skip = 3) %>% 
    select(-c(`Gene Type`, Description, `Gene Symbol`)) %>% 
    left_join(hugo_df, by = c("Gene ID" = "ensembl_gene_id")) %>% 
    dplyr::rename("Hugo" = hgnc_symbol) %>% 
    select(one_of(c("Hugo", anno_df$sample))) %>% 
    group_by(Hugo) %>% 
    summarise_all(sum) %>% 
    filter(!Hugo == "") %>% 
    ungroup 

log_tmm_df <- tmm_df %>% 
    df_to_matrix("Hugo") %>% 
    add(1) %>% 
    log10 %>% 
    matrix_to_df("Hugo")
    
write_tsv(tmm_df, "expression_linear.tsv")
write_tsv(log_tmm_df, "expression_log.tsv")
write_tsv(anno_df, "annotation.tsv")


expression_manifest_df <- tibble(
    path = c("expression_log.tsv", "expression_linear.tsv"),
    parent = upload_id,
    used = str_c(hugo_id, anno_id, expr_id, sep = ";"),
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    file_type = "expression",
    expression_type = "RNASeq", 
    rnaseq_type = "paired-end",
    rnaseq_normalization = "TMM",
    expression_space = c("log10", "linear")
)

write_tsv(expression_manifest_df, "expression_manifest.tsv")

syncToSynapse("expression_manifest.tsv")




activity_obj <- Activity(
    name = "create",
    description = "process GEO files into usable tables",
    used = list(hugo_id, anno_id, expr_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE64655/create_processed_tables.R")
)


upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)

