library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)


hugo_id <- "syn11536071"
expr_id <- "syn12667056"
anno_id <- "syn12667068"

upload_id  <- "syn12667035"


script_url <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE93722/create_processed_tables.R"
dataset <- "GSE93722"
activity_name <- "process GEO files into usable tables"


source("../../scripts/utils.R")
synLogin()

hugo_df <-  create_df_from_synapse_id(hugo_id)

anno_df <- anno_id %>% 
    create_df_from_synapse_id(unzip = T, skip = 31, nrow = 33) %>%
    rename("title" = `!Sample_title`) %>% 
    .[c(7,9,10),] %>% 
    mutate(title = c("tissue", "gender", "age")) %>% 
    transpose_df("title", "sample") %>% 
    mutate(gender = str_remove(gender, "gender: ")) %>% 
    mutate(age = str_remove(age, "age: "))


path <- download_from_synapse(expr_id)
system(str_c("cp ", path, " ."))
system(str_c("tar -xvf ", basename(path)))
samples <- 
    list.files() %>% 
    keep(str_detect(., "genes.results.txt.gz")) %>% 
    str_match("[:alnum:]+_(LAU[0-9]+).genes.results.txt.gz") %>% 
    .[,2]

tpm_df <- 
    list.files() %>% 
    keep(str_detect(., "genes.results.txt.gz")) %>%
    str_c("zcat ", .) %>% 
    map(fread) %>% 
    map(select, gene_id, TPM) %>% 
    reduce(full_join, by = "gene_id") %>% 
    set_colnames(c("ensembl_gene_id", samples)) %>% 
    left_join(hugo_df) %>% 
    select(hgnc_symbol, everything()) %>% 
    select(-ensembl_gene_id) %>% 
    set_colnames(c("Hugo", samples)) %>% 
    group_by(Hugo) %>% 
    summarise_all(sum) %>% 
    filter(!Hugo == "") %>% 
    ungroup %>% 
    .[,order(colnames(.))] %>% 
    select(Hugo, everything()) %>% 
    arrange(Hugo)

log_tpm_df <- tpm_df %>% 
    df_to_matrix("Hugo") %>% 
    add(1) %>% 
    log10 %>% 
    matrix_to_df("Hugo")

write_tsv(tpm_df, "expression_linear.tsv")
write_tsv(log_tpm_df, "expression_log.tsv")
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
    rnaseq_normalization = "TPM",
    expression_space = c("log10", "linear")
)

annotation_manifest_df <- tibble(
    path = "annotation.tsv",
    parent = upload_id,
    used = str_c(hugo_id, anno_id, expr_id, sep = ";"),
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    file_type = "annotations",
    annotations = str_c(colnames(anno_df)[-1], collapse = ";")
)

write_tsv(expression_manifest_df, "expression_manifest.tsv")
write_tsv(annotation_manifest_df, "annotation_manifest.tsv")

syncToSynapse("expression_manifest.tsv")
syncToSynapse("annotation_manifest.tsv")


