library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE93722/"


hugo_id <- "syn11536071"
expr_id <- "syn12667056"
anno_id <- "syn12667068"

upload_id  <- "syn12667035"


setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

hugo_df <-  create_df_from_synapse_id(hugo_id)

anno_df <- anno_id %>% 
    create_df_from_synapse_id(unzip = T, skip = 17, nrow = 33) %>%
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
    set_colnames(c("Hugo", "Ensembl", samples))


activity_obj <- Activity(
    name = "create",
    description = "process GEO files into usable tables",
    used = list(hugo_id, anno_id, expr_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE93722/create_processed_tables.R")
)

write_tsv(tpm_df, "expression.tsv")
write_tsv(anno_df, "annotation.tsv")

upload_file_to_synapse("expression.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)

