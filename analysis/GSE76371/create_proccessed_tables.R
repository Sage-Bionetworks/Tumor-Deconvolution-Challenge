library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE76371/"

expr_id    <- "syn12036837"
anno_id    <- "syn12036846"
upload_id  <- "syn12036836"


setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

anno_df <- anno_id %>% 
    download_from_synapse %>% 
    str_c("zcat ", .) %>% 
    fread(skip = 33, nrows = 39) %>% 
    rename("values" = "!Sample_geo_accession") %>% 
    .[15,] %>% 
    transpose_df("values", "sample") %>% 
    set_colnames(c("sample", "cell_type")) %>% 
    inset("cdna_method", value = c(rep("NuGEN", 13), rep("Clontech", 4)))

expr_id %>% 
    download_from_synapse %>%
    str_c("tar -xvf ", .) %>% 
    system(intern = T) %>% 
    str_c("gunzip ", .) %>% 
    walk(system)

files <- list.files() %>% 
    keep(str_detect(., "RPKMforgenes.txt$"))

sample_names <- files %>% 
    str_split("_") %>% 
    map_chr(extract2, 1)

count_dfs <- files %>% 
    map(fread) %>% 
    map(as_data_frame) %>% 
    map(select, c(1,4)) %>% 
    map(set_colnames, c("Hugo", "Count")) 

genes <- count_dfs[[1]]$Hugo

count_df <- count_dfs %>% 
    map(select, 2) %>% 
    bind_cols %>% 
    set_colnames(sample_names) %>% 
    inset("Hugo", value = genes) %>% 
    select(Hugo, everything())



activity_obj <- Activity(
    name = "create",
    description = "process GEO file into usable table",
    used = list(expr_id, anno_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE763715/create_processed_tables.R")
)

write_tsv(count_df, "counts.tsv")
upload_file_to_synapse("counts.tsv", upload_id, activity_obj = activity_obj)

write_tsv(anno_df, "annotation.tsv")
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)

