library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE84697/"

expr_id   <- "syn12299847"
anno_id   <- "syn12299846"
upload_id <- "syn12299844"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()


anno_df <- anno_id %>% 
    create_df_from_synapse_id(unzip = T, skip = 20, nrows = 30) %>% 
    .[c(11,13), -1] %>% 
    inset("attribute", value = c("cell_type", "patient")) %>% 
    transpose_df("attribute", "sample") %>% 
    mutate(cell_type = str_remove_all(cell_type, "cell type: ")) %>% 
    mutate(patient = str_remove_all(patient, "patient identifier: ")) %>% 
    filter(cell_type != "CD45- non-immune cells") %>% 
    arrange(sample)

expr_df <- create_df_from_synapse_id(expr_id, unzip = T) %>% 
    select(-CLASS) %>% 
    rename("Hugo" = GENE) %>% 
    filter(!str_detect(Hugo, "^[0-9]{1,2}-[:alnum:]{3}")) %>% 
    select(c("Hugo", anno_df$sample))

activity_obj <- Activity(
    name = "create",
    description = "process GEO file into usable table",
    used = list(expr_id, anno_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE84697/create_processed_tables.R")
)

write_tsv(expr_df, "expression.tsv")
upload_file_to_synapse("expression.tsv", upload_id, activity_obj = activity_obj)

write_tsv(anno_df, "annotation.tsv")
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)

