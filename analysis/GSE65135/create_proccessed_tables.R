library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(GEOquery)
library(biomaRt)


home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE65135/"

upload_id    <- "syn15664986"
gt_upload_id <- "syn15664985"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

gse <- getGEO("GSE65135", GSEMatrix = TRUE)

series_df <- 
    pData(phenoData(gse[[1]])) %>% 
    rownames_to_column("sample") %>% 
    as_data_frame %>% 
    dplyr::select(sample, title, `flow cytometry cell subset proportions:ch1`) %>% 
    set_colnames(c("sample", "id", "cell_types")) %>% 
    filter(cell_types != "NA") %>% 
    separate(cell_types, sep = "; ", into = as.character(1:20), fill = "right") %>%
    gather(key = "key", value = "value", -c(sample, id)) %>% 
    dplyr::select(-key) %>% 
    drop_na %>% 
    separate(value, sep = " = ", into = c("cell_type", "percent")) %>% 
    mutate(cell_type = str_replace_all(cell_type, " ", "_")) %>% 
    mutate(cell_type = str_replace_all(cell_type, "ï", "i")) %>% 
    mutate(percent = str_remove_all(percent, "%")) %>% 
    spread(key = "cell_type", value = "percent")


anno_df <-         dplyr::select(series_df, sample, id)
ground_truth_df <- dplyr::select(series_df, -id)


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

query_df <- 
    getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol"), mart = ensembl) %>% 
    drop_na() %>% 
    set_colnames(c("Affy", "Hugo")) %>% 
    filter(Affy != "") %>% 
    filter(Hugo != "")

expr_df <- 
    assayDataElement(gse$GSE65135_series_matrix.txt.gz@assayData, 'exprs') %>% 
    matrix_to_df("Affy") %>% 
    dplyr::select(c("Affy", anno_df$sample)) %>% 
    inner_join(query_df) %>% 
    dplyr::select(Hugo, everything()) %>% 
    dplyr::select(-Affy) %>% 
    group_by(Hugo) %>% 
    summarise_all(max) 


activity_obj <- Activity(
    name = "create",
    description = "process GEO data into usable tables",
    used = list(),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE65135/create_processed_tables.R")
)

write_tsv(expr_df, "expression_affy.tsv")
write_tsv(anno_df, "annotation.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")

upload_file_to_synapse("expression_affy.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("ground_truth.tsv", gt_upload_id, activity_obj = activity_obj)

