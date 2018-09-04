library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(GEOquery)
library(AnnotationDbi)
library(illuminaHumanv4.db)


home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE65133/"

upload_id    <- "syn15664932"
gt_upload_id <- "syn15664931"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

gse <- getGEO("GSE65133", GSEMatrix = TRUE)

series_df <- 
    pData(phenoData(gse[[1]])) %>% 
    rownames_to_column("sample") %>% 
    as_data_frame %>% 
    dplyr::select(sample, title, `flow cytometry cell subset proportions:ch1`) %>% 
    set_colnames(c("sample", "id", "cell_types")) %>% 
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

query_df <- 
    AnnotationDbi::select(
        illuminaHumanv4.db, 
        keys=keys(illuminaHumanv4.db,keytype="PROBEID"),
        columns=c("SYMBOL"), 
        keytype="PROBEID") %>% 
    as_data_frame() %>% 
    set_colnames(c("Illum", "Hugo")) %>% 
    drop_na()


log_expr_df <- 
    assayDataElement(gse$GSE65133_series_matrix.txt.gz@assayData, 'exprs') %>% 
    matrix_to_df("Illum") %>% 
    inner_join(query_df) %>% 
    dplyr::select(Hugo, everything()) %>% 
    dplyr::select(-Illum) %>% 
    group_by(Hugo) %>% 
    summarise_all(max) %>% 
    drop_na() %>% 
    filter(Hugo != "")

expr_df <- log_expr_df %>% 
    gather(key = "sample", value = "expr", -Hugo) %>% 
    mutate(expr = 2^expr) %>% 
    spread(key = "sample", value = "expr")
    
activity_obj <- Activity(
    name = "create",
    description = "process GEO data into usable tables",
    used = list(),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE65133/create_processed_tables.R")
)

write_tsv(log_expr_df, "log_expression_illumina.tsv")
write_tsv(expr_df, "expression_illumina.tsv")
write_tsv(anno_df, "annotation.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")

upload_file_to_synapse("log_expression_illumina.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("expression_illumina.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("ground_truth.tsv", gt_upload_id, activity_obj = activity_obj)

