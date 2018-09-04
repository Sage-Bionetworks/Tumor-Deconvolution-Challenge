library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(GEOquery)
library(AnnotationDbi)
library(hgu133a.db)
library(preprocessCore)


home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE65134/"

upload_id    <- "syn15664979"
gt_upload_id <- "syn15664978"

expr_gse65133_id <- "syn15667753"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()


m_expr_gse65133_v <- expr_gse65133_id %>% 
    create_df_from_synapse_id %>% 
    df_to_matrix("Hugo") %>% 
    rowMeans() 

gse <- getGEO("GSE65134", GSEMatrix = TRUE)

series_df <- 
    pData(phenoData(gse[[1]])) %>% 
    rownames_to_column("sample") %>% 
    as_data_frame %>% 
    filter(`disease status:ch1` == "Healthy") %>% 
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
        hgu133a.db, 
        keys=keys(hgu133a.db,keytype="PROBEID"),
        columns=c("SYMBOL"), 
        keytype="PROBEID") %>% 
    as_data_frame() %>% 
    set_colnames(c("Affy", "Hugo")) %>% 
    drop_na()

expr_df <- 
    assayDataElement(gse$GSE65134_series_matrix.txt.gz@assayData, 'exprs') %>% 
    matrix_to_df("Affy") %>% 
    dplyr::select(c("Affy", anno_df$sample)) %>% 
    inner_join(query_df) %>% 
    arrange(Hugo) %>% 
    dplyr::select(Hugo, everything()) %>% 
    dplyr::select(-Affy) %>% 
    group_by(Hugo) %>% 
    summarise_all(max) 

log_expr_df <- expr_df %>% 
    gather(key = "sample", value = "expr", -Hugo) %>% 
    mutate(expr = log2(expr)) %>% 
    spread(key = "sample", value = "expr")

norm_log_expr_df <- expr_df %>% 
    df_to_matrix("Hugo") %>% 
    normalize.quantiles.use.target(m_expr_gse65133_v) %>% 
    set_rownames(expr_df$Hugo) %>% 
    matrix_to_df("Hugo") %>% 
    set_colnames(colnames(expr_df))


activity_obj <- Activity(
    name = "create",
    description = "process GEO data into usable tables",
    used = list(),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE65134/create_processed_tables.R")
)

activity_obj2 <- Activity(
    name = "create",
    description = "process GEO data into usable tables, normalize with GSE65133",
    used = list(expr_gse65133_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE65134/create_processed_tables.R")
)

write_tsv(log_expr_df, "log_expression_affy.tsv")
write_tsv(log_expr_df, "normalized_log_expression_affy.tsv")
write_tsv(expr_df, "expression_affy.tsv")
write_tsv(anno_df, "annotation.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")

upload_file_to_synapse("log_expression_affy.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("normalized_log_expression_affy.tsv", upload_id, activity_obj = activity_obj2)
upload_file_to_synapse("expression_affy.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("ground_truth.tsv", gt_upload_id, activity_obj = activity_obj)

