library(tidyverse)
library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)
library(GEOquery)
library(AnnotationDbi)
library(illuminaHumanv4.db)


upload_id    <- "syn15664932"
gt_upload_id <- "syn15664931"

dataset       <- "GSE65133"
script_url    <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE65133/create_processed_tables.R"
activity_name <- "process_files_from_GEO"

source("../../scripts/utils.R")
synLogin()

gse <- getGEO(dataset, GSEMatrix = TRUE)

series_df <- 
    pData(phenoData(gse[[1]])) %>% 
    rownames_to_column("sample") %>% 
    as_tibble %>% 
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
    as_tibble() %>% 
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

linear_expr_df <- log_expr_df %>% 
    gather(key = "sample", value = "expr", -Hugo) %>% 
    mutate(expr = 2^expr) %>% 
    spread(key = "sample", value = "expr")


write_tsv(log_expr_df, "expression_log.tsv")
write_tsv(linear_expr_df, "expression_linear.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")
write_tsv(anno_df, "annotation.tsv")


expression_manifest_df <- tibble(
    path = c("expression_log.tsv", "expression_linear.tsv"),
    parent = upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    file_type = "expression",
    expression_type = "microarray", 
    microarray_type = "Illumina Human HT-12 V4 BeadChip",
    expression_space = c("log2", "linear")
)

ground_truth_manifest_df <- tibble(
    path = "ground_truth.tsv",
    parent = gt_upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    file_type = "ground truth",
    unit = "fraction",
    cell_types = str_c(colnames(ground_truth_df)[-1], collapse = ";")
)


write_tsv(expression_manifest_df, "expression_manifest.tsv")
write_tsv(ground_truth_manifest_df, "ground_truth_manifest.tsv")

syncToSynapse("expression_manifest.tsv")
syncToSynapse("ground_truth_manifest.tsv")



    
activity_obj <- Activity(
    name = "create",
    description = "process GEO data into usable tables",
    used = list(),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE65133/create_processed_tables.R")
)

 

upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)

