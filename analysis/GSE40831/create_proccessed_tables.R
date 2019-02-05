library(tidyverse)
library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)
library(GEOquery)
library(AnnotationDbi)
library(hgu133plus2.db)
library(xlsx)

upload_id     <- "syn18139812"
gt_id         <- "syn18139845"

source("../../scripts/utils.R")
synLogin()

script_url    <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE40831/create_processed_tables.R"
dataset       <- "GSE40831"
activity_name <- "process GEO files into usable tables"
used          <- gt_id



gse_object <- getGEO(dataset, GSEMatrix = TRUE)

ground_truth_df <- gt_id %>% 
    download_from_synapse() %>% 
    read.xlsx(sheetIndex = 1, startRow = 2, stringsAsFactors = F) %>% 
    as_tibble() %>% 
    dplyr::select(Reference.populations, Flow.cytometry) %>% 
    set_colnames(c("cell_type", "percent")) %>% 
    mutate(percent = as.double(percent)) %>% 
    filter(!is.na(percent))
    



df <- gse_object %>% 
    extract2(2) %>% 
    phenoData() %>% 
    pData() %>%
    rownames_to_column("sample") %>%
    as_tibble()




# %>% 
#     filter(characteristics_ch1.1 == "mixed cells") %>% 
#     dplyr::select(sample, characteristics_ch1)
# %>% 
#     inner_join(mix_df, by = c("characteristics_ch1" = "Mix")) %>% 
#     dplyr::select(-characteristics_ch1)


query_df <-
    AnnotationDbi::select(
        hgu133plus2.db,
        keys=keys(hgu133plus2.db,keytype="PROBEID"),
        columns=c("SYMBOL"),
        keytype="PROBEID") %>%
    as_tibble() %>%
    set_colnames(c("AFFY", "Hugo")) %>%
    drop_na()

expr_df <- gse_object$`GSE40831-GPL571_series_matrix.txt.gz`@assayData %>%
    assayDataElement('exprs') %>%
    matrix_to_df("AFFY") %>%
    inner_join(query_df) %>%
    dplyr::select(Hugo, everything()) %>%
    dplyr::select(-AFFY) %>%
    group_by(Hugo) %>%
    summarise_all(mean) %>%
    drop_na() %>%
    filter(Hugo != "") %>%
    gather(key = "sample", value = "expr", -Hugo)
# 
# samples_in_common <- intersect(expr_df$sample, ground_truth_df$sample)
# 
# expr_df <- expr_df %>% 
#     filter(sample %in% samples_in_common) %>% 
#     spread(key = "sample", value = "expr")
# 
# ground_truth_df <- ground_truth_df %>% 
#     filter(sample %in% samples_in_common)
#   
# activity_obj <- Activity(
#     name = "create",
#     description = "process GEO data into usable tables",
#     used = list("https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0006098"),
#     executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE11058/create_processed_tables.R")
# )
# 
# write_tsv(expr_df, "expression_illumina.tsv")
# write_tsv(ground_truth_df, "ground_truth.tsv")
# 
# upload_file_to_synapse("expression_illumina.tsv", upload_id, activity_obj = activity_obj)
# upload_file_to_synapse("ground_truth.tsv", gt_upload_id, activity_obj = activity_obj)
# 
