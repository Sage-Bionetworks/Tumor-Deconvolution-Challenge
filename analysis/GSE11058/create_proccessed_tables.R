library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(GEOquery)
library(AnnotationDbi)
library(hgu133plus2.db)

upload_id     <- "syn17038053"
gt_upload_id  <- "syn17038052"

source("../../scripts/utils.R")
synLogin()



# from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0006098
mix_df <- data_frame(
    "Mix" = c("MixA", "MixB", "MixC", "MixD"),
    "TCell_Jurkat" = c(25, 5, 1, 0.2),
    "BCell_IM9" = c(12.5, 31.7, 49.5, 33.27),
    "BCell_Raji" = c(25, 47.5, 16.5, 33.27),
    "Monocyte_THP1" = c(37.5, 15.8, 33.0, 33.27)
)

gse_object <- getGEO("GSE11058", GSEMatrix = TRUE)

ground_truth_df <- gse_object %>% 
    extract2(1) %>% 
    phenoData() %>% 
    pData() %>%
    rownames_to_column("sample") %>%
    as_data_frame %>% 
    filter(characteristics_ch1.1 == "mixed cells") %>% 
    dplyr::select(sample, characteristics_ch1) %>% 
    inner_join(mix_df, by = c("characteristics_ch1" = "Mix")) %>% 
    dplyr::select(-characteristics_ch1)


query_df <-
    AnnotationDbi::select(
        hgu133plus2.db,
        keys=keys(hgu133plus2.db,keytype="PROBEID"),
        columns=c("SYMBOL"),
        keytype="PROBEID") %>%
    as_data_frame() %>%
    set_colnames(c("Illum", "Hugo")) %>%
    drop_na()

expr_df <- gse_object$GSE11058_series_matrix.txt.gz@assayData %>% 
    assayDataElement('exprs') %>%
    matrix_to_df("Illum") %>%
    inner_join(query_df) %>%
    dplyr::select(Hugo, everything()) %>%
    dplyr::select(-Illum) %>%
    group_by(Hugo) %>%
    summarise_all(mean) %>%
    drop_na() %>%
    filter(Hugo != "") %>% 
    gather(key = "sample", value = "expr", -Hugo)

samples_in_common <- intersect(expr_df$sample, ground_truth_df$sample)

expr_df <- expr_df %>% 
    filter(sample %in% samples_in_common) %>% 
    spread(key = "sample", value = "expr")

ground_truth_df <- ground_truth_df %>% 
    filter(sample %in% samples_in_common)
  
activity_obj <- Activity(
    name = "create",
    description = "process GEO data into usable tables",
    used = list("https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0006098"),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE11058/create_processed_tables.R")
)

write_tsv(expr_df, "expression_illumina.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")

upload_file_to_synapse("expression_illumina.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("ground_truth.tsv", gt_upload_id, activity_obj = activity_obj)

