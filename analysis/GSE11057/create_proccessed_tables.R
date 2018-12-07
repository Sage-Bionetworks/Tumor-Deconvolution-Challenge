library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(GEOquery)
library(AnnotationDbi)
library(hgu133plus2.db)

upload_id     <- "syn17091411"
gt_upload_id  <- "syn17091410"

source("../../scripts/utils.R")
synLogin()


gse_object <- getGEO("GSE11057", GSEMatrix = TRUE)

geo_df <- gse_object %>% 
    extract2(1) %>% 
    phenoData() %>% 
    pData() %>%
    rownames_to_column("sample") %>%
    as_data_frame %>% 
    dplyr::select(sample, title, `fractional abundance:ch1`) %>% 
    set_colnames(c("sample", "title", "fraction")) %>% 
    separate(col = title, into = c("cell_type", "donor"), sep = ", donor ")

translation_df <- geo_df %>% 
    filter(cell_type == "PBMCs") %>% 
    dplyr::select("sample", "donor")

ground_truth_df <- geo_df %>% 
    dplyr::select(-sample) %>% 
    spread(key = "cell_type", value = "fraction") %>% 
    drop_na() %>% 
    inner_join(translation_df, .) %>% 
    dplyr::select(- c(donor, PBMCs))

query_df <-
    AnnotationDbi::select(
        hgu133plus2.db,
        keys=keys(hgu133plus2.db,keytype="PROBEID"),
        columns=c("SYMBOL"),
        keytype="PROBEID") %>%
    as_data_frame() %>%
    set_colnames(c("Probe", "Hugo")) %>%
    drop_na()

expr_df <- gse_object$GSE11057_series_matrix.txt.gz@assayData %>% 
    assayDataElement('exprs') %>%
    matrix_to_df("Probe") %>%
    inner_join(query_df) %>%
    dplyr::select(Hugo, everything()) %>%
    dplyr::select(-Probe) %>%
    group_by(Hugo) %>%
    summarise_all(mean) %>%
    drop_na() %>%
    filter(Hugo != "") %>% 
    gather(key = "sample", value = "expr", -Hugo) %>% 
    mutate(sample = str_remove_all(sample, "X"))

samples_in_common <- intersect(expr_df$sample, ground_truth_df$sample)

expr_df <- expr_df %>% 
    filter(sample %in% samples_in_common) %>% 
    spread(key = "sample", value = "expr")

ground_truth_df <- ground_truth_df %>% 
    filter(sample %in% samples_in_common)
  
activity_obj <- Activity(
    name = "create",
    description = "process GEO data into usable tables",
    used = list(),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE11057/create_processed_tables.R")
)

write_tsv(expr_df, "expression_affy.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")

upload_file_to_synapse("expression_affy.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("ground_truth.tsv", gt_upload_id, activity_obj = activity_obj)

