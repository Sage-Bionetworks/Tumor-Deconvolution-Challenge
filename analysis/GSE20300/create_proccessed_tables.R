library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(GEOquery)
library(xlsx)
library(AnnotationDbi)
library(hgu133plus2.db)

ground_truth_id <- "syn17089680"

upload_id     <- "syn17026189"
gt_upload_id  <- "syn17026188"

source("../../scripts/utils.R")
synLogin()


ground_truth_df <- ground_truth_id %>% 
    download_from_synapse() %>% 
    xlsx::read.xlsx(1)

gse <- getGEO("GSE20300", GSEMatrix = TRUE)

series_df <-
    pData(phenoData(gse[[1]])) %>%
    rownames_to_column("sample") %>%
    as_data_frame %>% 
    arrange(`transplant state:ch1`)

translation_df <- 
    bind_cols(ground_truth_df, series_df) %>% 
    dplyr::select(sample, title, Sample.ID)

ground_truth_df <- ground_truth_df %>% 
    left_join(translation_df) %>% 
    dplyr::select(-c(Sample.ID, Patient.Group, Total, title)) %>% 
    dplyr::select(sample, everything())

annotation_df <- series_df %>% 
    dplyr::select(sample, `transplant state:ch1`) %>% 
    dplyr::rename(transplant_state = `transplant state:ch1`)



query_df <-
    AnnotationDbi::select(
        hgu133plus2.db,
        keys=keys(hgu133plus2.db,keytype="PROBEID"),
        columns=c("SYMBOL"),
        keytype="PROBEID") %>%
    as_data_frame() %>%
    set_colnames(c("PROBE", "Hugo")) %>%
    drop_na()

expr_df <-
    assayDataElement(gse$GSE20300_series_matrix.txt.gz@assayData, 'exprs') %>% 
    matrix_to_df("PROBE") %>%
    inner_join(query_df) %>%
    dplyr::select(Hugo, everything()) %>%
    dplyr::select(-PROBE) %>%
    group_by(Hugo) %>%
    summarise_all(mean) %>%
    drop_na() %>%
    filter(Hugo != "") %>% 
    gather(key = "sample", value = "expr", -Hugo)

samples_in_common <- 
    purrr::reduce(list(expr_df$sample, ground_truth_df$sample, annotation_df$sample, translation_df$sample),
           intersect)


activity_obj <- Activity(
    name = "create",
    description = "process GEO data into usable tables",
    used = list(ground_truth_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE20300/create_processed_tables.R")
)


write_tsv(expr_df, "expression_affy.tsv")
write_tsv(annotation_df, "annotation.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")
write_tsv(translation_df, "translation.tsv")

upload_file_to_synapse("expression_affy.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("ground_truth.tsv", gt_upload_id, activity_obj = activity_obj)
upload_file_to_synapse("translation.tsv", upload_id, activity_obj = activity_obj)

