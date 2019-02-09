library(tidyverse)
library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)
library(GEOquery)
library(AnnotationDbi)
library(hugene11sttranscriptcluster.db)

upload_id     <- "syn17091886"
gt_upload_id  <- "syn17091885"

source("../../scripts/utils.R")
synLogin()

dataset       <- "GSE40240"
github.path   <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/"
dataset.path  <- paste0(github.path, dataset)
script_url    <- paste0(dataset.path, "/create_processed_tables.R")
used          <- NA
activity_name <- "process GEO data into usable tables"

gse_object <- getGEO(dataset, GSEMatrix = TRUE)

geo_df <- gse_object %>% 
    extract2(1) %>% 
    phenoData() %>% 
    pData() %>%
    rownames_to_column("sample") %>%
    as_data_frame %>%
    dplyr::select(
        "sample", 
        "age" = "age:ch1", 
        "gender" = "Sex:ch1",
        "race" = "ethnicity:ch1",
        "basophils" =  "relative.basophils:ch1",                
        "eosinophils" = "relative.eosinophils:ch1",
        "lymphocytes" = "relative.lymphocytes:ch1",         
        "monocytes" = "relative.monocytes:ch1",
        "neutrophils" = "relative.neutrophils:ch1") 


query_df <-
    AnnotationDbi::select(
        hugene11sttranscriptcluster.db,
        keys=keys(hugene11sttranscriptcluster.db,keytype="PROBEID"),
        columns=c("SYMBOL"),
        keytype="PROBEID") %>%
    as_data_frame() %>%
    set_colnames(c("Probe", "Hugo")) %>%
    drop_na()

## Confirm that these data are RMA-normalized data--i.e., in log2 space
data.processing <-
    gse_object %>%
    extract2(1) %>%
    phenoData() %>%
    pData() %>%
    dplyr::pull(data_processing) %>%
    unique() %>%
    as.character()

cat(paste0(dataset, " data processing: ", data.processing))
if(grepl(data.processing, pattern="RMA")) {
  cat(paste0("Based on RMA processing, assuming data are in log2 space\n"))
} else {
  stop("Data were not processed via RMA.  Can not confirm in log2 space\n")
}

expr_df <- gse_object$GSE40240_series_matrix.txt.gz@assayData %>% 
    assayDataElement('exprs') %>%
    matrix_to_df("Probe") %>%
    inner_join(query_df) %>%
    dplyr::select(Hugo, everything()) %>%
    dplyr::select(-Probe) %>%
    drop_na() %>%
    filter(Hugo != "") %>% 
    group_by(Hugo) %>%
    summarise_all(mean) %>%
    gather(key = "sample", value = "expr", -Hugo) %>% 
    mutate(sample = str_remove_all(sample, "X"))

samples_in_common <- intersect(expr_df$sample, geo_df$sample)

expr_df <- expr_df %>% 
    filter(sample %in% samples_in_common) %>% 
    spread(key = "sample", value = "expr")

linear_expr_df <- expr_df %>%
    df_to_matrix("Hugo") %>%
    raise_to_power(x=2, power=.) %>%
    matrix_to_df("Hugo")

ground_truth_df <- geo_df %>%
    dplyr::select(-c(age, gender, race)) %>% 
    filter(sample %in% samples_in_common)

annotation_df <- geo_df %>% 
    dplyr::select(sample, age, gender, race) %>% 
    filter(sample %in% samples_in_common)

gpl <- getGEO(gse_object[[1]]@annotation, destdir=".")
microarray_type <- NA
if(Meta(gpl)$title == "[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array") {
  microarray_type <- "Affymetrix HG-U133 Plus 2.0"
} else if(Meta(gpl)$title == "[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]") {
  microarray_type <- "Affymetrix Human Gene 1.0 ST Array"
} else {
  stop(paste0("Unknown array type ", Meta(gpl)$title))
}

expression_manifest_df <- tibble(
    path = c("expression_log.tsv", "expression_linear.tsv"),
    parent = upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "expression",
    expression_type = "microarray",
    microarray_type = microarray_type,
    expression_space = c("log2", "linear")
)

annotation_manifest_df <- tibble(
    path = "annotation.tsv",
    parent = upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "annotations",
    annotations = str_c(colnames(annotation_df)[-1], collapse = ";")
)

ground_truth_manifest_df <- tibble(
    path = "ground_truth.tsv",
    parent = gt_upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "ground truth",
    unit = "fraction",
    cell_types = str_c(colnames(ground_truth_df)[-1], collapse = ";")
)

write_tsv(expr_df, "expression_log.tsv")
write_tsv(linear_expr_df, "expression_linear.tsv")
write_tsv(annotation_df, "annotation.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")

write_tsv(expression_manifest_df, "expression_manifest.tsv")
write_tsv(annotation_manifest_df, "annotation_manifest.tsv")
write_tsv(ground_truth_manifest_df, "ground_truth_manifest.tsv")

syncToSynapse("expression_manifest.tsv")
syncToSynapse("annotation_manifest.tsv")
syncToSynapse("ground_truth_manifest.tsv")

