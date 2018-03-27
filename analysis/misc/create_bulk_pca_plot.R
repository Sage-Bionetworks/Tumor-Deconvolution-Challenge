library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(ggfortify)


home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/misc/"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)

hugo_id       <- "syn11536071"
hugo_df <-  create_df_from_synapse_id(hugo_id)


# MTAB ------------------------------------------------------------------------

MTAB_count_id <- "syn11958709"
MTAB_anno_id  <- "syn11968317"

MTAB_anno_df <- MTAB_anno_id %>% 
    download_from_synapse %>% 
    fread(select = c("Comment[ENA_RUN]", "Factor Value[cell type]")) %>% 
    as_data_frame %>% 
    set_colnames(c("sample", "cell_type")) %>% 
    mutate(sample = str_c(sample, ".bam")) %>% 
    distinct %>% 
    arrange(sample) %>% 
    inset("dataset", value = "MTAB") 

MTAB_cpm_df <- MTAB_count_id %>%
    create_df_from_synapse_id(unzip = T) %>% 
    select(-c(Chr, Start, End, Strand, Length)) %>% 
    df_to_matrix("Geneid") %>% 
    add(1) %>% 
    apply(2, calculate_cpm) %>% 
    matrix_to_df("ensembl_gene_id") %>% 
    inner_join(hugo_df) %>% 
    select(hgnc_symbol, everything()) %>% 
    filter(hgnc_symbol != "") %>% 
    select(-ensembl_gene_id) %>% 
    group_by(hgnc_symbol) %>% 
    summarise_all(.funs = sum) %>% 
    ungroup 

# GSE62408 --------------------------------------------------------------------

GSE62408_id <- "syn11915424"

GSE62408_rpkm_df <- GSE62408_id %>% 
    create_df_from_synapse_id(unzip = T) %>% 
    dplyr::rename("hgnc_symbol" = `-`) %>% 
    group_by(hgnc_symbol) %>% 
    summarise_all(.funs = sum) %>% 
    ungroup 

GSE62408_anno_df <- data_frame(
    "sample" = colnames(GSE62408_rpkm_df)[-1],
    "cell_type" = colnames(GSE62408_rpkm_df)[-1],
    "dataset" = "GSE62408")

# GSE64655 --------------------------------------------------------------------

GSE64655_expr_id <- "syn11969378"
GSE64655_anno_id <- "syn11969387"

GSE64655_expr_df <- GSE64655_expr_id %>%
    create_df_from_synapse_id(unzip = T, skip = 3) %>% 
    select(-c(`Gene Type`, Description, `Gene Symbol`)) %>% 
    left_join(hugo_df, by = c("Gene ID" = "ensembl_gene_id")) %>% 
    select(-`Gene ID`) %>% 
    group_by(hgnc_symbol) %>% 
    summarise_all(sum) %>% 
    ungroup %>% 
    filter(hgnc_symbol != "") 

GSE64655_anno_df <- GSE64655_anno_id %>% 
    create_df_from_synapse_id(unzip = T, skip = 30, nrow = 7) %>%
    rename("title" = `!Sample_title`) %>% 
    filter(title == "!Sample_source_name_ch1") %>% 
    transpose_df("title", "sample") %>% 
    set_colnames(c("sample", "cell_type")) %>% 
    mutate(patient = str_sub(sample, end = 4)) %>% 
    mutate(days = str_sub(sample, start = -9, end = -9)) %>% 
    mutate(ABV1 = str_sub(sample, start = 6, end = -11)) %>% 
    mutate(ABV2 = ifelse(ABV1 == "DC", "mDC", 
                         ifelse(ABV1 == "MO", "Mono", 
                                ifelse(ABV1 == "NEU", "Neut", 
                                       ifelse(ABV1 == "PMBC", "PBMC", ABV1))))) %>% 
    mutate(sample = str_c(patient, "_", ABV2, "_d", days)) %>% 
    select(-c(ABV1, ABV2)) %>% 
    inset("dataset", value = "GSE64655") %>% 
    arrange(sample)


# combine ---------------------------------------------------------------------

exp_dfs <- list(MTAB_cpm_df, GSE62408_rpkm_df, GSE64655_expr_df)
anno_dfs <- list(MTAB_anno_df, GSE62408_anno_df, GSE64655_anno_df)

genes <- exp_dfs %>% 
    map(extract2, "hgnc_symbol") %>% 
    reduce(intersect)

annotation_df <- bind_rows(anno_dfs)

combined_matrix <- exp_dfs %>% 
    map(filter, hgnc_symbol %in% genes) %>% 
    map(transpose_df, "hgnc_symbol", "sample") %>% 
    map(arrange, sample) %>% 
    bind_rows %>% 
    df_to_matrix("sample") %>% 
    add(1) %>% 
    log10

pc_obj <- prcomp(combined_matrix)

png('bulk_rna_pca.png', height = 500, width = 500)
autoplot(pc_obj, data = annotation_df, colour = 'cell_type', shape = "dataset", main = "Bulk RNA PCA")
dev.off()

