library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(ggfortify)


home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/misc/"

synLogin()
setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)

# MTAB ------------------------------------------------------------------------

MTAB_count_id <- "syn12031262"
MTAB_anno_id  <- "syn11968317"

MTAB_anno_df <- MTAB_anno_id %>% 
    create_df_from_synapse_id %>% 
    arrange(sample) %>% 
    inset("dataset", value = "MTAB") 

MTAB_cpm_df <- MTAB_count_id %>%
    create_df_from_synapse_id %>% 
    select(-ensembl_gene_id) %>% 
    group_by(hgnc_symbol) %>% 
    summarise_all(.funs = sum) %>% 
    filter(hgnc_symbol != "") %>% 
    df_to_matrix("hgnc_symbol") %>% 
    add(1) %>% 
    apply(2, calculate_cpm) %>% 
    matrix_to_df("Hugo")


# GSE62408 --------------------------------------------------------------------

GSE62408_id <- "syn11915424"

GSE62408_rpkm_df <- GSE62408_id %>% 
    create_df_from_synapse_id(unzip = T) %>% 
    dplyr::rename("Hugo" = `-`) %>% 
    group_by(Hugo) %>% 
    summarise_all(.funs = sum) %>% 
    ungroup %>% 
    .[,order(colnames(.))] %>% 
    select(Hugo, everything())

GSE62408_anno_df <- data_frame(
    "sample" = colnames(GSE62408_rpkm_df)[-1],
    "cell_type" = colnames(GSE62408_rpkm_df)[-1],
    "dataset" = "GSE62408") %>% 
    arrange(sample) %>% 
    inset("dataset", value = "GSE62408") 

# GSE64655 --------------------------------------------------------------------

GSE64655_expr_id <- "syn11973634"
GSE64655_anno_id <- "syn11973635"

GSE64655_expr_df <- GSE64655_expr_id %>%
    create_df_from_synapse_id %>% 
    select(-Ensembl) %>% 
    group_by(Hugo) %>% 
    summarise_all(.funs = sum) %>% 
    filter(Hugo != "")

GSE64655_anno_df <- GSE64655_anno_id %>% 
    create_df_from_synapse_id %>%
    arrange(sample) %>% 
    inset("dataset", value = "GSE64655") 


# combine ---------------------------------------------------------------------

exp_dfs <- list(MTAB_cpm_df, GSE62408_rpkm_df, GSE64655_expr_df)
anno_dfs <- list(MTAB_anno_df, GSE62408_anno_df, GSE64655_anno_df)

genes <- exp_dfs %>% 
    map(extract2, "Hugo") %>% 
    reduce(intersect)

annotation_df <- bind_rows(anno_dfs)

combined_matrix <- exp_dfs %>% 
    map(filter, Hugo %in% genes) %>% 
    map(transpose_df, "Hugo", "sample") %>% 
    map(arrange, sample) %>% 
    bind_rows %>% 
    df_to_matrix("sample") %>% 
    add(1) %>% 
    log10

pc_obj <- prcomp(combined_matrix)

png('bulk_rna_pca.png', height = 500, width = 500)
autoplot(pc_obj, data = annotation_df, colour = 'cell_type', shape = "dataset", main = "Bulk RNA PCA")
dev.off()

