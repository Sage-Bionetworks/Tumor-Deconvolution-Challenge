library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(pheatmap)
library(preprocessCore)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE62408/"

count_id      <- "syn11915424"
genes_id      <- "syn11918430"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores())

log_mat <- create_df_from_synapse_id(count_id, unzip = T) %>% 
    dplyr::rename("Hugo" = `-`) %>% 
    group_by(Hugo) %>% 
    summarise_all(mean) %>% 
    df_to_matrix("Hugo") %>%
    add(1) %>% 
    log10 

gene_df <-  create_df_from_synapse_id(genes_id)


mcp_genes <- gene_df %>% 
    filter(Method == "mcpcounter") %>%
    use_series("Hugo")

cs_genes <- gene_df %>% 
    filter(Method == "cibersort") %>%
    use_series("Hugo")

mcp_log_matrix <- log_mat %>% 
    .[rownames(.) %in% mcp_genes,]

cs_log_matrix <- log_mat %>% 
    .[rownames(.) %in% cs_genes,]

mcp_heatmap_row_df <- gene_df %>% 
    filter(Method == "mcpcounter") %>% 
    filter(Hugo %in% rownames(mcp_log_matrix)) %>% 
    select(-Method) %>% 
    arrange(cell_type) %>% 
    data.frame %>% 
    column_to_rownames("Hugo") %>% 
    set_names("Cell Type")

cs_zscore_matrix <- normalize.quantiles(cs_log_matrix) %>%
    set_colnames(colnames(cs_log_matrix)) %>%
    set_rownames(rownames(cs_log_matrix)) %>% 
    zscore_matrix

mcp_zscore_matrix <- normalize.quantiles(mcp_log_matrix) %>%
    set_colnames(colnames(mcp_log_matrix)) %>%
    set_rownames(rownames(mcp_log_matrix)) %>%
    .[ rownames(mcp_heatmap_row_df),] %>% 
    zscore_matrix

png('GSE62408_mcpcounter_genes_heatmap.png', width = 4000, height = 4000)
pheatmap(
    mcp_zscore_matrix,
    main = "MCPCounter GSE62408",
    annotation_row = mcp_heatmap_row_df,
    cluster_rows = F,
    scale = "none")
dev.off()

png('GSE62408_mcpcounter_genes_rows_clustered_heatmap.png', width = 4000, height = 4000)
pheatmap(
    mcp_zscore_matrix,
    main = "MCPCounter GSE62408",
    annotation_row = mcp_heatmap_row_df,
    scale = "none")
dev.off()

png('GSE62408_cibersort_genes_heatmap.png', width = 4000, height = 4000)
pheatmap(
    cs_zscore_matrix,
    main = "Cibersort GSE62408",
    scale = "none")
dev.off()

# -----------------------------------------------------------------------------

cs_result_df <- "cs_results.tsv" %>%
    str_c(tmp_dir, .) %>%
    read_tsv %>%
    select(-c(`P.value`, `RMSE`, Correlation)) %>% 
    df_to_matrix("cell_type") %>% 
    .[,colSums(.) > 0] %>% 
    matrix_to_df("cell_type") %>% 
    set_colnames(str_replace_all(colnames(.), "\\.", " ")) %>% 
    gather("cibersort_cell_type", "predicted_fraction", `B cells naive`:Neutrophils)
    

png('GSE62408_cibersort_facet_scatterplot.png', height = 1000)
ggplot(cs_result_df, aes(x = cibersort_cell_type, y = predicted_fraction)) +
    geom_point() +
    facet_grid(cell_type ~ .) +
    ylab("Predicted fraction") +
    xlab("Cibersort cell type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(strip.text.y = element_text(size = 16)) +
    ggtitle("Cibersort GSE62408")
dev.off()

mcp_result_df <- "mcp_results.tsv" %>%
    str_c(tmp_dir, .) %>%
    read_tsv %>% 
    transpose_df("cell_type", "cell_type") %>% 
    df_to_matrix("cell_type") %>% 
    .[,colSums(.) > 0] %>% 
    matrix_to_df("cell_type") %>% 
    set_colnames(str_replace_all(colnames(.), "\\.", " ")) %>% 
    gather("mpcounter_cell_type", "predicted_score", `T cells`:Fibroblasts) 

png('GSE62408_mcpcounter_facet_scatterplot.png', height = 1000)
ggplot(mcp_result_df, aes(x = mpcounter_cell_type, y = predicted_score)) +
    geom_point() +
    ylab("Predicted score") +
    xlab("MCPCounter cell type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    facet_grid(cell_type ~ .) +
    theme(strip.text.y = element_text(size = 16)) +
    ggtitle("MCPCounter GSE62408")
dev.off()

# pca plots -------------------------------------------------------------------

pca_matrix <- t(log_mat)
annotation_df <- pca_matrix %>% 
    matrix_to_df("sample") %>% 
    select(sample)


png('GSE62408_PCA.png', height = 1000)
autoplot(
    prcomp(pca_matrix), 
    data = annotation_df, 
    shape = "sample", 
    size = 4,
    main = "GSE62408") +
    scale_shape_manual(values = 16:22) +
    theme_bw()
dev.off()