library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(pheatmap)
library(preprocessCore)
library(ggfortify)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE64655/"

expr_id        <- "syn11973634"
annotation_id  <- "syn11973635"
genes_id       <- "syn11918430"
cs_results_id  <- "syn11969430"
mcp_results_id <- "syn11969431"



setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores())

annotation_df <- create_df_from_synapse_id(annotation_id)
gene_df <-  create_df_from_synapse_id(genes_id)

log_mat <- expr_id %>%
    create_df_from_synapse_id() %>% 
    select(- Ensembl) %>% 
    group_by(Hugo) %>% 
    summarise_all(sum) %>% 
    ungroup %>% 
    filter(Hugo != "") %>% 
    df_to_matrix("Hugo") %>% 
    .[rowSums(.) > 0,] %>% 
    add(1) %>% 
    log10

mcp_genes <- gene_df %>% 
    filter(Method == "mcpcounter") %>%
    use_series("Hugo")

cs_genes <- gene_df %>% 
    filter(Method == "cibersort") %>%
    use_series("Hugo")


# heatmaps --------------------------------------------------------------------

zscore_mat <- log_mat %>% 
    quantile_normalize_matrix %>% 
    zscore_matrix

mcp_zscore_matrix <- zscore_mat %>% 
    .[rownames(.) %in% mcp_genes,]

cs_zscore_matrix <- zscore_mat %>% 
    .[rownames(.) %in% cs_genes,]

mcp_heatmap_row_df <- gene_df %>% 
    filter(Method == "mcpcounter") %>% 
    filter(Hugo %in% rownames(mcp_zscore_matrix)) %>% 
    select(-Method) %>% 
    arrange(cell_type) %>% 
    data.frame %>% 
    column_to_rownames("Hugo") %>% 
    set_names("Cell Type")

mcp_zscore_matrix <-  mcp_zscore_matrix[rownames(mcp_heatmap_row_df),]

heatmap_col_df <- annotation_df %>% 
    data.frame %>% 
    column_to_rownames("sample")

png('GSE64655_mcpcounter_genes_heatmap.png', width = 4000, height = 4000)
pheatmap(
    mcp_zscore_matrix,
    main = "MCPCounter GSE64655",
    annotation_row = mcp_heatmap_row_df,
    annotation_col = heatmap_col_df,
    cluster_rows = F,
    scale = "none")
dev.off()

png('GSE64655_mcpcounter_genes_rows_clustered_heatmap.png', width = 4000, height = 4000)
pheatmap(
    mcp_zscore_matrix,
    main = "MCPCounter GSE64655",
    annotation_row = mcp_heatmap_row_df,
    annotation_col = heatmap_col_df,
    scale = "none")
dev.off()

png('GSE64655_cibersort_genes_heatmap.png', width = 4000, height = 4000)
pheatmap(
    cs_zscore_matrix,
    main = "Cibersort GSE64655",
    annotation_col = heatmap_col_df,
    scale = "none")
dev.off()

# scatter plots ---------------------------------------------------------------

cs_result_df <- cs_results_id %>%
    create_df_from_synapse_id %>% 
    select(-c(`P-value`, `RMSE`, Correlation)) %>% 
    df_to_matrix("sample") %>% 
    .[,colSums(.) > 0] %>% 
    matrix_to_df("sample") %>% 
    set_colnames(str_replace_all(colnames(.), "\\.", " ")) %>% 
    gather("cibersort_cell_type", "predicted_fraction", `B cells naive`:Neutrophils) %>% 
    full_join(annotation_df, by = c("sample"))

mcp_result_df <- mcp_results_id %>%
    download_from_synapse %>% 
    read.table %>% 
    t %>% 
    .[,colSums(.) > 0] %>% 
    matrix_to_df("sample") %>% 
    set_colnames(str_replace_all(colnames(.), "\\.", " ")) %>% 
    gather("mcpcounter_cell_type", "predicted_score", `T cells`:Fibroblasts) %>% 
    full_join(annotation_df, by = c("sample")) 

png('GSE64655_cibersort_facet_scatterplot.png', height = 1000)
ggplot(cs_result_df, aes(x = cibersort_cell_type, y = predicted_fraction)) +
    geom_point() +
    facet_grid(cell_type ~ .) +
    ylab("Predicted fraction") +
    xlab("Cibersort cell type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(strip.text.y = element_text(size = 10, angle = 0)) +
    ggtitle("Cibersort GSE64655")
dev.off()

png('GSE64655_mcpcounter_facet_scatterplot.png', height = 1000)
ggplot(mcp_result_df, aes(x = mcpcounter_cell_type, y = predicted_score)) +
    geom_point() +
    facet_grid(cell_type ~ .) +
    ylab("Predicted score") +
    xlab("MCPCounter cell type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(strip.text.y = element_text(size = 10, angle = 0)) +
    ggtitle("MCPCounter GSE64655")
dev.off()


# pca plots -------------------------------------------------------------------

pca_matrix <- t(log_mat)

png('GSE64655_PCA.png', height = 1000)
autoplot(
    prcomp(pca_matrix), 
    data = annotation_df, 
    shape = "cell_type", 
    colour = "patient",
    size = 3,
    main = "GSE64655") +
    scale_shape_manual(values = 16:22) +
    theme_bw()
dev.off()












