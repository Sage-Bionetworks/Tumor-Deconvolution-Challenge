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
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE83115/"

expr_id        <- "syn11977964"
genes_id       <- "syn11918430"
cs_results_id  <- "syn11978093"
mcp_results_id <- "syn11978096"



setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores())


gene_df <-  create_df_from_synapse_id(genes_id)

expr_mat <- expr_id %>%
    create_df_from_synapse_id() %>% 
    select(- Ensembl) %>% 
    group_by(Hugo) %>% 
    summarise_all(sum) %>% 
    ungroup %>% 
    filter(Hugo != "") %>% 
    df_to_matrix("Hugo") 

log_mat <- expr_mat %>% 
    .[rowSums(.) > 0,] %>% 
    add(1) %>% 
    log10

annotation_df <- expr_mat %>% 
    t %>% 
    matrix_to_df("sample") %>% 
    select(sample) %>% 
    mutate(cell_type = str_sub(sample, end = -5)) %>% 
    mutate(batch = str_sub(sample, start = -3))

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

png('GSE83115_mcpcounter_genes_heatmap.png', width = 1000, height = 1000)
pheatmap(
    mcp_zscore_matrix,
    main = "MCPCounter GSE83115",
    annotation_row = mcp_heatmap_row_df,
    annotation_col = heatmap_col_df,
    cluster_rows = F,
    scale = "none")
dev.off()

png('GSE83115_mcpcounter_genes_rows_clustered_heatmap.png', width = 1000, height = 1000)
pheatmap(
    mcp_zscore_matrix,
    main = "MCPCounter GSE83115",
    annotation_row = mcp_heatmap_row_df,
    annotation_col = heatmap_col_df,
    scale = "none")
dev.off()

png('GSE83115_cibersort_genes_heatmap.png', width = 1000, height = 1000)
pheatmap(
    cs_zscore_matrix,
    main = "Cibersort GSE83115",
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
    gather("cibersort_cell_type", "predicted_fraction", `B cells naive`:Eosinophils) %>% 
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

png('GSE83115_cibersort_facet_scatterplot.png', height = 1000)
ggplot(cs_result_df, aes(x = cibersort_cell_type, y = predicted_fraction)) +
    geom_point() +
    facet_grid(cell_type ~ .) +
    ylab("Predicted fraction") +
    xlab("Cibersort cell type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(strip.text.y = element_text(size = 10, angle = 0)) +
    ggtitle("Cibersort GSE83115")
dev.off()

png('GSE83115_mcpcounter_facet_scatterplot.png', height = 1000)
ggplot(mcp_result_df, aes(x = mcpcounter_cell_type, y = predicted_score)) +
    geom_point() +
    facet_grid(cell_type ~ .) +
    ylab("Predicted score") +
    xlab("MCPCounter cell type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(strip.text.y = element_text(size = 10, angle = 0)) +
    ggtitle("MCPCounter GSE83115")
dev.off()


# pca plots -------------------------------------------------------------------

pca_matrix <- t(log_mat)

png('GSE83115_PCA.png', height = 700, width = 700)
autoplot(
    prcomp(pca_matrix), 
    data = annotation_df, 
    shape = "cell_type", 
    colour = "batch",
    size = 3,
    main = "GSE83115") +
    scale_shape_manual(values = 16:19) +
    theme_bw()
dev.off()












