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
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE84697/"

expr_id        <- "syn12299856"
anno_id        <- "syn12299857"
genes_id       <- "syn11918430"
cs_results_id  <- "syn12299862"
mcp_results_id <- "syn12299863"



setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores())


gene_df <- create_df_from_synapse_id(genes_id)

mcp_genes <- gene_df %>% 
    filter(Method == "mcpcounter") %>%
    use_series("Hugo")

cs_genes <- gene_df %>% 
    filter(Method == "cibersort") %>%
    use_series("Hugo")

anno_df <- create_df_from_synapse_id(anno_id) %>% 
    arrange(sample)

log_matrix <- expr_id %>%
    create_df_from_synapse_id() %>% 
    df_to_matrix("Hugo") %>% 
    .[,order(colnames(.))] %>% 
    .[rowSums(.) > 0,]




# heatmaps --------------------------------------------------------------------

zscore_mat <- log_matrix %>% 
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

heatmap_col_df <- anno_df %>% 
    data.frame %>% 
    column_to_rownames("sample")

png('GSE84697_mcpcounter_genes_heatmap.png', width = 800, height = 1000)
pheatmap(
    mcp_zscore_matrix,
    main = "MCPCounter GSE84697",
    annotation_row = mcp_heatmap_row_df,
    annotation_col = heatmap_col_df,
    cluster_rows = F,
    scale = "none")
dev.off()

png('GSE84697_mcpcounter_genes_rows_clustered_heatmap.png', width = 800, height = 1000)
pheatmap(
    mcp_zscore_matrix,
    main = "MCPCounter GSE84697",
    annotation_row = mcp_heatmap_row_df,
    annotation_col = heatmap_col_df,
    scale = "none")
dev.off()

png('GSE84697_cibersort_genes_heatmap.png', width = 800, height = 3000)
pheatmap(
    cs_zscore_matrix,
    main = "Cibersort GSE84697",
    annotation_col = heatmap_col_df,
    scale = "none",
    fontsize = 15,
    fontsize_row = 5)
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
    full_join(anno_df, by = c("sample"))

mcp_result_df <- mcp_results_id %>%
    download_from_synapse %>% 
    read.table %>% 
    t %>% 
    .[,colSums(.) > 0] %>% 
    matrix_to_df("sample") %>% 
    set_colnames(str_replace_all(colnames(.), "\\.", " ")) %>% 
    gather("mcpcounter_cell_type", "predicted_score", `T cells`:Fibroblasts) %>% 
    full_join(anno_df, by = c("sample")) 

png('GSE84697_cibersort_facet_scatterplot.png', height = 1000)
ggplot(cs_result_df, aes(x = cibersort_cell_type, y = predicted_fraction)) +
    geom_point() +
    facet_grid(cell_type ~ .) +
    ylab("Predicted fraction") +
    xlab("Cibersort cell type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(strip.text.y = element_text(size = 10, angle = 0)) +
    ggtitle("Cibersort GSE84697")
dev.off()

png('GSE84697_mcpcounter_facet_scatterplot.png', height = 1000)
ggplot(mcp_result_df, aes(x = mcpcounter_cell_type, y = predicted_score)) +
    geom_point() +
    facet_grid(cell_type ~ .) +
    ylab("Predicted score") +
    xlab("MCPCounter cell type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(strip.text.y = element_text(size = 10, angle = 0)) +
    ggtitle("MCPCounter GSE84697")
dev.off()


# pca plots -------------------------------------------------------------------

pca_matrix <- t(log_matrix)

png('GSE84697_PCA.png', width = 800, height = 800)
autoplot(
    prcomp(pca_matrix), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "patient",
    size = 3,
    main = "GSE84697") +
    scale_shape_manual(values = 1:19) +
    theme_bw()
dev.off()

