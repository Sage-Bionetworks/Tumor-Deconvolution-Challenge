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
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE103322/"

expr_id        <- "syn11990376"
anno_id        <- "syn11990378"
genes_id       <- "syn11918430"


setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores())

annotation_df <- create_df_from_synapse_id(anno_id) %>% 
    filter(cell_source != "Lymph node") %>% 
    group_by(cell_type) %>% 
    slice(1:15) %>% 
    ungroup

gene_df <-  create_df_from_synapse_id(genes_id)

mcp_genes <- gene_df %>% 
    filter(Method == "mcpcounter") %>%
    use_series("Hugo")

cs_genes <- gene_df %>% 
    filter(Method == "cibersort") %>%
    use_series("Hugo")

expr_mat <- expr_id %>%
    create_df_from_synapse_id() %>%
    df_to_matrix("Hugo") %>% 
    .[,colnames(.) %in% anno_df$sample] %>% 
    .[rownames(.) %in% mcp_genes,]

mcp_heatmap_row_df <- gene_df %>% 
    filter(Method == "mcpcounter") %>% 
    filter(Hugo %in% rownames(expr_mat)) %>% 
    select(-Method) %>% 
    arrange(cell_type) %>% 
    data.frame %>% 
    column_to_rownames("Hugo") %>% 
    set_names("Cell Type")

heatmap_col_df <- annotation_df %>% 
    data.frame %>% 
    column_to_rownames("sample")

expr_mat <- expr_mat %>% 
    .[rownames(mcp_heatmap_row_df),] %>% 
    .[,rownames(heatmap_col_df)]


png('GSE103322_mcpcounter_genes_heatmap.png', width = 4000, height = 4000)
pheatmap(
    expr_mat,
    main = "MCPCounter GSE103322",
    annotation_row = mcp_heatmap_row_df,
    annotation_col = heatmap_col_df,
    cluster_rows = F,
    cluster_cols = F,
    scale = "none")
dev.off()

# heatmaps --------------------------------------------------------------------

zscore_mat <- log_mat %>% 
    quantile_normalize_matrix %>% 
    zscore_matrix

mcp_zscore_matrix <- zscore_mat %>% 
    .[rownames(.) %in% mcp_genes,]

cs_zscore_matrix <- zscore_mat %>% 
    .[rownames(.) %in% cs_genes,]



mcp_zscore_matrix <-  mcp_zscore_matrix[rownames(mcp_heatmap_row_df),]





png('GSE103322_mcpcounter_genes_rows_clustered_heatmap.png', width = 4000, height = 4000)
pheatmap(
    mcp_zscore_matrix,
    main = "MCPCounter GSE103322",
    annotation_row = mcp_heatmap_row_df,
    annotation_col = heatmap_col_df,
    scale = "none")
dev.off()

png('GSE103322_cibersort_genes_heatmap.png', width = 4000, height = 4000)
pheatmap(
    cs_zscore_matrix,
    main = "Cibersort GSE103322",
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

png('GSE103322_cibersort_facet_scatterplot.png', height = 1000)
ggplot(cs_result_df, aes(x = cibersort_cell_type, y = predicted_fraction)) +
    geom_point() +
    facet_grid(cell_type ~ .) +
    ylab("Predicted fraction") +
    xlab("Cibersort cell type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(strip.text.y = element_text(size = 10, angle = 0)) +
    ggtitle("Cibersort GSE103322")
dev.off()

png('GSE103322_mcpcounter_facet_scatterplot.png', height = 1000)
ggplot(mcp_result_df, aes(x = mcpcounter_cell_type, y = predicted_score)) +
    geom_point() +
    facet_grid(cell_type ~ .) +
    ylab("Predicted score") +
    xlab("MCPCounter cell type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(strip.text.y = element_text(size = 10, angle = 0)) +
    ggtitle("MCPCounter GSE103322")
dev.off()


# pca plots -------------------------------------------------------------------

pca_matrix <- t(log_mat)

png('GSE103322_PCA.png')
autoplot(
    prcomp(pca_matrix), 
    data = annotation_df, 
    shape = "cell_type", 
    colour = "batch",
    size = 3,
    main = "GSE103322") +
    scale_shape_manual(values = 16:19) +
    theme_bw()
dev.off()












