library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(preprocessCore)
library(ggfortify)
library(sva)
library(BiocParallel)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE89225/"

expr_id        <- "syn12063105"
anno_id        <- "syn12063109"


setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

anno_df <- create_df_from_synapse_id(anno_id) %>% 
    arrange(sample)

gene_matrix <- expr_id %>%
    create_df_from_synapse_id() %>% 
    group_by(Hugo) %>% 
    summarise_all(sum) %>% 
    ungroup %>% 
    df_to_matrix("Hugo") %>% 
    .[,order(colnames(.))] %>% 
    .[!rowSums(.) == 0,] %>% 
    calculate_cpm %>% 
    add(1) %>% 
    log10 

par <- MulticoreParam(workers = 7)

combat_platform <-  ComBat(
    gene_matrix,
    batch = anno_df$platform, 
    mod = model.matrix(~1, data = anno_df),
    BPPARAM = par)

combat_source <-  ComBat(
    gene_matrix,
    batch = anno_df$source, 
    mod = model.matrix(~1, data = anno_df),
    BPPARAM = par)

autoplot(
    prcomp(t(gene_matrix)), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "platform",
    size = 3,
    main = "GSE89225") +
    scale_shape_manual(values = 21:25) +
    theme_bw()

autoplot(
    prcomp(t(combat_platform)), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "platform",
    size = 3,
    main = "GSE89225") +
    scale_shape_manual(values = 1:19) +
    theme_bw()

autoplot(
    prcomp(t(combat_source)), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "platform",
    size = 3,
    main = "GSE89225") +
    scale_shape_manual(values = 1:19) +
    theme_bw()








combat_platform_patient <-  ComBat(
    combat_platform,
    batch = anno_df$patient, 
    mod = model.matrix(~1, data = anno_df),
    BPPARAM = par)





log_matrix_platform <- combat_platform %>% 
    calculate_cpm %>% 
    add(1) %>% 
    log10 %>% 
    t

log_matrix_platform[is.nan(log_matrix_platform)] <- 0


# pc_obj <- prcomp(log_matrix)
# 
# pc_df <- pc_obj %>%
#     use_series(x) %>%
#     matrix_to_df("sample") %>%
#     select(sample, PC1, PC2) %>%
#     left_join(anno_df)
# 
# x %>% ggplot() +
#     aes(x = PC1, y = PC2) +
#     geom_point(size = 5, aes(color = platform, shape = patient)) +
#     geom_point(size = 2, aes(color = source)) +
#     theme_bw() +
#     scale_shape_manual(values = 0:20) +
#     scale_colour_manual(values = c("red", "purple", "black", "blue", "orange"))
    

autoplot(
    prcomp(log_matrix), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "platform",
    size = 3,
    main = "GSE89225") +
    scale_shape_manual(values = 21:25) +
    theme_bw()

autoplot(
    prcomp(log_matrix_platform), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "platform",
    size = 3,
    main = "GSE89225") +
    scale_shape_manual(values = 1:19) +
    theme_bw()

autoplot(
    prcomp(log_matrix), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "source",
    size = 3,
    main = "GSE89225") +
    scale_shape_manual(values = 1:19) +
    theme_bw()

autoplot(
    prcomp(log_matrix_platform), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "source",
    size = 3,
    main = "GSE89225") +
    scale_shape_manual(values = 1:19) +
    theme_bw()

autoplot(
    prcomp(log_matrix), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "patient",
    size = 3,
    main = "GSE89225") +
    scale_shape_manual(values = 1:19) +
    theme_bw()

autoplot(
    prcomp(log_matrix_platform), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "patient",
    size = 3,
    main = "GSE89225") +
    scale_shape_manual(values = 1:19) +
    theme_bw()





log_matrix_platform_patient <- combat_platform_patient %>% 
    calculate_cpm %>% 
    add(1) %>% 
    log10 %>% 
    t

log_matrix_platform_patient[is.nan(log_matrix_platform_patient)] <- 0

autoplot(
    prcomp(log_matrix_platform_patient), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "platform",
    size = 3,
    main = "GSE89225") +
    scale_shape_manual(values = 21:25) +
    theme_bw()

autoplot(
    prcomp(log_matrix_platform), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "platform",
    size = 3,
    main = "GSE89225") +
    scale_shape_manual(values = 1:19) +
    theme_bw()

autoplot(
    prcomp(log_matrix), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "source",
    size = 3,
    main = "GSE89225") +
    scale_shape_manual(values = 1:19) +
    theme_bw()

autoplot(
    prcomp(log_matrix_platform), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "source",
    size = 3,
    main = "GSE89225") +
    scale_shape_manual(values = 1:19) +
    theme_bw()

autoplot(
    prcomp(log_matrix), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "patient",
    size = 3,
    main = "GSE89225") +
    scale_shape_manual(values = 1:19) +
    theme_bw()

autoplot(
    prcomp(log_matrix_platform), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "patient",
    size = 3,
    main = "GSE89225") +
    scale_shape_manual(values = 1:19) +
    theme_bw()




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

png('GSE89225_mcpcounter_genes_heatmap.png', width = 800, height = 1000)
pheatmap(
    mcp_zscore_matrix,
    main = "MCPCounter GSE89225",
    annotation_row = mcp_heatmap_row_df,
    annotation_col = heatmap_col_df,
    cluster_rows = F,
    scale = "none")
dev.off()

png('GSE89225_mcpcounter_genes_rows_clustered_heatmap.png', width = 800, height = 1000)
pheatmap(
    mcp_zscore_matrix,
    main = "MCPCounter GSE89225",
    annotation_row = mcp_heatmap_row_df,
    annotation_col = heatmap_col_df,
    scale = "none")
dev.off()

png('GSE89225_cibersort_genes_heatmap.png', width = 800, height = 3000)
pheatmap(
    cs_zscore_matrix,
    main = "Cibersort GSE89225",
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

png('GSE89225_cibersort_facet_scatterplot.png', height = 1000)
ggplot(cs_result_df, aes(x = cibersort_cell_type, y = predicted_fraction)) +
    geom_point() +
    facet_grid(cell_type ~ .) +
    ylab("Predicted fraction") +
    xlab("Cibersort cell type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(strip.text.y = element_text(size = 10, angle = 0)) +
    ggtitle("Cibersort GSE89225")
dev.off()

png('GSE89225_mcpcounter_facet_scatterplot.png', height = 1000)
ggplot(mcp_result_df, aes(x = mcpcounter_cell_type, y = predicted_score)) +
    geom_point() +
    facet_grid(cell_type ~ .) +
    ylab("Predicted score") +
    xlab("MCPCounter cell type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(strip.text.y = element_text(size = 10, angle = 0)) +
    ggtitle("MCPCounter GSE89225")
dev.off()


# pca plots -------------------------------------------------------------------

pca_matrix <- t(log_matrix)

png('GSE89225_patient_PCA.png', width = 800, height = 800)
autoplot(
    prcomp(pca_matrix), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "patient",
    size = 3,
    main = "GSE89225") +
    scale_shape_manual(values = 1:19) +
    theme_bw()
dev.off()

png('GSE89225_source_PCA.png', width = 800, height = 800)
autoplot(
    prcomp(pca_matrix), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "source",
    size = 3,
    main = "GSE89225") +
    scale_shape_manual(values = 1:9) +
    theme_bw()
dev.off()


png('GSE89225_platform_PCA.png', width = 800, height = 800)
autoplot(
    prcomp(pca_matrix), 
    data = anno_df, 
    shape = "cell_type", 
    colour = "platform",
    size = 3,
    main = "GSE89225") +
    scale_shape_manual(values = 1:9) +
    theme_bw()
dev.off()