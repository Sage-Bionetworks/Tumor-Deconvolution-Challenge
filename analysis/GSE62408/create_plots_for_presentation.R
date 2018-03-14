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
cibersort_id  <- "syn11958451"
mcpcounter_id <- "syn11958450"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores())

rpkm_df <- create_df_from_synapse_id(count_id, unzip = T) %>% 
    dplyr::rename("Hugo" = `-`)

gene_df <-  create_df_from_synapse_id(genes_id)

cs_translation_df <- create_df_from_synapse_id(cibersort_id)
mcp_translation_df <- create_df_from_synapse_id(mcpcounter_id)

mcp_genes <- gene_df %>% 
    filter(Method == "mcpcounter") %>%
    use_series("Hugo")

cs_genes <- gene_df %>% 
    filter(Method == "cibersort") %>%
    use_series("Hugo")

mcp_zscore_matrix <- rpkm_df %>% 
    filter(Hugo %in% mcp_genes) %>% 
    group_by(Hugo) %>% 
    summarise_all(mean) %>% 
    df_to_matrix("Hugo") %>%
    add(1) %>% 
    log10 

cs_zscore_matrix <- rpkm_df %>% 
    filter(Hugo %in% cs_genes) %>% 
    group_by(Hugo) %>% 
    summarise_all(mean) %>% 
    df_to_matrix("Hugo") %>% 
    add(1) %>% 
    log10

mcp_heatmap_row_df <- gene_df %>% 
    filter(Method == "mcpcounter") %>% 
    filter(Hugo %in% rownames(mcp_zscore_matrix)) %>% 
    select(-Method) %>% 
    arrange(cell_type) %>% 
    data.frame %>% 
    column_to_rownames("Hugo") %>% 
    set_names("Cell Type")

cs_quant_matrix <- normalize.quantiles(cs_zscore_matrix) %>%
    set_colnames(colnames(cs_zscore_matrix)) %>%
    set_rownames(rownames(cs_zscore_matrix)) %>% 
    zscore_matrix

mcp_quant_matrix <- normalize.quantiles(mcp_zscore_matrix) %>%
    set_colnames(colnames(mcp_zscore_matrix)) %>%
    set_rownames(rownames(mcp_zscore_matrix)) %>%
    .[ rownames(mcp_heatmap_row_df),] %>% 
    zscore_matrix

png('GSE62408_mcpcounter_genes_heatmap.png', width = 4000, height = 4000)
pheatmap(
    mcp_quant_matrix,
    main = "MCPCounter GSE62408",
    annotation_row = mcp_heatmap_row_df,
    cluster_rows = F,
    scale = "none")
dev.off()

png('GSE62408_mcpcounter_genes_rows_clustered_heatmap.png', width = 4000, height = 4000)
pheatmap(
    mcp_quant_matrix,
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



# 
# 
# cs_result_df <- "cs_results.tsv" %>%
#     str_c(tmp_dir, .) %>%
#     read_tsv %>%
#     transpose_df("cell_type", "cibersort_cell_type") %>% 
#     gather("actual_cell_type", "result_percentage", NKcells:Granulocytes) %>% 
#     mutate(cibersort_cell_type = str_replace_all(cibersort_cell_type, "\\.", " ")) %>% 
#     left_join(cs_translation_df) %>% 
#     .[complete.cases(.),] %>% 
#     split(.$actual_cell_type) %>% 
#     map(group_by, actual_cell_type, translated_cibersort_cell_type) %>% 
#     map(dplyr::summarise, predicted_fraction = sum(result_percentage)) %>% 
#     bind_rows %>% 
#     ungroup %>% 
#     mutate(annotated_fraction = ifelse(actual_cell_type == translated_cibersort_cell_type, 1.0, 0.0)) 
# 
# mcp_annotated_df <- mcp_translation_df %>% 
#     .[complete.cases(.),] %>% 
#     inset("annotated_fraction", value = 1)
# 
# mcp_result_df <- "mcp_results.tsv" %>%
#     str_c(tmp_dir, .) %>%
#     read_tsv %>% 
#     df_to_matrix("cell_type") %>% 
#     apply(1, function(row) row / max(row)) %>% 
#     t %>% 
#     matrix_to_df("mcpcounter_cell_type") %>% 
#     gather("translated_mpcounter_cell_type", "nomalized_predicted_score", NKcells:Granulocytes) %>% 
#     left_join(mcp_annotated_df) %>% 
#     mutate(annotated_fraction = ifelse(is.na(annotated_fraction), 0, 1)) 
# 
# 
# 
# 
# png('GSE62408_mcpcounter_facet_scatterplot.png', height = 1000)
# ggplot(mcp_result_df, aes(x = mcpcounter_cell_type, y = nomalized_predicted_score)) +
#     geom_point() +
#     ylab("Nomalized predicted score") +
#     xlab("MCPCounter cell type") +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90, size = 12)) +
#     theme(axis.text.y = element_text(size = 12)) +
#     facet_grid(translated_mpcounter_cell_type ~ .) +
#     theme(strip.text.y = element_text(size = 16)) +
#     ggtitle("MCPCounter GSE62408")
# dev.off()
# 
# png('GSE62408_cibersort_facet_scatterplot.png', height = 1000)
# ggplot(cs_result_df, aes(x = translated_cibersort_cell_type, y = predicted_fraction)) +
#     geom_point() +
#     facet_grid(actual_cell_type ~ .) +
#     ylab("Predicted fraction") +
#     xlab("Translated Cibersort cell type") +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90, size = 12)) +
#     theme(axis.text.y = element_text(size = 12)) +
#     theme(strip.text.y = element_text(size = 16)) +
#     ggtitle("Cibersort GSE62408")
# dev.off()
# 
# 
# 
# png('GSE62408_cibersort_scatterplot.png')
# ggplot(cs_result_df, aes(y = predicted_fraction, x = as.factor(annotated_fraction))) +
#     geom_boxplot(outlier.shape = NA) +
#     geom_jitter(
#         aes(color = actual_cell_type, 
#             shape = translated_cibersort_cell_type,
#             size = as.factor(annotated_fraction))) +
#     scale_shape_manual(values = 0:6) +
#     ylab("Cibersort predicted fraction") +
#     xlab("Annotated raction") +
#     theme_bw() +
#     ggtitle("Cibersort GSE62408")
# dev.off()
# 
# 
# png('GSE62408_mcpcounter_scatterplot.png')
# ggplot(mcp_result_df, aes(y = nomalized_predicted_score, x = as.factor(annotated_fraction))) +
#     geom_boxplot(outlier.shape = NA) +
#     geom_jitter(
#         aes(color = mcpcounter_cell_type, 
#             shape = translated_mpcounter_cell_type,
#             size = as.factor(annotated_fraction))) +
#     scale_shape_manual(values = 0:6) +
#     ylab("Cibersort predicted fraction") +
#     xlab("Annotated fraction") +
#     theme_bw() +
#     ggtitle("MCPCounter GSE62408")
# dev.off()
