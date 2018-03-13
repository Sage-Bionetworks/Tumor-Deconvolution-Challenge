library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(pheatmap)
library(preprocessCore)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/E-MTAB-2319/"

count_id      <- "syn11958709"
annotation_id <- "syn11958711"
hugo_id       <- "syn11536071"
count_id      <- "syn11915424"
genes_id      <- "syn11918430"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores())

gene_df <-  create_df_from_synapse_id(genes_id)
hugo_df <-  create_df_from_synapse_id(hugo_id)
count_df <- count_id %>% 
    create_df_from_synapse_id(unzip = T) %>% 
    select(-c(Chr, Start, End, Strand, Length))
    # inner_join(hugo_df, by = c("Geneid" = "ensembl_gene_id")) %>% 
    
    # select(hgnc_symbol, everything()) %>% 
    # filter(hgnc_symbol != "")

# group_by(hgnc_symbol) %>% 
#     summarise_all(.funs = sum) %>% 
#     ungroup %>% 

cpm_matrix <- count_df %>% 
    df_to_matrix("Geneid") %>% 
    .[rowSums(.) > 0,] %>%  
    add(1) %>% 
    apply(2, calculate_cpm) %>% 
    matrix_to_df("ensembl_gene_id") %>% 
    inner_join(hugo_df) %>% 
    select(hgnc_symbol, everything()) %>% 
    filter(hgnc_symbol != "") %>% 
    select(-ensembl_gene_id) %>% 
    group_by(hgnc_symbol) %>% 
    summarise_all(.funs = sum) %>% 
    ungroup %>% 
    df_to_matrix("hgnc_symbol")


%>% 
    log10 %>% 
    quantile_normalize_matrix %>% 
    zscore_matrix

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

png('E-MTAB-2319_mcpcounter_genes_heatmap.png', width = 4000, height = 4000)
pheatmap(
    mcp_quant_matrix,
    main = "MCPCounter E-MTAB-2319",
    annotation_row = mcp_heatmap_row_df,
    cluster_rows = F,
    scale = "none")
dev.off()

png('E-MTAB-2319_mcpcounter_genes_rows_clustered_heatmap.png', width = 4000, height = 4000)
pheatmap(
    mcp_quant_matrix,
    main = "MCPCounter E-MTAB-2319",
    annotation_row = mcp_heatmap_row_df,
    scale = "none")
dev.off()

png('E-MTAB-2319_cibersort_genes_heatmap.png', width = 4000, height = 4000)
pheatmap(
    cs_zscore_matrix,
    main = "Cibersort E-MTAB-2319",
    scale = "none")
dev.off()

# -----------------------------------------------------------------------------

cs_result_df <- "cs_results.tsv" %>%
    str_c(tmp_dir, .) %>%
    read_tsv %>%
    transpose_df("cell_type", "cibersort_cell_type") %>% 
    gather("actual_cell_type", "result_percentage", NKcells:Granulocytes) %>% 
    mutate(cibersort_cell_type = str_replace_all(cibersort_cell_type, "\\.", " ")) %>% 
    left_join(cs_translation_df) %>% 
    .[complete.cases(.),] %>% 
    split(.$actual_cell_type) %>% 
    map(group_by, actual_cell_type, translated_cibersort_cell_type) %>% 
    map(dplyr::summarise, predicted_fraction = sum(result_percentage)) %>% 
    bind_rows %>% 
    ungroup %>% 
    mutate(annotated_fraction = ifelse(actual_cell_type == translated_cibersort_cell_type, 1.0, 0.0)) 

mcp_annotated_df <- mcp_translation_df %>% 
    .[complete.cases(.),] %>% 
    inset("annotated_fraction", value = 1)

mcp_result_df <- "mcp_results.tsv" %>%
    str_c(tmp_dir, .) %>%
    read_tsv %>% 
    df_to_matrix("cell_type") %>% 
    apply(1, function(row) row / max(row)) %>% 
    t %>% 
    matrix_to_df("mcpcounter_cell_type") %>% 
    gather("translated_mpcounter_cell_type", "nomalized_predicted_score", NKcells:Granulocytes) %>% 
    left_join(mcp_annotated_df) %>% 
    mutate(annotated_fraction = ifelse(is.na(annotated_fraction), 0, 1)) 




png('E-MTAB-2319_mcpcounter_facet_scatterplot.png', height = 1000)
ggplot(mcp_result_df, aes(x = mcpcounter_cell_type, y = nomalized_predicted_score)) +
    geom_point() +
    ylab("Nomalized predicted score") +
    xlab("MCPCounter cell type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    facet_grid(translated_mpcounter_cell_type ~ .) +
    theme(strip.text.y = element_text(size = 16)) +
    ggtitle("MCPCounter E-MTAB-2319")
dev.off()

png('E-MTAB-2319_cibersort_facet_scatterplot.png', height = 1000)
ggplot(cs_result_df, aes(x = translated_cibersort_cell_type, y = predicted_fraction)) +
    geom_point() +
    facet_grid(actual_cell_type ~ .) +
    ylab("Predicted fraction") +
    xlab("Translated Cibersort cell type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(strip.text.y = element_text(size = 16)) +
    ggtitle("Cibersort E-MTAB-2319")
dev.off()



png('E-MTAB-2319_cibersort_scatterplot.png')
ggplot(cs_result_df, aes(y = predicted_fraction, x = as.factor(annotated_fraction))) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        aes(color = actual_cell_type, 
            shape = translated_cibersort_cell_type,
            size = as.factor(annotated_fraction))) +
    scale_shape_manual(values = 0:6) +
    ylab("Cibersort predicted fraction") +
    xlab("Annotated raction") +
    theme_bw() +
    ggtitle("Cibersort E-MTAB-2319")
dev.off()


png('E-MTAB-2319_mcpcounter_scatterplot.png')
ggplot(mcp_result_df, aes(y = nomalized_predicted_score, x = as.factor(annotated_fraction))) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        aes(color = mcpcounter_cell_type, 
            shape = translated_mpcounter_cell_type,
            size = as.factor(annotated_fraction))) +
    scale_shape_manual(values = 0:6) +
    ylab("Cibersort predicted fraction") +
    xlab("Annotated fraction") +
    theme_bw() +
    ggtitle("MCPCounter E-MTAB-2319")
dev.off()
