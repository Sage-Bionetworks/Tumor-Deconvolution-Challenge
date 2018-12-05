
## @knitr functions

group_cell_types <- function(df, cell_groupings){
    if(!is.null(cell_groupings)){
        for(cell_grouping in cell_groupings){
            df <- group_by_cell_type(df, cell_grouping)
        }
    } 
    return(df)
}


group_by_cell_type <- function(df, cell_grouping){
    if(!(all(unlist(cell_grouping$old_cols) %in% colnames(df)))) {
        stop(paste0("Could not find ", paste(unlist(cell_grouping$old_cols), collapse=", "), " in:\n",
                    "columns ", paste(colnames(df), collapse=", "), "\n"))
    }
    df <- df %>%
        select(unlist(cell_grouping$old_cols)) %>%
        rowSums %>%
        inset(df, cell_grouping$new_col, value = .)
}

dowload_and_format_mcpcounter_df <- function(synapse_id){
    synapse_id %>% 
        download_from_synapse %>% 
        read.table %>% 
        t %>% 
        .[,colSums(.) > 0] %>% 
        matrix_to_df("sample") %>% 
        set_colnames(str_replace_all(colnames(.), "\\.", "_")) %>% 
        set_colnames(str_replace_all(colnames(.), "[:space:]", "_")) 
}

dowload_and_format_cibersort_df <- function(synapse_id){
    config$synapse_ids$cibersort_results %>% 
        create_df_from_synapse_id %>% 
        select(-c(`P-value`, `RMSE`, Correlation)) %>% 
        df_to_matrix("sample") %>% 
        .[,colSums(.) > 0] %>% 
        matrix_to_df("sample") %>% 
        set_colnames(str_replace_all(colnames(.), "\\.", "_")) %>% 
        set_colnames(str_replace_all(colnames(.), "[:space:]", "_")) 
}


## @knitr cibersort_vs_ground_truth

create_cs_scatter_plot <- function(type, plot_df){
    obj <- cor.test(
        plot_df$predicted_fraction, 
        plot_df$mean_fraction)
    p <- obj$p.value %>% 
        round(4)
    r <- obj$estimate %>% 
        round(4)
    p <- plot_df %>% 
        ggplot(aes(x = predicted_fraction, y = mean_fraction)) +
        geom_point(size = 4, aes(color = sample, shape = cell_type)) +
        geom_smooth(method = 'lm') +
        geom_abline() +
        geom_errorbar(aes(ymin = mean_fraction - sd_fraction,
                          ymax = mean_fraction + sd_fraction),
                      width = sd(plot_df$predicted_fraction) / 8) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, size = 12)) +
        theme(axis.text.y = element_text(size = 12)) +
        theme(strip.text.y = element_text(size = 10, angle = 0)) +
        ggtitle(str_c(type, ", ground truth vs Cibersort predictions, R=", r, " P=", p)) + 
        ylab("Ground truth fraction") +
        xlab("Cibersort predicted fraction")
    print(p)
}

make_cibersort_vs_ground_truth_plots <- function(config){
    
    results_df <- config$synapse_ids$cibersort_results%>%
        dowload_and_format_cibersort_df %>%  
        group_cell_types(config$cs_gt_groups) %>% 
        select(c("sample", config$cibersort_common_groups)) %>% 
        gather("cell_type", "predicted_fraction", -"sample")
    
    
    ground_truth_df <- config$synapse_ids$ground_truth %>%
        create_df_from_synapse_id %>%
        set_colnames(str_replace_all(colnames(.), "\\.", "_")) %>% 
        set_colnames(str_replace_all(colnames(.), "[:space:]", "_")) %>% 
        group_cell_types(config$gt_cs_groups) %>% 
        select(c("sample", config$cibersort_common_groups)) %>% 
        gather("cell_type", "fraction", -sample) %>%
        mutate(fraction = fraction / 100) %>% 
        .[complete.cases(.),] %>% 
        group_by(sample, cell_type) %>% 
        dplyr::summarise(sd_fraction = sd(fraction), mean_fraction = mean(fraction))
    
    plot_df <-
        inner_join(results_df, ground_truth_df)
    
    ## create_cs_scatter_plot("All_cells", plot_df)
    cell_types <- plot_df %>% 
        use_series(cell_type) %>% 
        unique %>% 
        sort
    plot_dfs <- plot_df %>% 
        split(.$cell_type)
    walk2(cell_types, plot_dfs, create_cs_scatter_plot)
    
}

make_cibersort_vs_ground_truth_plots(config)




## @knitr mcpcounter_vs_ground_truth

create_mcp_scatter_plot <- function(type, plot_df){
    obj <- cor.test(
        plot_df$score, 
        plot_df$mean_fraction)
    p <- obj$p.value %>% 
        round(4)
    r <- obj$estimate %>% 
        round(4)
    p <- plot_df %>% 
        ggplot(aes(x = score, y = mean_fraction)) +
        geom_point(size = 4, aes(color = sample, shape = cell_type)) +
        geom_smooth(method = 'lm') +
        geom_errorbar(aes(ymin = mean_fraction - sd_fraction,
                          ymax = mean_fraction + sd_fraction),
                      width = sd(plot_df$score) / 8) +
        theme(axis.text.x = element_text(angle = 90, size = 12)) +
        theme(axis.text.y = element_text(size = 12)) +
        theme(strip.text.y = element_text(size = 10, angle = 0)) +
        ggtitle(str_c(type, ", ground truth vs MCPcounter scores, R=", r, " P=", p)) + 
        ylab("Ground truth fraction") +
        xlab("MCPcounter score")
    print(p)
}


make_mcpcounter_vs_ground_truth_plots <- function(config){
    
    results_df <- config$synapse_ids$mcpcounter_results %>%
        dowload_and_format_mcpcounter_df %>% 
        group_cell_types(config$mcp_gt_groups) %>% 
        select(c("sample", config$mcpcounter_common_groups)) %>% 
        gather("cell_type", "score", -"sample")
    
    ground_truth_df <- config$synapse_ids$ground_truth %>%
        create_df_from_synapse_id %>%
        set_colnames(str_replace_all(colnames(.), "\\.", "_")) %>% 
        set_colnames(str_replace_all(colnames(.), "[:space:]", "_")) %>% 
        group_cell_types(config$gt_mcp_groups) %>% 
        select(c("sample", config$mcpcounter_common_groups)) %>% 
        gather("cell_type", "fraction", -sample) %>%
        mutate(fraction = fraction / 100) %>% 
        .[complete.cases(.),] %>% 
        group_by(sample, cell_type) %>% 
        dplyr::summarise(sd_fraction = sd(fraction), mean_fraction = mean(fraction))
    
    
    plot_df <- inner_join(results_df, ground_truth_df)
    cell_types <- plot_df %>% 
        use_series(cell_type) %>% 
        unique %>% 
        sort
    plot_dfs <- plot_df %>% 
        split(.$cell_type)
    walk2(cell_types, plot_dfs, create_mcp_scatter_plot)
}

make_mcpcounter_vs_ground_truth_plots(config)


## @knitr cibersort gene heatmaps

create_cibersort_gene_heatmaps <- function(config){
    
    
    heatmap_col_df <- config$synapse_ids$annotations %>% 
        create_df_from_synapse_id %>% 
        arrange(sample) %>% 
        data.frame %>% 
        column_to_rownames("sample")
    
    gene_df <- config$synapse_ids$cibersort_genes %>% 
        create_df_from_synapse_id %>% 
        filter(Method == "cibersort") 
    
    genes <- gene_df %>% 
        use_series("Hugo") %>% 
        unique %>% 
        sort
    
    zscore_m <- config$synapse_ids$log_tpm_expression %>% 
        create_df_from_synapse_id %>% 
        df_to_matrix("Hugo") %>% 
        .[rowSums(.) > 0,] %>% 
        quantile_normalize_matrix %>% 
        zscore_matrix %>% 
        .[rownames(.) %in% genes,] %>% 
        .[complete.cases(.),]
    
    pheatmap(
        zscore_m,
        main = "Cibersort genes",
        annotation_col = heatmap_col_df,
        scale = "none",
        fontsize = 15,
        fontsize_row = 5)

}

create_cibersort_gene_heatmaps(config)


## @knitr mcpcounter gene heatmaps

create_mcpcounter_gene_heatmaps <- function(
    annotations, mcpcounter_genes, log_tpm_expression){
    
    mcp_heatmap_col_df <- annotations %>% 
        create_df_from_synapse_id %>% 
        arrange(sample) %>% 
        data.frame %>% 
        column_to_rownames("sample")
    
    mcp_gene_df <- mcpcounter_genes %>% 
        create_df_from_synapse_id %>% 
        filter(Method == "mcpcounter") 
    
    mcp_genes <- mcp_gene_df %>% 
        use_series("Hugo") %>% 
        unique %>% 
        sort
    
    mcp_zscore_matrix <- log_tpm_expression %>% 
        create_df_from_synapse_id %>% 
        df_to_matrix("Hugo") %>% 
        .[rowSums(.) > 0,] %>% 
        quantile_normalize_matrix %>% 
        zscore_matrix %>% 
        .[rownames(.) %in% mcp_genes,] %>% 
        .[complete.cases(.),]
    
    mcp_heatmap_row_df <- mcp_gene_df %>% 
        filter(Method == "mcpcounter") %>% 
        filter(Hugo %in% rownames(mcp_zscore_matrix)) %>% 
        select(-Method) %>% 
        arrange(cell_type) %>% 
        data.frame %>% 
        column_to_rownames("Hugo") %>% 
        set_names("Cell Type")
    
    mcp_zscore_matrix <-  mcp_zscore_matrix[rownames(mcp_heatmap_row_df),]
    
    pheatmap(
        mcp_zscore_matrix,
        main = "MCPCounter genes",
        annotation_row = mcp_heatmap_row_df,
        annotation_col = mcp_heatmap_col_df,
        cluster_rows = F,
        scale = "none")
    
    pheatmap(
        mcp_zscore_matrix,
        main = "MCPCounter genes",
        annotation_row = mcp_heatmap_row_df,
        annotation_col = mcp_heatmap_col_df,
        scale = "none")
}

create_mcpcounter_gene_heatmaps(
    config$synapse_ids$annotations,
    config$synapse_ids$mcpcounter_genes,
    config$synapse_ids$log_tpm_expression
)


## @knitr cibersort_results

create_cibersort_scatterplots <- function(annotations, cibersort_results){
    
    anno_df <- annotations %>% 
        create_df_from_synapse_id %>% 
        arrange(sample) 
    
    cs_result_df <- cibersort_results %>%
        dowload_and_format_cibersort_df %>% 
        gather("cibersort_cell_type", "predicted_fraction", -sample) %>% 
        inner_join(anno_df, by = c("sample"))
    
    cs_plot <- ggplot(cs_result_df, aes(x = cibersort_cell_type, y = predicted_fraction)) +
        geom_point() +
        ylab("Predicted fraction") +
        xlab("Cibersort cell type") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, size = 12)) +
        theme(axis.text.y = element_text(size = 12)) +
        theme(strip.text.y = element_text(size = 10, angle = 0)) +
        ggtitle("Cibersort Results")
    
    if(!is.null(cs_result_df$cell_type)) cs_plot <- cs_plot + facet_grid(cell_type ~ .)
    print(cs_plot)
    
}


create_cibersort_scatterplots(
    config$synapse_ids$annotations,
    config$synapse_ids$cibersort_results
)

## @knitr mcpcounter_results

create_mcpcounter_scatterplots <- function(annotations, mcpcounter_results){
    
    anno_df <- annotations %>% 
        create_df_from_synapse_id %>% 
        arrange(sample) 
    
    mcp_result_df <- mcpcounter_results %>%
        dowload_and_format_mcpcounter_df %>% 
        gather("mcpcounter_cell_type", "predicted_score", -sample) %>% 
        inner_join(anno_df, by = c("sample")) 
    
    mcp_plot <- ggplot(mcp_result_df, aes(x = mcpcounter_cell_type, y = predicted_score)) +
        geom_point() +
        ylab("Predicted score") +
        xlab("MCPCounter cell type") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, size = 12)) +
        theme(axis.text.y = element_text(size = 12)) +
        theme(strip.text.y = element_text(size = 10, angle = 0)) +
        ggtitle("MCPCounter results")
    
    if(!is.null(mcp_result_df$cell_type)) mcp_plot <- mcp_plot + facet_grid(cell_type ~ .)
    print(mcp_plot)
}

create_mcpcounter_scatterplots(
    config$synapse_ids$annotations,
    config$synapse_ids$mcpcounter_results
)


## @knitr pca_plots

create_pca_plot <- function(config){
    
    anno_df <- config$synapse_ids$annotations %>% 
        create_df_from_synapse_id %>% 
        arrange(sample) 
    
    pca_matrix <- config$synapse_ids$log_tpm_expression %>% 
        create_df_from_synapse_id %>% 
        df_to_matrix("Hugo") %>% 
        .[rowSums(.) > 0,] %>% 
        t
    
    if(is.null(config$pca_plot_aes$size)) size = 4
    else{size = config$pca_plot_aes$size}
    
    p <- autoplot(
        prcomp(pca_matrix), 
        data = anno_df, 
        shape = config$pca_plot_aes$shape, 
        size = size,
        colour = config$pca_plot_aes$color,
        main = "PC 1 vs 2") +
        scale_shape_manual(values = 1:19) +
        theme_bw()
    print(p)
}

create_pca_plot(config)


## @knitr cibersort_gsea

create_cibersort_gsea_plot <- function(config){
    
    anno_df <- config$synapse_ids$annotations %>% 
        create_df_from_synapse_id %>% 
        arrange(sample) 
    
    cs_genes <- config$synapse_ids$cibersort_genes %>% 
        create_df_from_synapse_id %>%  
        filter(Method == "cibersort") %>%
        split(.$cell_type) %>%
        map(use_series, Hugo)
    
    cs_ssgsea_df <- config$synapse_ids$log_tpm_expression %>% 
        create_df_from_synapse_id %>% 
        df_to_matrix("Hugo") %>% 
        .[rowSums(.) > 0,] %>% 
        gsva(cs_genes, method = "ssgsea", verbose = F) %>%
        matrix_to_df("CS_cell_type") %>%
        gather(key = "sample", value = "enrichment" , -CS_cell_type) %>%
        left_join(anno_df) 
    
    plot <- ggplot(cs_ssgsea_df, aes(x = CS_cell_type, y = enrichment)) +
        geom_point() +
        ylab("GSEA enrichment score") +
        xlab("Cibersort cell type") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, size = 12)) +
        theme(axis.text.y = element_text(size = 12)) +
        theme(strip.text.y = element_text(size = 10, angle = 0)) +
        ggtitle("GSEA with Cibersort genes")
    
    if(!is.null(anno_df$cell_type)) plot <- plot + facet_grid(cell_type ~ .)
    print(plot)
}

create_cibersort_gsea_plot(config)

## @knitr mcpcounter_gsea
create_mcpcounter_gsea_plot <- function(config){
    anno_df <- config$synapse_ids$annotations %>% 
        create_df_from_synapse_id %>% 
        arrange(sample) 
    
    mcp_genes <- config$synapse_ids$mcpcounter_genes %>% 
        create_df_from_synapse_id %>%  
        filter(Method == "mcpcounter") %>%
        split(.$cell_type) %>%
        map(use_series, Hugo)
    
    mcp_ssgsea_df <- config$synapse_ids$log_tpm_expression %>% 
        create_df_from_synapse_id %>% 
        df_to_matrix("Hugo") %>% 
        .[rowSums(.) > 0,] %>% 
        gsva(mcp_genes, method = "ssgsea", verbose = F) %>%
        matrix_to_df("MCP_cell_type") %>%
        gather(key = "sample", value = "enrichment" , -MCP_cell_type) %>%
        left_join(anno_df) 
    
    plot <- ggplot(mcp_ssgsea_df, aes(x = MCP_cell_type, y = enrichment)) +
        geom_point() +
        ylab("GSEA enrichment score") +
        xlab("MCPcounter cell type") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, size = 12)) +
        theme(axis.text.y = element_text(size = 12)) +
        theme(strip.text.y = element_text(size = 10, angle = 0)) +
        ggtitle("GSEA with MCPcounter genes")
    
    if(!is.null(anno_df$cell_type)) plot <- plot + facet_grid(cell_type ~ .)
    print(plot)
}

create_mcpcounter_gsea_plot(config)

