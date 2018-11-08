
library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)


home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/Documents/Presentations/Tumor_d/9_12_18/"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores())

# fastq mixing
cs_results_25_id  <- "syn12385636"
cs_results_80_id  <- "syn16801065"
mcp_results_25_id <- "syn12385637"
mcp_results_80_id <- "syn16801066"


cs_result_df_25 <- cs_results_25_id %>%
    create_df_from_synapse_id %>%
    dplyr::select(sample, `T cells CD8`) %>%
    gather("cibersort_cell_type", "predicted_fraction", `T cells CD8`) %>% 
    mutate(reads = "25M") %>% 
    mutate(CD8_fraction = str_match(sample, "CD4_CD8_([:print:]+)_rep")[,2]) %>% 
    mutate(CD8_percent = as.character(100 * as.numeric(CD8_fraction))) %>% 
    mutate(predicted_percent = 100 * as.numeric(predicted_fraction)) 

           
cs_result_df_80 <- cs_results_80_id %>%
    create_df_from_synapse_id %>%
    dplyr::select(sample, `T cells CD8`) %>%
    gather("cibersort_cell_type", "predicted_fraction", `T cells CD8`) %>% 
    mutate(reads = "80M") %>% 
    mutate(sample = str_replace_all(sample, "5e-04", "0.0005")) %>% 
    mutate(CD8_fraction = str_match(sample, "CD4_CD8_([:print:]+)_rep")[,2]) %>% 
    mutate(CD8_percent = as.character(100 * as.numeric(CD8_fraction))) %>% 
    mutate(predicted_percent = 100 * as.numeric(predicted_fraction))

cs_result_df <- 
    bind_rows(cs_result_df_25, cs_result_df_80) %>% 
    group_by(cibersort_cell_type, CD8_percent, reads) %>%
    summarise(stdev = sd(predicted_percent), mean = mean(predicted_percent)) 




mcp_result_25_df <- mcp_results_25_id %>%
    download_from_synapse %>%
    read.table %>%
    t %>%
    matrix_to_df("sample") %>%
    set_colnames(str_replace_all(colnames(.), "\\.", " ")) %>%
    dplyr::select(sample, `CD8 T cells`, `Cytotoxic lymphocytes`) %>%
    gather("mcpcounter_cell_type", "predicted_score", -sample) %>%
    mutate(CD8_fraction = str_match(sample, "CD4_CD8_([:print:]+)_rep")[,2]) %>% 
    mutate(CD8_percent = as.character(100 * as.numeric(CD8_fraction)))

mcp_result_80_df <- mcp_results_80_id %>%
    download_from_synapse %>%
    read.table %>%
    t %>%
    matrix_to_df("sample") %>%
    set_colnames(str_replace_all(colnames(.), "\\.", " ")) %>%
    dplyr::select(sample, `CD8 T cells`, `Cytotoxic lymphocytes`) %>%
    gather("mcpcounter_cell_type", "predicted_score", -sample) %>%
    mutate(sample = str_replace_all(sample, "5e.04", "0.0005")) %>% 
    mutate(CD8_fraction = str_match(sample, "CD4_CD8_([:print:]+)_rep")[,2]) %>% 
    mutate(CD8_percent = as.character(100 * as.numeric(CD8_fraction))) 
    

svg("cibesort_results_25M.svg", width = 7, height = 10)
cs_result_df_25 %>%
    group_by(cibersort_cell_type, CD8_percent) %>%
    summarise(stdev = sd(predicted_percent), mean = mean(predicted_percent)) %>%
    ggplot(aes(x = CD8_percent, y = mean)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    xlab("CD8+ percent") +
    ylab("Cibersort percent") +
    ggtitle("Cibersort CD4+ with CD8+ spike in, 25M reads") +
    geom_errorbar(aes(ymin = mean - stdev,
                      ymax = mean + stdev),
                  width = 1,
                  position = position_dodge(0.05))
dev.off()

svg("cibesort_results_80M.svg", width = 7, height = 10)
cs_result_df_80 %>%
    group_by(cibersort_cell_type, CD8_percent) %>%
    summarise(stdev = sd(predicted_percent), mean = mean(predicted_percent)) %>%
    ggplot(aes(x = CD8_percent, y = mean)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    xlab("CD8+ percent") +
    ylab("Cibersort percent") +
    ggtitle("Cibersort CD4+ with CD8+ spike in, 80M reads") +
    geom_errorbar(aes(ymin = mean - stdev,
                      ymax = mean + stdev),
                  width = 1,
                  position = position_dodge(0.05))
dev.off()

svg("mcpcounter_results_25M.svg", width = 7, height = 10)
mcp_result_25_df %>%
    group_by(mcpcounter_cell_type, CD8_percent) %>%
    summarise(stdev = sd(predicted_score), mean = mean(predicted_score)) %>%
    ggplot(aes(x = CD8_percent, y = mean,  color = mcpcounter_cell_type)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    xlab("CD8+ percent") +
    ylab("MCPcounter scores") +
    ggtitle("MCPCounter CD4+ with CD8+ spike in, 25M reads") +
    geom_errorbar(aes(ymin = mean - stdev,
                      ymax = mean + stdev),
                  width = 1,
                  position = position_dodge(0.05))
dev.off()

svg("mcpcounter_results_80M.svg", width = 7, height = 10)
mcp_result_80_df %>%
    group_by(mcpcounter_cell_type, CD8_percent) %>%
    summarise(stdev = sd(predicted_score), mean = mean(predicted_score)) %>%
    ggplot(aes(x = CD8_percent, y = mean,  color = mcpcounter_cell_type)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    xlab("CD8+ percent") +
    ylab("MCPcounter scores") +
    ggtitle("MCPCounter CD4+ with CD8+ spike in, 80M reads") +
    geom_errorbar(aes(ymin = mean - stdev,
                      ymax = mean + stdev),
                  width = 1,
                  position = position_dodge(0.05))
dev.off()


# cibersort data sets
gt_GSE65135_df <- create_df_from_synapse_id("syn15664994") %>% 
    gather(key = "cell_type", value = "ground_truth_fraction", -sample) %>% 
    mutate(ground_truth_percent = ground_truth_fraction)

cs_GSE65135_df <- create_df_from_synapse_id("syn16784404") %>% 
    select(-c(`P-value`, Correlation, RMSE)) %>%
    set_colnames(., str_replace_all(colnames(.), " ", "_")) %>% 
    mutate(B_cells = B_cells_naive + B_cells_memory + Plasma_cells) %>% 
    mutate(CD4_T_cells = T_cells_CD4_naive + T_cells_CD4_memory_resting + T_cells_CD4_memory_activated + T_cells_follicular_helper + `T_cells_regulatory_(Tregs)`) %>% 
    mutate(CD8_T_cells = T_cells_CD8) %>% 
    select(sample, B_cells, CD4_T_cells, CD8_T_cells) %>% 
    mutate(total = B_cells + CD4_T_cells + CD8_T_cells) %>% 
    gather(key = "cell_type", value = "cibersort_predicted_fraction", -c(sample, total)) %>% 
    mutate(cibersort_predicted_fraction_normalized = cibersort_predicted_fraction / total) %>% 
    mutate(cibersort_predicted_percent = cibersort_predicted_fraction_normalized * 100)


GSE65135_df <- inner_join(gt_GSE65135_df, cs_GSE65135_df)

cell_types <- GSE65135_df %>% 
    use_series(cell_type) %>% 
    unique %>% 
    sort 

GSE65135_dfs <- split(GSE65135_df, GSE65135_df$cell_type) 




create_scatter_plot <- function(df, cell_type){
    obj <- cor.test(df$cibersort_predicted_percent, df$ground_truth_percent)
    p <- obj$p.value %>% 
        round(4)
    r <- obj$estimate %>% 
        round(4)
    title <- str_c(cell_type, " R=", r, " P=", p)
    p <- df %>% 
        ggplot(aes(x = ground_truth_percent, y = cibersort_predicted_percent)) +
        geom_point(size = 2) +
        geom_smooth(method = 'lm') +
        geom_abline() +
        theme_bw() +
        ggtitle(title) +
        ylab("Cibersort predicted percent") +
        xlab("Flow cytometry percent")
    print(p)
}



svg("b_cells.svg", width = 4, height = 6)
create_scatter_plot(GSE65135_dfs[[1]], "B_cells")
dev.off()

svg("cd4t.svg", width = 4, height = 6)
create_scatter_plot(GSE65135_dfs[[2]], "CD4_T_cells")
dev.off()

svg("cd8t.svg", width = 4, height = 6)
create_scatter_plot(GSE65135_dfs[[3]], "CD8_T_cells")
dev.off()
