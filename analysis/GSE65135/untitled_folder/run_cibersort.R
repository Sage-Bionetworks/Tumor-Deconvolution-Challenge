source("./CIBERSORT.R")

library(tidyverse)
library(magrittr)

result_matrix <- CIBERSORT("LM22.txt", "expression_affy_from_cel.tsv") 

result_df <- result_matrix %>% 
    as.data.frame %>% 
    rownames_to_column("sample") %>% 
    as_data_frame() %>% 
    select(-c(`P-value`, Correlation, RMSE)) %>%
    set_colnames(., str_replace_all(colnames(.), " ", "_")) %>% 
    mutate(B_cells = B_cells_naive + B_cells_memory + Plasma_cells) %>% 
    mutate(CD4_T_cells = 
               T_cells_CD4_naive + 
               T_cells_CD4_memory_resting + 
               T_cells_CD4_memory_activated + 
               T_cells_follicular_helper + 
               `T_cells_regulatory_(Tregs)`) %>% 
    mutate(CD8_T_cells = T_cells_CD8) %>% 
    select(sample, B_cells, CD4_T_cells, CD8_T_cells) %>% 
    mutate(total = B_cells + CD4_T_cells + CD8_T_cells) %>% 
    gather(key = "cell_type", value = "cibersort_predicted_fraction", -c(sample, total)) %>% 
    mutate(cibersort_predicted_fraction_normalized = cibersort_predicted_fraction / total) %>% 
    mutate(cibersort_predicted_percent = cibersort_predicted_fraction_normalized * 100)

gt_df <- read_tsv("ground_truth.tsv") %>% 
    gather(key = "cell_type", value = "ground_truth_fraction", -sample) %>% 
    rename(ground_truth_percent = ground_truth_fraction)


combined_df <- inner_join(gt_df, result_df)

cell_types <- result_df %>% 
    use_series(cell_type) %>% 
    unique %>% 
    sort 

dfs <- split(combined_df, combined_df$cell_type) 


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
create_scatter_plot(dfs[[1]], "B_cells")
dev.off()

svg("cd4t.svg", width = 4, height = 6)
create_scatter_plot(dfs[[2]], "CD4_T_cells")
dev.off()

svg("cd8t.svg", width = 4, height = 6)
create_scatter_plot(dfs[[3]], "CD8_T_cells")
dev.off()

