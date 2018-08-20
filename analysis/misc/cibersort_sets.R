library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE65135/"


setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()


create_scatter_plot <- function(df, title){
    p <- df %>% 
        ggplot(aes(x = ground_truth_fraction, y = cibersort_predicted_fraction)) +
        geom_point(size = 2, aes(shape = cell_type, color = sample )) +
        geom_smooth(method = 'lm') +
        geom_abline() +
        theme_bw() +
        scale_shape_manual(values = 0:9) +
        ggtitle(df$cell_type[[1]])
    print(p)
}


# ---------------------------------------
gt1 <- create_df_from_synapse_id("syn15664977") %>% 
    gather(key = "cell_type", value = "ground_truth_fraction", -sample) %>% 
    mutate(ground_truth_fraction = ground_truth_fraction / 100)

cs1 <- create_df_from_synapse_id("syn15665004") %>% 
    select(-c(`P-value`, Correlation, RMSE)) %>%
    set_colnames(., str_replace_all(colnames(.), " ", "_")) %>% 
    select(-c(`T_cells_regulatory_(Tregs)`)) %>% 
    mutate(NK_cells = NK_cells_resting + NK_cells_activated) %>% 
    select(-c(NK_cells_resting,NK_cells_activated)) %>% 
    dplyr::rename(Gamma_delta_T_cells = T_cells_gamma_delta) %>% 
    dplyr::rename(Naive_B_cells = B_cells_naive) %>% 
    dplyr::rename(Memory_B_cells = B_cells_memory) %>% 
    dplyr::rename(CD8_T_cells = T_cells_CD8) %>% 
    dplyr::rename(Activated_memory_CD4_T_cells = T_cells_CD4_memory_activated) %>% 
    dplyr::rename(Naive_CD4_T_cells = T_cells_CD4_naive) %>% 
    dplyr::rename(Resting_memory_CD4_T_cells = T_cells_CD4_memory_resting) %>% 
    select(-c(Plasma_cells, Macrophages_M0, Macrophages_M1, Macrophages_M2, 
              Dendritic_cells_resting, Dendritic_cells_activated, 
              Mast_cells_resting, Mast_cells_activated, Eosinophils,
              Neutrophils, T_cells_follicular_helper)) %>% 
    gather(key = "cell_type", value = "cibersort_predicted_fraction", -sample)


combined_df1 <- inner_join(gt1, cs1)

cor_df1 <- combined_df1 %>% 
    select(-sample) %>% 
    group_by(cell_type) %>%
    dplyr::summarise(c = cor(ground_truth_fraction, cibersort_predicted_fraction)) 

# combined_df1 %>%
#     split(.$cell_type) %>% 
#     walk(create_scatter_plot)

combined_df1 %>% 
    ggplot(aes(x = ground_truth_fraction, y = cibersort_predicted_fraction)) +
    geom_point(size = 2, aes(color = sample, shape = cell_type)) +
    geom_smooth(method = 'lm') +
    geom_abline() +
    theme_bw() +
    scale_shape_manual(values = 0:9) +
    ggtitle("GSE65133")




# ---------------------------------------
gt2 <- create_df_from_synapse_id("syn15664984") %>% 
    gather(key = "cell_type", value = "ground_truth_fraction", -sample) %>% 
    mutate(ground_truth_fraction = ground_truth_fraction / 100)

cs2 <- create_df_from_synapse_id("syn15665021") %>% 
    select(-c(`P-value`, Correlation, RMSE)) %>%
    set_colnames(., str_replace_all(colnames(.), " ", "_")) %>% 
    dplyr::rename(Tregs = `T_cells_regulatory_(Tregs)`) %>% 
    select(sample, Tregs) %>% 
    gather(key = "cell_type", value = "cibersort_predicted_fraction", -sample)


combined_df2 <- inner_join(gt2, cs2)

cor_df2 <- combined_df2 %>% 
    select(-sample) %>% 
    group_by(cell_type) %>%
    dplyr::summarise(c = cor(ground_truth_fraction, cibersort_predicted_fraction)) 

combined_df2 %>% 
    ggplot(aes(x = ground_truth_fraction, y = cibersort_predicted_fraction)) +
    geom_point(size = 2, aes(shape = cell_type, color = sample )) +
    geom_smooth(method = 'lm') +
    geom_abline() +
    theme_bw() +
    scale_shape_manual(values = 0:9) +
    ggtitle("GSE65134")



    



# ---------------------------------------
gt3 <- create_df_from_synapse_id("syn15664994") %>% 
    gather(key = "cell_type", value = "ground_truth_fraction", -sample) %>% 
    mutate(ground_truth_fraction = ground_truth_fraction / 100)

cs3 <- create_df_from_synapse_id("syn15665023") %>% 
    select(-c(`P-value`, Correlation, RMSE)) %>%
    set_colnames(., str_replace_all(colnames(.), " ", "_")) %>% 
    mutate(B_cells = B_cells_naive + B_cells_naive) %>% 
    mutate(CD4_T_cells = T_cells_CD4_naive + T_cells_CD4_memory_resting + T_cells_CD4_memory_activated) %>% 
    mutate(CD8_T_cells = T_cells_CD8) %>% 
    select(sample, B_cells, CD4_T_cells, CD8_T_cells) %>% 
    gather(key = "cell_type", value = "cibersort_predicted_fraction", -sample)


combined_df3 <- inner_join(gt3, cs3)

cor_df3 <- combined_df3 %>% 
    select(-sample) %>% 
    group_by(cell_type) %>%
    dplyr::summarise(c = cor(ground_truth_fraction, cibersort_predicted_fraction)) 

combined_df3 %>% 
    ggplot(aes(x = ground_truth_fraction, y = cibersort_predicted_fraction)) +
    geom_point(size = 2, aes(shape = cell_type, color = sample )) +
    geom_smooth(method = 'lm') +
    geom_abline() +
    theme_bw() +
    scale_shape_manual(values = 0:9) +
    ggtitle("GSE65135")