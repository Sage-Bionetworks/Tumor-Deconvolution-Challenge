library(tidyverse)
library(magrittr)
library(synapser)
library(MCPcounter)

devtools::source_url("https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R")
source("../../challenge_models/cibersort_coarse/docker_files/CIBERSORT.R")
source("dataset-setup.R")
source("../../scripts/utils.R")

synapser::synLogin()


cibersort_coarse_tbl <- tibble::tribble(
    ~cell.type, ~cibersort.cell.type,
    "B.cells", "B cells naive",
    "B.cells", "B cells memory",
    "CD4.T.cells", "T cells CD4 naive", 
    "CD4.T.cells", "T cells CD4 memory resting", 
    "CD4.T.cells", "T cells CD4 memory activated",
    "CD4.T.cells", "T cells regulatory (Tregs)", 
    "CD4.T.cells", "T cells follicular helper",
    "CD8.T.cells", "T cells CD8",
    "T.cells.gamma.delta", "T cells CD8",    
    "NK.cells", "NK cells resting", 
    "NK.cells", "NK cells activated",
    "neutrophils", "Neutrophils",
    "monocytic.lineage", "Monocytes",
    "monocytic.lineage", "Macrophages M0",
    "monocytic.lineage", "Macrophages M1",
    "monocytic.lineage", "Macrophages M2",
    "monocytic.lineage", "Dendritic cells resting",
    "monocytic.lineage", "Dendritic cells activated",
    "fibroblasts", "non.immune",
    "endothelial.cells", "non.immune"
)

cibersort_fine_tbl <- tibble::tribble(
    ~cell.type, ~cibersort.cell.type,
    "memory.B.cells", "B cells memory",
    "naive.B.cells", "B cells naive",
    "memory.CD4.T.cells", "T cells CD4 memory resting",
    "memory.CD4.T.cells", "T cells CD4 memory activated",
    "naive.CD4.T.cells", "T cells CD4 naive",
    "regulatory.T.cells", "T cells regulatory (Tregs)",
    "memory.CD8.T.cells", "T cells CD8",
    "naive.CD8.T.cells", "T cells CD8",
    "NK.cells", "NK cells resting", 
    "NK.cells", "NK cells activated",
    "neutrophils", "Neutrophils",
    "monocytes", "Monocytes",
    "myeloid.dendritic.cells", "Dendritic cells resting", 
    "myeloid.dendritic.cells", "Dendritic cells activated",
    "macrophages", "Macrophages M0",
    "macrophages", "Macrophages M1",
    "macrophages", "Macrophages M2",
    "fibroblasts", "non.immune",
    "endothelial.cells", "non.immune"
)

mcpcounter_coarse_tbl <- tibble::tribble(
    ~cell.type, ~mcpcounter.cell.type,
    "B.cells", "B lineage",
    "CD8.T.cells", "CD8 T cells",
    "NK.cells", "NK cells",
    "neutrophils", "Neutrophils",
    "monocytic.lineage", "Monocytic lineage",
    "fibroblasts", "Fibroblasts",
    "endothelial.cells", "Endothelial cells"
)

mcpcounter_fine_tbl <- tibble::tribble(
    ~cell.type, ~mcpcounter.cell.type,
    "NK.cells", "NK cells",
    "neutrophils", "Neutrophils",
    "myeloid.dendritic.cells", "Myeloid dendritic cells",
    "fibroblasts", "Fibroblasts",
    "endothelial.cells", "Endothelial cells",
)

coarse_mix_tbl <- coarse_mix_id %>% 
    create_df_from_synapse_id() %>% 
    dplyr::mutate(sample_name = stringr::str_c("s", 1:dplyr::n()))

obfuscated.dataset <- coarse_mix_tbl$dataset.name[[1]]
fine_datatset_name <- stringr::str_c("F", obfuscated.dataset)
coarse_datatset_name <- stringr::str_c("C", obfuscated.dataset)

file_tbl <- create_entity_tbl(dataset_id)

coarse_gt_tbl <- file_tbl %>% 
    dplyr::filter(name == "coarse_gt.csv") %>% 
    dplyr::pull(id) %>% 
    synapse_file_to_tbl(delim = ",") %>% 
    tidyr::drop_na()

fine_gt_tbl <- file_tbl %>% 
    dplyr::filter(name == "fine_gt.csv") %>% 
    dplyr::pull(id) %>% 
    synapse_file_to_tbl(delim = ",") %>% 
    tidyr::drop_na()

coarse_expr_tbl <- file_tbl %>% 
    dplyr::filter(name == stringr::str_c(
        coarse_datatset_name, 
        "-hugo-gene-expr.csv"
    )) %>% 
    dplyr::pull(id) %>% 
    synapse_file_to_tbl(delim = ",")

fine_expr_tbl <- file_tbl %>% 
    dplyr::filter(name == stringr::str_c(
        fine_datatset_name, 
        "-hugo-gene-expr.csv"
    )) %>% 
    dplyr::pull(id) %>% 
    synapse_file_to_tbl(delim = ",")

MCP_coarse_res_tbl <- coarse_expr_tbl %>% 
    data.frame() %>% 
    tibble::column_to_rownames("Gene") %>% 
    as.matrix() %>% 
    MCPcounter.estimate(featuresType = "HUGO_symbols") %>% 
    data.frame() %>% 
    tibble::rownames_to_column("mcpcounter.cell.type") %>% 
    dplyr::as_tibble() %>% 
    tidyr::pivot_longer(
        -mcpcounter.cell.type, 
        names_to = "sample.id",
        values_to = "predicted"
    ) %>% 
    dplyr::inner_join(mcpcounter_coarse_tbl) %>% 
    dplyr::group_by(cell.type, sample.id) %>%
    dplyr::summarise(predicted = mean(predicted)) %>% 
    dplyr::ungroup() %>% 
    dplyr::inner_join(coarse_gt_tbl) %>% 
    dplyr::mutate(model = "mcpcounter_coarse")
    

MCP_fine_res_tbl <- fine_expr_tbl %>% 
    data.frame() %>% 
    tibble::column_to_rownames("Gene") %>% 
    as.matrix() %>% 
    MCPcounter.estimate(featuresType = "HUGO_symbols") %>% 
    data.frame() %>% 
    tibble::rownames_to_column("mcpcounter.cell.type") %>% 
    dplyr::as_tibble() %>% 
    tidyr::pivot_longer(
        -mcpcounter.cell.type, 
        names_to = "sample.id",
        values_to = "predicted"
    ) %>% 
    dplyr::inner_join(mcpcounter_fine_tbl) %>% 
    dplyr::group_by(cell.type, sample.id) %>%
    dplyr::summarise(predicted = mean(predicted)) %>% 
    dplyr::ungroup() %>% 
    dplyr::inner_join(fine_gt_tbl) %>% 
    dplyr::mutate(model = "mcpcounter_fine")


readr::write_tsv(coarse_expr_tbl, "coarse.tsv")
readr::write_tsv(fine_expr_tbl, "fine.tsv")


Cib_coarse_res_tbl <- 
    CIBERSORT(
        "coarse.tsv",
        "../../challenge_models/cibersort_coarse/docker_files/LM22.tsv",
        ##        QN = F
        absolute = TRUE,
        abs_method = "sig.score"
        # absmean = TRUE
    ) %>% 
    data.frame() %>% 
    tibble::rownames_to_column("cibersort.cell.type") %>% 
    dplyr::as_tibble() %>% 
    dplyr::select(-c(P.value, Correlation, RMSE)) %>% 
    tidyr::pivot_longer(
        -cibersort.cell.type, 
        names_to = "sample.id",
        values_to = "predicted"
    ) %>% 
    dplyr::inner_join(cibersort_coarse_tbl) %>% 
    dplyr::group_by(cell.type, sample.id) %>%
    dplyr::summarise(predicted = mean(predicted)) %>% 
    dplyr::ungroup() %>% 
    dplyr::inner_join(coarse_gt_tbl) %>% 
    dplyr::mutate(model = "cibersort_coarse")

Cib_fine_res_tbl <- 
    CIBERSORT(
        "fine.tsv",
        "../../challenge_models/cibersort_coarse/docker_files/LM22.tsv",
##        QN = F
        absolute = TRUE,
        abs_method = "sig.score"
        # absmean = TRUE
    ) %>% 
    data.frame() %>% 
    tibble::rownames_to_column("cibersort.cell.type") %>% 
    dplyr::as_tibble() %>% 
    dplyr::select(-c(P.value, Correlation, RMSE)) %>% 
    tidyr::pivot_longer(
        -cibersort.cell.type, 
        names_to = "sample.id",
        values_to = "predicted"
    ) %>% 
    dplyr::inner_join(cibersort_fine_tbl) %>% 
    dplyr::group_by(cell.type, sample.id) %>%
    dplyr::summarise(predicted = mean(predicted)) %>% 
    dplyr::ungroup() %>% 
    dplyr::inner_join(fine_gt_tbl) %>% 
    dplyr::mutate(model = "cibersort_fine")

result_tbl <- 
    list(
        Cib_coarse_res_tbl, 
        Cib_fine_res_tbl, 
        MCP_fine_res_tbl, 
        MCP_coarse_res_tbl
    ) %>% 
    dplyr::bind_rows() %>% 
    dplyr::group_by(cell.type, model) %>% 
    dplyr::summarise(
        pearson = cor(predicted, measured, method = "pearson"),
        spearman = cor(predicted, measured, method = "spearman")
    ) %>% 
    tidyr::separate(model, into = c("method", "subchallenge"), sep = "_")

upload_tbl_to_synapse <- function(tbl, file_name, id, delim){
    readr::write_delim(tbl, file_name, delim)
    file_entity <- synapser::File(path = file_name, parent = id)
    synapser::synStore(file_entity)
}

write_tbl <- function(tbl, file_name, id, delim){
    readr::write_delim(tbl, file_name, delim)
}

## upload_tbl_to_synapse(result_tbl, "model_correlations.csv", dataset_id, ",")
write_tbl(result_tbl, "model_correlations.csv", dataset_id, ",")


create_fit_plot <- function(title, data){
    p <- data %>%
        ggplot(aes(x = measured, y = predicted)) +
        geom_point() +
        geom_smooth(method = 'lm') +
        ggtitle(title)
    print(p)
}


plot_table <- 
    list(
        Cib_coarse_res_tbl, 
        Cib_fine_res_tbl, 
        MCP_fine_res_tbl, 
        MCP_coarse_res_tbl
    ) %>% 
    dplyr::bind_rows() %>% 
    dplyr::group_by(cell.type, model) %>% 
    dplyr::mutate(
        pearson = cor(predicted, measured, method = "pearson"),
        title = stringr::str_c(model, cell.type, pearson, sep = "; pearson: ")
    ) %>% 
    dplyr::ungroup() %>%
    dplyr::arrange(cell.type) %>%
    dplyr::select(title, predicted, measured) %>% 
    dplyr::group_by(title) %>% 
    tidyr::nest()

pdf("all-fits.pdf", onefile = TRUE)
purrr::pmap(plot_table, create_fit_plot)
d <- dev.off()


