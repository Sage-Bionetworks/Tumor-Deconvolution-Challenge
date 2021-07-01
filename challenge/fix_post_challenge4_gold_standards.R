synapser::synLogin()

devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)

fine_cell_types <- c(
    "myeloid.dendritic.cells",
    "endothelial.cells",
    "fibroblasts",
    "macrophages",
    "memory.CD4.T.cells",
    "memory.CD8.T.cells",
    "monocytes",
    "naive.B.cells",
    "naive.CD4.T.cells",
    "naive.CD8.T.cells",
    "neutrophils",
    "NK.cells",
    "regulatory.T.cells",
    "memory.B.cells"
)

coarse_cell_types <- c(
    "B.cells",
    "CD4.T.cells",
    "CD8.T.cells",
    "NK.cells",
    "neutrophils",
    "monocytic.lineage",
    "fibroblasts",
    "endothelial.cells"
)

fine <- "syn23019512" %>% 
    synapse_file_to_tbl(delim = ",") %>% 
    dplyr::mutate("cell.type" = factor(.data$cell.type, levels = fine_cell_types)) %>% 
    dplyr::group_split(dataset.name) %>% 
    purrr::map(tidyr::complete, dataset.name, sample.id, cell.type) %>% 
    dplyr::bind_rows() %>% 
    readr::write_csv("ball-fine-gold-standard.csv")

synapser::File("ball-fine-gold-standard.csv", parent = "syn23019061") %>% 
    synapser::synStore()

coarse <- "syn23019511" %>% 
    synapse_file_to_tbl(delim = ",") %>% 
    dplyr::mutate("cell.type" = factor(.data$cell.type, levels = coarse_cell_types)) %>% 
    dplyr::group_split(dataset.name) %>% 
    purrr::map(tidyr::complete, dataset.name, sample.id, cell.type) %>% 
    dplyr::bind_rows() %>% 
    readr::write_csv("ball-coarse-gold-standard.csv")

synapser::File("ball-coarse-gold-standard.csv", parent = "syn23019061") %>% 
    synapser::synStore()
