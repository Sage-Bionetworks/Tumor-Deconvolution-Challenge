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

fine <- "syn21752551" %>% 
    synapse_file_to_tbl(delim = ",") %>% 
    dplyr::mutate("cell.type" = factor(.data$cell.type, levels = fine_cell_types)) %>% 
    dplyr::group_split(dataset.name) %>% 
    purrr::map(tidyr::complete, dataset.name, sample.id, cell.type) %>% 
    dplyr::bind_rows() %>% 
    readr::write_csv("in-silico-val-fine.csv")

synapser::File("in-silico-val-fine.csv", parent = "syn22361008") %>% 
    synapser::synStore()


"syn21752552" %>% 
    synapse_file_to_tbl(delim = ",") %>% 
    dplyr::rename("cell.type" = "sample") %>% 
    dplyr::group_split(dataset.name) %>% 
    purrr::map(tidyr::complete, dataset.name, sample.id, cell.type) %>% 
    dplyr::bind_rows() %>% 
    readr::write_csv("in-silico-val-coarse.csv")

synapser::File("in-silico-val-coarse.csv", parent = "syn22361008") %>% 
    synapser::synStore()
