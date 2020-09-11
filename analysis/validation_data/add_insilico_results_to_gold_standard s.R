library(magrittr)

synapser::synLogin()
devtools::source_url("https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R")

parent_dir <- "syn21820011"

file_to_synapse <- function(file, id){
    file %>% 
        synapser::File(parent = id) %>% 
        synapser::synStore()
    file.remove(file)
}

cells <- list(
    "coarse" = c(
        "B.cells",
        "CD4.T.cells",
        "CD8.T.cells",
        "NK.cells",
        "neutrophils",
        "monocytic.lineage",
        "fibroblasts",
        "endothelial.cells"
    ),
    "fine" = c(
        "memory.B.cells",
        "naive.B.cells",
        "memory.CD4.T.cells",
        "naive.CD4.T.cells",
        "regulatory.T.cells",
        "memory.CD8.T.cells",
        "naive.CD8.T.cells",
        "NK.cells",
        "neutrophils",
        "monocytes",
        "myeloid.dendritic.cells",
        "macrophages",
        "fibroblasts",
        "endothelial.cells"
    )
)

ids <- list(
    "coarse_invitro" = "syn21820375",
    "fine_invitro"   = "syn21820376",
    "coarse_insilico" = "syn22013527",
    "fine_ininsilico" = "syn22013526"
)

file_names <- list(
    "coarse" = "coarse.csv",
    "fine"   = "fine.csv"
)

coarse_insilico <- "syn22013527" %>% 
    synapse_file_to_tbl(delim = ",") %>% 
    dplyr::mutate(cell.type = factor(cell.type, levels = cells$coarse)) %>% 
    dplyr::group_by(dataset.name) %>%
    dplyr::group_split() %>%
    purrr::map(tidyr::complete, dataset.name, sample.id, cell.type) %>%
    dplyr::bind_rows()

"syn21820375" %>% 
    synapse_file_to_tbl(delim = ",") %>% 
    dplyr::filter(.data$dataset.name %in% c("DS1", "DS2", "DS3", "DS4")) %>% 
    dplyr::bind_rows(coarse_insilico) %>% 
    readr::write_csv("coarse.csv")
    
file_to_synapse("coarse.csv", parent_dir)
    

fine_insilico <- "syn22013526" %>% 
    synapse_file_to_tbl(delim = ",") %>% 
    dplyr::mutate(cell.type = factor(cell.type, levels = cells$fine)) %>% 
    dplyr::group_by(dataset.name) %>%
    dplyr::group_split() %>%
    purrr::map(tidyr::complete, dataset.name, sample.id, cell.type) %>%
    dplyr::bind_rows()

"syn21820376" %>% 
    synapse_file_to_tbl(delim = ",") %>% 
    dplyr::filter(.data$dataset.name %in% c("DS1", "DS2", "DS3", "DS4")) %>% 
    dplyr::bind_rows(fine_insilico) %>% 
    readr::write_csv("fine.csv")

file_to_synapse("fine.csv", parent_dir)
