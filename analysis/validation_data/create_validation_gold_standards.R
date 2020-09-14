library(magrittr)

synapser::synLogin()
devtools::source_url("https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R")

parent_dir <- "syn21820011"

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
    "coarse" = "syn21590364",
    "fine"   = "syn21590365"
)

file_names <- list(
    "coarse" = "coarse.csv",
    "fine"   = "fine.csv"
)

format_gold_standard <- function(id, celltypes, filename){
    id %>% 
        synapse_file_to_tbl(delim = ",") %>% 
        dplyr::filter(stringr::str_detect(sample.id, "^[RB]M")) %>% 
        dplyr::mutate(cell.type = factor(cell.type, levels = celltypes)) %>% 
        dplyr::group_by(dataset.name) %>%
        dplyr::group_split() %>%
        purrr::map(tidyr::complete, dataset.name, sample.id, cell.type) %>%
        dplyr::bind_rows() %>% 
        readr::write_csv(filename)
}

file_to_synapse <- function(file, id){
    file %>% 
        synapser::File(parent = id) %>% 
        synapser::synStore()
    file.remove(file)
}

purrr::pmap(list(ids, cells, file_names), format_gold_standard)
purrr::walk(file_names, file_to_synapse, parent_dir)




    





