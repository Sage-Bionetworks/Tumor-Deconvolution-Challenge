library(synapser)
library(tidyverse)

devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)

synapser::synLogin()

fastq_dir_id   <- "syn21557721"
md5sum_file_id <- "syn21570548"
fileview_id    <- "syn21570575"

md5sum_tbl <-  md5sum_file_id %>% 
    synapse_file_to_tbl(delim = "   ", col_names = c("md5", "name")) %>%
    dplyr::mutate(
        md5 = stringr::str_remove_all(md5, " "),
        name = stringr::str_remove_all(name, " ")
    )


fileview_tbl <- "SELECT id, name FROM syn21570575" %>% 
    query_synapse_table() %>% 
    dplyr::mutate(
        entity = purrr::map(id, synapser::synGet, downloadFile = F),
        size = purrr::map_dbl(entity, ~ .x$get('fileSize')),
        md5 = purrr::map_chr(entity, ~ .x$get('md5')),
    ) 

check_md5sum_tbl <- dplyr::full_join(fileview_tbl, md5sum_tbl)

if (!nrow(fileview_tbl) == nrow(check_md5sum_tbl)) {
    stop("filenames or md5s don't match")
}

annotation_tbl <- fileview_tbl %>% 
    dplyr::select(entity = id, name, size) %>% 
    dplyr::mutate(file = stringr::str_remove_all(name, ".fastq.gz")) %>% 
    dplyr::select(-name) %>% 
    dplyr::mutate(
        pair = stringr::str_match(file, "_(R[12]$)")[,2],
        sample = stringr::str_remove_all(file, "_(R[12]$)")
    ) %>% 
    dplyr::select(-file) %>% 
    dplyr::mutate(
        sample2   = stringr::str_remove_all(sample, "_([12]$)"),
        replicate = stringr::str_match(sample, "_([12]$)")[,2],
        type      = dplyr::if_else(
            str_detect(sample2, "^BM"),
            "Biological Mix",
            dplyr::if_else(
                str_detect(sample2, "^RM"),
                "Random Mix",
                "Purified Cell"
            )
        ),
        celltype = dplyr::if_else(
            str_detect(sample2, "^[BR]M"),
            NA_character_,
            sample2
        )
    ) %>% 
    dplyr::select(-sample2) %>% 
    tidyr::nest(annotations = -entity)

purrr::pmap(annotation_tbl, synapser::synSetAnnotations)



