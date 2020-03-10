library(synapser)
library(tidyverse)

devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)

synapser::synLogin()

fastq_dir_id   <- "syn21557721"
md5sum_file_id <- "syn21570548"
fileview_id    <- "syn21570575"

rm_mix_reads <- "syn21614530" %>% 
    synapse_file_to_tbl(delim = " ", col_names = c("lines", "name")) %>% 
    dplyr::mutate(
        sample = stringr::str_match(name, "^(RM[0-9]+)_R")[,2],
        pair = stringr::str_match(name, "_(R[12]).")[,2],
        reads = as.integer(lines / 4)
    ) %>% 
    dplyr::select(-c(lines, name)) %>% 
    distinct()

bm_mix_reads <- "syn21614539" %>% 
    synapse_file_to_tbl(delim = " ", col_names = c("lines", "name")) %>% 
    dplyr::mutate(
        sample = stringr::str_match(name, "^(BM[0-9]+)_R")[,2],
        pair = stringr::str_match(name, "_(R[12]).")[,2],
        reads = as.integer(lines / 4)
    ) %>% 
    dplyr::select(-c(lines, name)) %>% 
    distinct()
    

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
    dplyr::select(-file)

pc_table <- annotation_tbl %>% 
    dplyr::filter(!str_detect(sample, "^[RB]M")) %>% 
    dplyr::mutate(type = "Purified Cell")

pc_table1 <- pc_table %>% 
    dplyr::filter(!stringr::str_detect(sample, "_[12]$")) %>% 
    dplyr::mutate(cell_type = sample) %>% 
    tidyr::nest(annotations = -entity)

pc_table2 <- pc_table %>% 
    dplyr::filter(stringr::str_detect(sample, "_[12]$")) %>% 
    dplyr::mutate(
        replicate = stringr::str_match(sample, "_([12]$)")[,2],
        cell_type =  stringr::str_remove_all(sample, "_([12]$)"),
    ) %>% 
    tidyr::nest(annotations = -entity)

bm_mix_table <- annotation_tbl %>% 
    dplyr::filter(str_detect(sample, "^BM")) %>% 
    dplyr::mutate(type = "Biological Mix") %>% 
    dplyr::left_join(bm_mix_reads) %>% 
    tidyr::nest(annotations = -entity)

rm_mix_table <- annotation_tbl %>% 
    dplyr::filter(str_detect(sample, "^RM")) %>% 
    dplyr::mutate(type = "Random Mix") %>% 
    dplyr::left_join(rm_mix_reads) %>% 
    tidyr::nest(annotations = -entity)

res <-
    list(pc_table1, pc_table2, bm_mix_table, rm_mix_table) %>% 
    purrr::map(purrr::pmap, synapser::synSetAnnotations)




