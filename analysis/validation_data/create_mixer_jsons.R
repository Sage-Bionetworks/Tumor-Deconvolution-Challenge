library(magrittr)

synapser::synLogin()
devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)

parameter_list <- list(
    "kallisto_index" = list(
        "path" = "Homo_sapiens.GRCh38.cdna.all.idx",
        "class" = "File"),
    "synapse_config" = list(
        "path" = ".synapseConfig",
        "class" = "File"),
    "fastq_destination_id" = "syn21608643",
    "kallisto_destination_id" = "syn21608711"
)


uploaded_samples <- "select sample from syn21609795" %>%
    query_synapse_table() %>%
    dplyr::pull(sample)

cell_fastq_tbl <- 
    paste(
        "SELECT name, pair, sample FROM syn21570575",
        "WHERE type = 'Purified Cell'"
    ) %>% 
    query_synapse_table(.) %>% 
    tidyr::pivot_wider(., names_from = pair, values_from = name) %>% 
    dplyr::rename(cell_type = sample)

invitro_mixture_tbl <- "SELECT sample, size, reads FROM syn21570575" %>%
    query_synapse_table(.) %>%  
    dplyr::group_by(sample) %>% 
    dplyr::summarise(
        size = sum(size),
        reads = max(reads)
    ) %>% 
    dplyr::ungroup()


insilico_mixture_tbl <- "syn21608888" %>% 
    synapse_file_to_tbl(.) %>% 
    dplyr::filter(actual > 0) %>% 
    dplyr::rename(percent = actual, cell_type = cell.type) 

param_tbl <- invitro_mixture_tbl %>% 
    dplyr::arrange(sample) %>% 
    dplyr::filter(sample %in% insilico_mixture_tbl$sample) %>% 
    dplyr::filter(!sample %in% uploaded_samples) %>% 
    dplyr::filter(!is.na(reads)) %>% 
    dplyr::mutate(cum_size = cumsum(size)) %>% 
    dplyr::filter(cum_size < 100000000000) %>% 
    dplyr::select(-c(size, cum_size)) %>% 
    dplyr::inner_join(insilico_mixture_tbl, by = "sample") %>% 
    dplyr::inner_join(cell_fastq_tbl, by = "cell_type") %>% 
    dplyr::mutate(
        total = as.integer(reads * percent),
        R1 = paste0("fastqs/", R1),
        R2 = paste0("fastqs/", R2)
    ) %>% 
    dplyr::mutate(
        R1 = stringr::str_remove_all(R1, ".gz"),
        R2 = stringr::str_remove_all(R2, ".gz")
    ) %>% 
    dplyr::select(
        p1_fastq_arrays   = R1,
        p2_fastq_arrays   = R2,
        total_read_arrays = total,
        sample
    ) %>% 
    dplyr::mutate(
        p1_fastq_arrays = purrr::map(
            p1_fastq_arrays, ~list("class" = "File", "path" = .x)
        ),
        p2_fastq_arrays = purrr::map(
            p2_fastq_arrays, ~list("class" = "File", "path" = .x)
        ),
    ) %>% 
    dplyr::group_by(sample) %>% 
    dplyr::summarise_all(list) %>% 
    dplyr::mutate(
        p1_fastq_output_names = paste0(sample, "_p1.fastq"),
        p2_fastq_output_names = paste0(sample, "_p2.fastq"),
        kallisto_output_names = paste0(sample, ".tsv")
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-sample) %>% 
    as.list() %>% 
    c(parameter_list) %>% 
    RJSONIO::toJSON() %>% 
    writeLines("mixer_bm2.json")
 