library(synapser)
library(tidyverse)

synapser::synLogin()

devtools::source_url(
    "https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R"
)

insillico_dir_id        <- "syn22013292"
upload_dir_id           <- "syn21821096"
dataset_file_id         <- "syn21590362"
invitro_input_file_id   <- "syn21782475"
insilico_input_file_id  <- "syn22013441"
translation_id          <- "syn21574276"

input_cols <- c(
    "dataset.name",
    "cancer.type",
    "platform",
    "scale", 
    "normalization", 
    "native.probe.type", 
    "native.expr.file",
    "hugo.expr.file", 
    "ensg.expr.file", 
    "symbol.compression.function", 
    "ensg.compression.function",
    "native.expr.est.counts.file",
    "hugo.expr.est.counts.file", 
    "ensg.expr.est.counts.file", 
    "symbol.compression.est.counts.function", 
    "ensg.compression.est.counts.function", 
    "symbol.to.native.mapping.file", 
    "ensg.to.native.mapping.file", 
    "fastq1.files", 
    "fastq2.files", 
    "fastq.samples"
)

expr_file_tbl <- dplyr::tibble(
    type = c("hugo_tpm", "hugo_counts", "ensg_tpm", "ensg_counts"),
    id   = c("syn21576632", "syn21576630", "syn21576631", "syn21576629")
)

dataset_tbl   <- dataset_file_id %>% 
    synapse_file_to_tbl() %>% 
    dplyr::filter(dataset %in%  c("DS1", "DS2", "DS3", "DS4")) %>% 
    dplyr::select(dataset, id) %>% 
    dplyr::group_by(dataset) %>% 
    dplyr::summarize(samples = list(id)) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(dataset) 

combined_tbl <- 
    merge(dataset_tbl, expr_file_tbl) %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(filename = paste0(type, "_", dataset, ".csv"))

upload_expression_file <- function(id, samples, filename){
    id %>% 
        synapse_file_to_tbl(delim = ",") %>% 
        dplyr::select(c("Gene", samples)) %>% 
        readr::write_csv(filename)
    
    filename %>% 
        synapser::File(parent = upload_dir_id) %>% 
        synapser::synStore()
    
}


# invitro expression files ----------------------------------------------------

# param_tbl <- combined_tbl %>%
#     dplyr::select(id, samples, filename) %>%
#     pmap(upload_expression_file)

# trnaslation files------------------------------------------------------------

# translation_tbl <- "syn21574276" %>% 
#     synapse_file_to_tbl() %>% 
#     dplyr::select(ensg = ensembl_gene_id, hugo = hgnc_symbol) %>% 
#     dplyr::distinct() 
# 
# translation_tbl %>%
#     tidyr::drop_na() %>%
#     dplyr::select(
#         from = ensg,
#         to = hugo
#     ) %>%
#     readr::write_tsv("native_to_hugo.tsv")
# 
# "native_to_hugo.tsv" %>% 
#     synapser::File(parent = upload_dir_id) %>% 
#     synapser::synStore()
# 
# translation_tbl %>% 
#     dplyr::select(from = ensg) %>% 
#     dplyr::mutate(to = from) %>% 
#     dplyr::distinct() %>% 
#     tidyr::drop_na() %>% 
#     readr::write_tsv("native_to_ensg.tsv")
# 
# "native_to_ensg.tsv" %>% 
#     synapser::File(parent = upload_dir_id) %>% 
#     synapser::synStore()
    
# ------------------------------------------------------------


expr_file_tbl2 <- combined_tbl %>% 
    dplyr::select(dataset, type, filename, samples) %>% 
    tidyr::pivot_wider(names_from = type, values_from = filename) %>% 
    tidyr::unnest(cols = c(samples)) %>% 
    dplyr::mutate(
        fastq1.files = paste0(samples, "_R1.fastq.gz "),
        fastq2.files = paste0(samples, "_R2.fastq.gz ")
    ) %>% 
    dplyr::group_by(dataset, hugo_tpm, hugo_counts, ensg_tpm, ensg_counts) %>% 
    dplyr::summarise(
        fastq.samples = paste0(samples, collapse = ","),
        fastq1.files = paste0(fastq1.files, collapse = ","),
        fastq2.files = paste0(fastq2.files, collapse = ",")
    ) %>% 
    dplyr::rename(
        dataset.name = dataset,
        hugo.expr.file = hugo_tpm,
        ensg.expr.file = ensg_tpm,
        hugo.expr.est.counts.file = hugo_counts,
        ensg.expr.est.counts.file = ensg_counts
    ) %>% 
    dplyr::mutate(
        native.expr.file = ensg.expr.file,
        native.expr.est.counts.file = ensg.expr.est.counts.file
    )
    
invitro_input_tbl <- invitro_input_file_id %>% 
    synapse_file_to_tbl(delim = ",") %>% 
    dplyr::filter(dataset.name %in% c("DS1", "DS2", "DS3", "DS4")) %>% 
    dplyr::select(-c(
        "native.expr.file", "hugo.expr.file", "ensg.expr.file",
        "fastq.samples", "fastq1.files", "fastq2.files"
    )) %>% 
    dplyr::inner_join(expr_file_tbl2) %>% 
    dplyr::mutate(
        symbol.compression.est.counts.function = "colMeans",
        ensg.compression.est.counts.function = "colMeans",
        symbol.to.native.mapping.file = "native_to_hugo.tsv",
        ensg.to.native.mapping.file = "native_to_ensg.tsv"
    )

insilico_input_tbl <- insilico_input_file_id %>% 
    synapse_file_to_tbl(delim = ",") 

dplyr::bind_rows(invitro_input_tbl, insilico_input_tbl) %>% 
    readr::write_csv("input.csv")

"input.csv" %>% 
    synapser::File(parent = upload_dir_id) %>% 
    synapser::synStore()

input_cols[!input_cols %in% colnames(input_tbl)]
colnames(input_tbl)[!colnames(input_tbl) %in% input_cols]

# insilico expression files ----------------------------------------------------

insilico_expr_files <- insilico_input_tbl %>% 
    dplyr::select("dataset.name", dplyr::contains("expr")) %>% 
    tidyr::pivot_longer(., -"dataset.name", names_to = "x", values_to = "y") %>% 
    dplyr::select("y") %>% 
    dplyr::distinct() %>% 
    dplyr::pull("y")

insillico_dir_id %>% 
    create_entity_tbl(.) %>% 
    dplyr::filter(.data$name %in% insilico_expr_files) %>% 
    dplyr::pull("id") %>% 
    purrr::map(synapserutils::copy, upload_dir_id)


