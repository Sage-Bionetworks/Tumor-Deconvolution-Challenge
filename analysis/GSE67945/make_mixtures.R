library(magrittr)
synapser::synLogin()
devtools::source_url("https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R")


mixture_id     <- "syn20742216"
count_dir_id   <- "syn20742192"
destination_id <- "syn20742190"
dataset_name   <- "DS500"

params <- list(
    "dataset.name" = dataset_name,
    "cancer.type" = "CRC",
    "platform" = "Illumina",
    "scale" = "Linear",
    "normalization" = "CPM",
    "native.probe.type" = "ENSG",
    "hugo.expr.file" = stringr::str_c(
        dataset_name,
        "-hugo-gene-expr.csv"
    ),
    "ensg.expr.file" = stringr::str_c(
        dataset_name, 
        "-ensg-gene-expr.csv"
    ),
    "native.expr.file" = stringr::str_c(
        dataset_name, 
        "-ensg-gene-expr.csv"
    ),
    "symbol.compression.function" = "colMeans",
    "ensg.compression.function" = "identity",
    "symbol.to.native.mapping.file" = stringr::str_c(
        dataset_name, 
        "-symbol-to-native-mapping.tsv"
    ),
    "ensg.to.native.mapping.file" = stringr::str_c(
        dataset_name, 
        "-ensg-to-native-mapping.tsv"
    ),
    "fastq1.files" = NA,
    "fastq2.files" = NA,
    "fastq.samples" = NA
)


coarse_celltypes <- c(
    "B.cells",
    "CD4.T.cells",
    "CD8.T.cells",
    "NK.cells",
    "neutrophils",
    "monocytic.lineage",
    "fibroblasts",
    "endothelial.cells"
)

fine_celltypes <- c(
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

upload_file_to_synapse <- function(path, destination_id){
    path %>% 
        synapser::File(destination_id) %>% 
        synapser::synStore()
    
    file.remove(path)
}


gene_df <- 
    biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl") %>% 
    biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), mart = .)

gene_df %>% 
    dplyr::select(from = ensembl_gene_id) %>% 
    dplyr::mutate(to = from) %>% 
    readr::write_tsv(params$ensg.to.native.mapping.file)

upload_file_to_synapse(params$ensg.to.native.mapping.file, destination_id)

gene_df %>% 
    dplyr::select(from = ensembl_gene_id, to = hgnc_symbol) %>% 
    dplyr::filter(!to == "") %>% 
    readr::write_tsv(params$symbol.to.native.mapping.file)

upload_file_to_synapse(params$symbol.to.native.mapping.file, destination_id)

make_expr <- function(samples, freqs){
    exp <-  
        stringr::str_c(
            unlist(stringr::str_split(samples, pattern = "; ")),
            unlist(stringr::str_split(freqs, pattern = "; ")),
            sep = " * ", 
            collapse = " + "
    ) %>% 
        rlang::parse_expr()
}

mixture_df <- mixture_id %>% 
    synapse_file_to_tbl() %>% 
    dplyr::select(freqs, samples) %>% 
    dplyr::mutate(mix = stringr::str_c("S", 1: dplyr::n())) %>% 
    dplyr::mutate(exp = purrr:::map2(samples, freqs, make_expr))

file_names <- count_dir_id %>% 
    create_entity_tbl() %>% 
    dplyr::pull(id) %>% 
    purrr::map(synapser::synGet) %>% 
    purrr::map_chr("path")

sample_names <- file_names %>% 
    basename() %>% 
    stringr::str_remove_all("^GSM[0-9]+_") %>% 
    stringr::str_remove_all(".gene_counts.txt$")

mutate_func <- function(df, col, exp){
    df %>%
        dplyr::mutate(!!col := !!exp) %>%
        dplyr::select("ensembl_gene_id", col)
}

cpm_df <-
    purrr::map2(
        file_names,
        sample_names,
        ~readr::read_tsv(.x, col_names = c("ensembl_gene_id", .y))
    ) %>%
    purrr::reduce(dplyr::full_join) %>% 
    dplyr::filter(stringr::str_detect(ensembl_gene_id, "^ENSG[:digit:]+")) %>%
    tidyr::gather(sample, counts, -ensembl_gene_id) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(cpm = counts * 10^6 / sum(counts)) %>%
    dplyr::select(-counts) %>% 
    dplyr::ungroup() %>%
    tidyr::spread(sample, cpm) %>% 
    dplyr::rename()

ensembl_df <- 
    purrr::map2(
        mixture_df$mix, 
        mixture_df$exp,
        ~mutate_func(cpm_df, .x, .y)
    ) %>%
    purrr::reduce(dplyr::inner_join) 

ensembl_df %>%
    dplyr::rename(Gene = ensembl_gene_id) %>% 
    readr::write_csv(params$ensg.expr.file)

upload_file_to_synapse(params$ensg.expr.file, destination_id)

hugo_df <- ensembl_df %>% 
    dplyr::inner_join(gene_df) %>% 
    dplyr::select(-ensembl_gene_id) %>%
    dplyr::select(hgnc_symbol, dplyr::everything()) %>% 
    dplyr::filter(!hgnc_symbol == "") %>% 
    dplyr::group_by(hgnc_symbol) %>%
    dplyr::summarise_all(mean) %>%
    dplyr::ungroup() %>% 
    dplyr::rename(Gene = hgnc_symbol) %>% 
    readr::write_csv(params$hugo.expr.file)

upload_file_to_synapse(params$hugo.expr.file, destination_id)


gt_df <- mixture_df %>%
    dplyr::select(-exp) %>%
    tidyr::separate_rows(freqs, samples, sep = "; ") %>%
    dplyr::mutate(cell.type = dplyr::if_else(
        stringr::str_detect(samples, "HU"),
        "endothelial.cells",
        "fibroblasts"
    )) %>%
    dplyr::select(-samples) %>%
    dplyr::rename(sample.id = mix, measured = freqs) %>% 
    dplyr::mutate(dataset.name = params$dataset.name) %>% 
    dplyr::select(
        dataset.name,
        sample.id,
        cell.type,
        measured) %>% 
    dplyr::mutate(cell.type = as.factor(cell.type))
    

gt_df %>% 
    dplyr::mutate(
        cell.type = forcats::fct_expand(cell.type, coarse_celltypes)
    ) %>% 
    tidyr::complete(
        dataset.name,
        sample.id,
        cell.type,
        fill = list("measured" = NA)
    ) %>% 
    readr::write_csv("coarse_gt.csv")

upload_file_to_synapse("coarse_gt.csv", destination_id)

gt_df %>% 
    dplyr::mutate(
        cell.type = forcats::fct_expand(cell.type, fine_celltypes)
    ) %>% 
    tidyr::complete(
        dataset.name,
        sample.id,
        cell.type,
        fill = list("measured" = NA)
    ) %>% 
    readr::write_csv("fine_gt.csv")

upload_file_to_synapse("fine_gt.csv", destination_id)

params %>% 
    data.frame() %>% 
    readr::write_csv("input.csv") 

upload_file_to_synapse("input.csv", destination_id)





