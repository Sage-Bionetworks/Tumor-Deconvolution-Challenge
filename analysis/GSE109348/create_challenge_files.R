library(magrittr)


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


## Update the following:
source("dataset-setup.R")

suppressPackageStartupMessages(p_load(tidyverse))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(synapserutils))
suppressPackageStartupMessages(p_load(data.table))
suppressPackageStartupMessages(p_load(magrittr))
suppressPackageStartupMessages(p_load(Biobase))
suppressPackageStartupMessages(p_load(ImmuneSpaceR))
suppressPackageStartupMessages(p_load(reshape2))
suppressPackageStartupMessages(p_load(GEOquery))



## Begin configuration

file <- "create_challenge_files.R"

## End configuration

script_url <- paste0(url.base, "/", dataset, "/", file)
activity_name <- "create Challenge expression and ground truth files"

source("../../scripts/utils.R")
synLogin()

# mixture files
fine_mix_tbl <- fine_mix_id %>% 
    create_df_from_synapse_id() %>% 
    dplyr::mutate(sample_name = stringr::str_c("s", 1:dplyr::n()))

coarse_mix_tbl <- coarse_mix_id %>% 
    create_df_from_synapse_id() %>% 
    dplyr::mutate(sample_name = stringr::str_c("s", 1:dplyr::n()))



output.folder.synId <- create.folder(staging.folder.synId, dataset)
metadata.file.name <- paste0(dataset, "-metadata.tsv")

if(!grepl(dataset, pattern="GSE")) {
  stop("Was expecting an GSE dataset\n")
}

gse.ids <- list(dataset)
names(gse.ids) <- gse.ids

gses <- llply(gse.ids, .parallel = TRUE, .fun = function(gse.id) getGEO(gse.id, GSEMatrix=TRUE))
gses <- unlist(gses)

## Process ground truth

## Create a data frame holding the ground truth with columns sample, cell.type, and measured
gt.df <- gses %>%
  get.geo.metadata.tbl %>%
  rownames_to_column("sample") %>% 
  as.data.frame()
  # dplyr::select(sample, title, `flow cytometry cell subset proportions:ch1`) %>%
  # set_colnames(c("sample", "id", "cell_types")) %>%
  # filter(cell_types != "NA") %>%
  # separate(cell_types, sep = "; ", into = as.character(1:20), fill = "right") %>%
  # gather(key = "key", value = "value", -c(sample, id)) %>%
  # dplyr::select(-key) %>%
  # dplyr::select(-id) %>%
  # drop_na %>%
  # separate(value, sep = " = ", into = c("cell.type", "measured")) %>%
  # mutate(cell.type = str_replace_all(cell.type, " ", "_")) %>%
  # mutate(cell.type = str_replace_all(cell.type, "ï", "i")) %>%
  # mutate(measured = str_remove_all(measured, "%")) %>%
  # mutate(measured = as.numeric(measured))

## End processing ground truth

## Process GEO expression



## May need to update the following get.geo.platform.name function
platform <- coarse_mix_tbl$platform[[1]]
cancer.type <- coarse_mix_tbl$cancer.type[[1]]
normalization <- coarse_mix_tbl$normalization[[1]]
scale <- coarse_mix_tbl$scale[[1]]
native.probe.type <- coarse_mix_tbl$native.probe.type[[1]]
obfuscated.dataset <- coarse_mix_tbl$dataset.name[[1]]
fine_datatset_name <- stringr::str_c("F", obfuscated.dataset)
coarse_datatset_name <- stringr::str_c("C", obfuscated.dataset)



## Get the GEO expression matrices, combining into one matrix
## First column is 'Gene'
expr.mats <- llply(gses, .fun = function(gse) exprs(gse) %>% as.data.frame)
expr.mat  <- Reduce("cbind", expr.mats)
expr.mat  <- expr.mat %>% 
    drop.duplicate.columns() %>% 
    rownames_to_column(var = "Gene") %>% 
    tidyr::pivot_longer(-Gene) %>% 
    dplyr::arrange(name) %>% 
    dplyr::group_by(name)
    

if (coarse_mix_tbl$normalization.before.summing[[1]] == "divide.by.sample.total"){
    expr.mat <- dplyr::mutate(expr.mat, value = value / sum(value, na.rm = T))
} else {
    stop("no normalization method")
}

expr.mat <- expr.mat %>% 
    dplyr::ungroup() %>% 
    tidyr::pivot_wider(names_from = name, values_from = value)


data.processing <- unlist(get.geo.data.processing(gses))
print(data.processing)

## Extract mappings from probe to gene symbol and Ensembl ID.
probe.to.symbol.map <- get.probe.to.symbol.map(gses)
probe.to.ensg.map <- get.probe.to.ensg.map(gses)

compression.fun <- "colMeans"
symbol.compression.fun <- compression.fun
ensg.compression.fun <- compression.fun
expr.mat.symbol <- translate.genes(expr.mat, probe.to.symbol.map, fun = symbol.compression.fun)
expr.mat.ensg <- translate.genes(expr.mat, probe.to.ensg.map, fun = ensg.compression.fun)

coarse_input_tbl <- dplyr::tibble(
    dataset.name = obfuscated.dataset,
    cancer.type = cancer.type,
    platform = platform,
    scale = scale,
    normalization = normalization,
    native.probe.type = native.probe.type, 
    native.expr.file = stringr::str_c(coarse_datatset_name, "-native-gene-expr.csv"),
    hugo.expr.file = stringr::str_c(coarse_datatset_name, "-hugo-gene-expr.csv"),
    ensg.expr.file = stringr::str_c(coarse_datatset_name, "-ensg-gene-expr.csv"),
    symbol.compression.function = symbol.compression.fun,
    ensg.compression.function = ensg.compression.fun,
    symbol.to.native.mapping.file = stringr::str_c(coarse_datatset_name, "-symbol-to-native-mapping.tsv"),
    ensg.to.native.mapping.file = stringr::str_c(coarse_datatset_name, "-ensg-to-native-mapping.tsv")
)

fine_input_tbl <- dplyr::tibble(
    dataset.name = obfuscated.dataset,
    cancer.type = cancer.type,
    platform = platform,
    scale = scale,
    normalization = normalization,
    native.probe.type = native.probe.type, 
    native.expr.file = stringr::str_c(fine_datatset_name, "-native-gene-expr.csv"),
    hugo.expr.file = stringr::str_c(fine_datatset_name, "-hugo-gene-expr.csv"),
    ensg.expr.file = stringr::str_c(fine_datatset_name, "-ensg-gene-expr.csv"),
    symbol.compression.function = symbol.compression.fun,
    ensg.compression.function = ensg.compression.fun,
    symbol.to.native.mapping.file = stringr::str_c(fine_datatset_name, "-symbol-to-native-mapping.tsv"),
    ensg.to.native.mapping.file = stringr::str_c(fine_datatset_name, "-ensg-to-native-mapping.tsv")
)


expr.mats <- list("native" = expr.mat, "ensg" = expr.mat.ensg, "hugo" = expr.mat.symbol)
mapping.mats <- list("symbol" = probe.to.symbol.map, "ensg" = probe.to.ensg.map)

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

mutate_func <- function(df, col, exp){
    df %>%
        dplyr::mutate(!!col := !!exp) %>%
        dplyr::select("Gene", col)
}

coarse_mix_tbl2 <- coarse_mix_tbl %>% 
    dplyr::select(sample_name, freqs, samples, cell.types) %>% 
    dplyr::mutate(mix = stringr::str_c("S", 1: dplyr::n())) %>% 
    dplyr::mutate(exp = purrr:::map2(samples, freqs, make_expr))

fine_mix_tbl2 <- fine_mix_tbl %>% 
    dplyr::select(sample_name, freqs, samples, cell.types) %>% 
    dplyr::mutate(mix = stringr::str_c("S", 1: dplyr::n())) %>% 
    dplyr::mutate(exp = purrr:::map2(samples, freqs, make_expr))

coarse_native_df <- 
    purrr::map2(
        coarse_mix_tbl2$mix, 
        coarse_mix_tbl2$exp,
        ~mutate_func(expr.mat, .x, .y)
    ) %>%
    purrr::reduce(dplyr::inner_join) 

coarse_ensembl_df <- 
    purrr::map2(
        coarse_mix_tbl2$mix, 
        coarse_mix_tbl2$exp,
        ~mutate_func(expr.mat.ensg, .x, .y)
    ) %>%
    purrr::reduce(dplyr::inner_join) 

coarse_hugo_df <- 
    purrr::map2(
        coarse_mix_tbl2$mix, 
        coarse_mix_tbl2$exp,
        ~mutate_func(expr.mat.symbol, .x, .y)
    ) %>%
    purrr::reduce(dplyr::inner_join) 

fine_native_df <- 
    purrr::map2(
        fine_mix_tbl2$mix, 
        fine_mix_tbl2$exp,
        ~mutate_func(expr.mat, .x, .y)
    ) %>%
    purrr::reduce(dplyr::inner_join) 

fine_ensembl_df <- 
    purrr::map2(
        fine_mix_tbl2$mix, 
        fine_mix_tbl2$exp,
        ~mutate_func(expr.mat.ensg, .x, .y)
    ) %>%
    purrr::reduce(dplyr::inner_join) 

fine_hugo_df <- 
    purrr::map2(
        fine_mix_tbl2$mix, 
        fine_mix_tbl2$exp,
        ~mutate_func(expr.mat.symbol, .x, .y)
    ) %>%
    purrr::reduce(dplyr::inner_join) 

gt_coarse <- coarse_mix_tbl2 %>% 
    dplyr::select(
        sample.id = mix, 
        measured = freqs, 
        cell.type = cell.types
    ) %>% 
    tidyr::separate_rows(measured, cell.type, sep = "; ") %>% 
    dplyr::filter(cell.type != "coarse.noise") %>% 
    dplyr::mutate(
        measured = as.double(measured),
        dataset.name = coarse_datatset_name
    ) %>% 
    dplyr::mutate(
        cell.type = forcats::fct_expand(cell.type, coarse_celltypes)
    ) %>% 
    tidyr::complete(
        dataset.name,
        sample.id,
        cell.type,
        fill = list("measured" = NA)
    ) 


gt_fine <- fine_mix_tbl2 %>% 
    dplyr::select(
        sample.id = mix, 
        measured = freqs, 
        cell.type = cell.types
    ) %>% 
    tidyr::separate_rows(measured, cell.type, sep = "; ") %>% 
    dplyr::filter(cell.type != "fine.noise") %>% 
    dplyr::mutate(
        measured = as.double(measured),
        dataset.name = fine_datatset_name
    ) %>% 
    dplyr::mutate(
        cell.type = forcats::fct_expand(cell.type, fine_celltypes)
    ) %>% 
    tidyr::complete(
        dataset.name,
        sample.id,
        cell.type,
        fill = list("measured" = NA)
    ) 

upload_tbl_to_synapse <- function(tbl, file_name, id, delim){
    readr::write_delim(tbl, file_name, delim)
    file_entity <- synapser::File(path = file_name, parent = id)
    synapser::synStore(file_entity)
}

upload_tbl_to_synapse(gt_fine, "fine_gt.csv", dataset_id, ",")
upload_tbl_to_synapse(gt_coarse, "coarse_gt.csv", dataset_id, ",")

upload_tbl_to_synapse(fine_input_tbl, "fine_input.csv", dataset_id, ",")
upload_tbl_to_synapse(coarse_input_tbl, "coarse_input.csv", dataset_id, ",")

upload_tbl_to_synapse(
    coarse_ensembl_df, 
    stringr::str_c(coarse_datatset_name, "-ensg-gene-expr.csv"), 
    dataset_id, 
    ","
)

upload_tbl_to_synapse(
    fine_ensembl_df, 
    stringr::str_c(fine_datatset_name, "-ensg-gene-expr.csv"), 
    dataset_id, 
    ","
)

upload_tbl_to_synapse(
    coarse_hugo_df, 
    stringr::str_c(coarse_datatset_name, "-hugo-gene-expr.csv"), 
    dataset_id, 
    ","
)

upload_tbl_to_synapse(
    fine_hugo_df, 
    stringr::str_c(fine_datatset_name, "-hugo-gene-expr.csv"), 
    dataset_id, 
    ","
)

upload_tbl_to_synapse(
    coarse_native_df, 
    stringr::str_c(coarse_datatset_name, "-native-gene-expr.csv"), 
    dataset_id, 
    ","
)

upload_tbl_to_synapse(
    fine_native_df, 
    stringr::str_c(fine_datatset_name, "-native-gene-expr.csv"), 
    dataset_id, 
    ","
)

upload_tbl_to_synapse(
    probe.to.symbol.map, 
    stringr::str_c(coarse_datatset_name, "-symbol-to-native-mapping.tsv"), 
    dataset_id, 
    "\t"
)

upload_tbl_to_synapse(
    probe.to.symbol.map, 
    stringr::str_c(fine_datatset_name, "-symbol-to-native-mapping.tsv"), 
    dataset_id, 
    "\t"
)

upload_tbl_to_synapse(
    probe.to.ensg.map, 
    stringr::str_c(coarse_datatset_name, "-ensg-to-native-mapping.tsv"), 
    dataset_id, 
    "\t"
)

upload_tbl_to_synapse(
    probe.to.ensg.map, 
    stringr::str_c(fine_datatset_name, "-ensg-to-native-mapping.tsv"), 
    dataset_id, 
    "\t"
)


