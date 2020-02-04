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
## NB: need to load mygene before synapser to avoid an error
suppressPackageStartupMessages(p_load(mygene))
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


expr.mat  <- expr_id %>% 
    create_df_from_synapse_id(unzip = T) %>% 
    dplyr::rename(Gene = V1) %>% 
    tidyr::pivot_longer(-Gene) %>% 
    dplyr::arrange(name) 

translation_df <- expr.mat %>% 
    dplyr::select(name) %>% 
    dplyr::distinct() %>% 
    dplyr::mutate(new_name = stringr::str_c("input", 1:dplyr::n()))

expr.mat <- expr.mat %>% 
    dplyr::inner_join(translation_df) %>% 
    dplyr::select(-name) %>% 
    dplyr::rename(name = new_name)

# mixture files
fine_mix_tbl <- fine_mix_id %>% 
    create_df_from_synapse_id() %>% 
    tidyr::separate_rows(samples, sep = "; ") %>% 
    dplyr::inner_join(translation_df, by = c("samples" = "name")) %>% 
    dplyr::select(-samples) %>% 
    dplyr::rename(samples = new_name) %>% 
    dplyr::group_by_at(vars(-samples)) %>% 
    dplyr::summarise(samples = str_c(samples, collapse = "; ")) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(sample_name = stringr::str_c("S", 1:dplyr::n()))

coarse_mix_tbl <- coarse_mix_id %>% 
    create_df_from_synapse_id() %>% 
    tidyr::separate_rows(samples, sep = "; ") %>% 
    dplyr::inner_join(translation_df, by = c("samples" = "name")) %>% 
    dplyr::select(-samples) %>% 
    dplyr::rename(samples = new_name) %>% 
    dplyr::group_by_at(vars(-samples)) %>% 
    dplyr::summarise(samples = str_c(samples, collapse = "; ")) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(sample_name = stringr::str_c("S", 1:dplyr::n()))


# if(!grepl(dataset, pattern="GSE")) {
#   stop("Was expecting an GSE dataset\n")
# }
# 
# gse.ids <- list(dataset)
# names(gse.ids) <- gse.ids
# 
# gses <- llply(gse.ids, .parallel = TRUE, .fun = function(gse.id) getGEO(gse.id, GSEMatrix=TRUE))
# gses <- unlist(gses)
# 




platform <- coarse_mix_tbl$platform[[1]]
cancer.type <- coarse_mix_tbl$cancer.type[[1]]
normalization <- coarse_mix_tbl$normalization[[1]]
scale <- coarse_mix_tbl$scale[[1]]
native.probe.type <- coarse_mix_tbl$native.probe.type[[1]]
obfuscated.dataset <- coarse_mix_tbl$dataset.name[[1]]
normalization_method <- coarse_mix_tbl$normalization.before.summing[[1]]
fine_datatset_name <- stringr::str_c("F", obfuscated.dataset)
coarse_datatset_name <- stringr::str_c("C", obfuscated.dataset)




    

if (normalization_method == "divide.by.sample.total"){
    expr.mat <- dplyr::mutate(expr.mat, value = value / sum(value, na.rm = T))
} else if (normalization_method == "divide.by.sample.total.multiply.by.million"){
    expr.mat <- dplyr::mutate(expr.mat, value = 1000000 * value / sum(value, na.rm = T))
} else if (normalization_method == "exponentiate.base.2"){
    expr.mat <- dplyr::mutate(expr.mat, value = 2^value)
} else if (normalization_method == "None"){
    
} else {
    stop("no normalization method")
}

expr.mat <- expr.mat %>% 
    dplyr::ungroup() %>% 
    tidyr::pivot_wider(names_from = name, values_from = value)


# data.processing <- unlist(get.geo.data.processing(gses))
# print(data.processing)

## Extract mappings from probe to gene symbol and Ensembl ID.
ensg.to.symbol.map <- get.ensg.to.sym.map(expr.mat$Gene)
ensg.to.ensg.map <- dplyr::mutate(ensg.to.symbol.map, to = from)

compression.fun <- "colMeans"
symbol.compression.fun <- compression.fun
ensg.compression.fun <- compression.fun

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



expr.mat.symbol <- translate.genes(expr.mat, ensg.to.symbol.map, fun = symbol.compression.fun)
expr.mat.ensg <- translate.genes(expr.mat, ensg.to.ensg.map, fun = ensg.compression.fun)


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

eps <- 10^-4

## Confirm all freqs sum to one
l_ply(1:nrow(coarse_mix_tbl),
      .fun = function(i) {
               freqs <- unlist(strsplit(as.character(coarse_mix_tbl[i, "freqs"]), split=";[ ]*"))
	       freqs <- as.numeric(freqs)
	       sm <- sum(freqs)
	       if(abs(sm - 1) > eps) { stop(paste0("coarse sum = ", sm, "\n")) }
             })
cat("Confirmed all coarse frequencies sum to one\n")

l_ply(1:nrow(fine_mix_tbl),
      .fun = function(i) {
               freqs <- unlist(strsplit(as.character(fine_mix_tbl[i, "freqs"]), split=";[ ]*"))
	       freqs <- as.numeric(freqs)
	       sm <- sum(freqs)
	       if(abs(sm - 1) > eps) { stop(paste0("fine sum = ", sm, "\n")) }
             })
cat("Confirmed all fine frequencies sum to one\n")



coarse_mix_tbl2 <- coarse_mix_tbl %>% 
    dplyr::select(sample_name, freqs, samples, cell.types) %>% 
    dplyr::mutate(mix = stringr::str_c("S", 1: dplyr::n())) %>% 
    dplyr::mutate(exp = purrr:::map2(samples, freqs, make_expr))

fine_mix_tbl2 <- fine_mix_tbl %>% 
    dplyr::select(sample_name, freqs, samples, cell.types) %>% 
    dplyr::mutate(mix = stringr::str_c("S", 1: dplyr::n())) %>% 
    dplyr::mutate(exp = purrr:::map2(samples, freqs, make_expr))

if (normalization_method == "exponentiate.base.2"){
    post_func <- log2
} else if (normalization_method %in% c("None", "divide.by.sample.total", "divide.by.sample.total.multiply.by.million")){
    post_func <- identity
} else {
    stop("no normalization method")
}


coarse_native_df <- 
    purrr::map2(
        coarse_mix_tbl2$mix, 
        coarse_mix_tbl2$exp,
        ~mutate_func(expr.mat, .x, .y)
    ) %>%
    purrr::reduce(dplyr::inner_join) %>% 
    dplyr::mutate_if(is.numeric, post_func)

coarse_ensembl_df <- 
    purrr::map2(
        coarse_mix_tbl2$mix, 
        coarse_mix_tbl2$exp,
        ~mutate_func(expr.mat.ensg, .x, .y)
    ) %>%
    purrr::reduce(dplyr::inner_join) %>% 
    dplyr::mutate_if(is.numeric, post_func)

coarse_hugo_df <- 
    purrr::map2(
        coarse_mix_tbl2$mix, 
        coarse_mix_tbl2$exp,
        ~mutate_func(expr.mat.symbol, .x, .y)
    ) %>%
    purrr::reduce(dplyr::inner_join) %>% 
    dplyr::mutate_if(is.numeric, post_func)

fine_native_df <- 
    purrr::map2(
        fine_mix_tbl2$mix, 
        fine_mix_tbl2$exp,
        ~mutate_func(expr.mat, .x, .y)
    ) %>%
    purrr::reduce(dplyr::inner_join) %>% 
    dplyr::mutate_if(is.numeric, post_func)

fine_ensembl_df <- 
    purrr::map2(
        fine_mix_tbl2$mix, 
        fine_mix_tbl2$exp,
        ~mutate_func(expr.mat.ensg, .x, .y)
    ) %>%
    purrr::reduce(dplyr::inner_join) %>% 
    dplyr::mutate_if(is.numeric, post_func)

fine_hugo_df <- 
    purrr::map2(
        fine_mix_tbl2$mix, 
        fine_mix_tbl2$exp,
        ~mutate_func(expr.mat.symbol, .x, .y) 
    ) %>%
    purrr::reduce(dplyr::inner_join) %>% 
    dplyr::mutate_if(is.numeric, post_func)

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
    plyr::ddply(.variables = c("dataset.name", "sample.id", "cell.type"),
    		.fun = function(df) data.frame(measured = sum(df$measured))
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
    plyr::ddply(.variables = c("dataset.name", "sample.id", "cell.type"),
    		.fun = function(df) data.frame(measured = sum(df$measured))
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

upload_tbl_to_synapse(fine_input_tbl, "fine_input.csv", dataset_id, ",")
upload_tbl_to_synapse(coarse_input_tbl, "coarse_input.csv", dataset_id, ",")

## Confirm ground truths sum to one and that no cell types are duplicated within a sample
d_ply(gt_coarse,
      .variables = c("dataset.name", "sample.id"),
      .fun = function(df) {
	       sm <- sum(na.omit(df$measured))
	       if(abs(sm - 1) > eps) { stop(paste0("coarse ground truth sum = ", sm, "\n")) }
	       dups <- duplicated(df$cell.type)
	       if(any(dups)) { stop(paste0("duplicated cell types: ", paste(df$cell.type[dupgs], collapse=", "), "\n")) }
             })
cat("Confirmed all coarse ground truth frequencies sum to one and no duplicated cell types\n")

d_ply(gt_fine,
      .variables = c("dataset.name", "sample.id"),
      .fun = function(df) {
	       sm <- sum(na.omit(df$measured))      
	       if(abs(sm - 1) > eps) { stop(paste0("fine ground truth sum = ", sm, "\n")) }
	       dups <- duplicated(df$cell.type)
	       if(any(dups)) { stop(paste0("duplicated cell types: ", paste(df$cell.type[dupgs], collapse=", "), "\n")) }
             })
cat("Confirmed all fine ground truth frequencies sum to one and no duplicated cell types\n")

upload_tbl_to_synapse(gt_fine, "fine_gt.csv", dataset_id, ",")
upload_tbl_to_synapse(gt_coarse, "coarse_gt.csv", dataset_id, ",")

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
    ensg.to.symbol.map, 
    stringr::str_c(coarse_datatset_name, "-symbol-to-native-mapping.tsv"), 
    dataset_id, 
    "\t"
)


upload_tbl_to_synapse(
    ensg.to.ensg.map, 
    stringr::str_c(coarse_datatset_name, "-ensg-to-native-mapping.tsv"), 
    dataset_id, 
    "\t"
)




