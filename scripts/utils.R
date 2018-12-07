

# synapse ---------------------------------------------------------------------

require(synapser)
require(data.table)
require(magrittr)
require(dplyr)
require(stringr)

create_df_from_synapse_id <- function(syn_id, location = NULL, unzip = F, ...){
    path <- download_from_synapse(syn_id, location)
    ## This doesn't work on OSX, where zcat expects the file to have .Z suffix
    ##    if(unzip) path <- stringr::str_c("zcat ", path)
    if(unzip) path <- stringr::str_c("zcat < ", path)
    path %>% 
        data.table::fread(...) %>% 
        dplyr::as_data_frame() 
}

download_from_synapse <- function(syn_id, location = NULL){
    path = synapser::synGet(syn_id, downloadLocation = location)$path
    return(path)
}

upload_file_to_synapse <- function(
    path, synapse_id, 
    annotation_list = NULL, 
    activity_obj = NULL, 
    ret = "entity"){
    
    entity <- synapser::File(
        path = path, 
        parent = synapse_id, 
        annotations = annotation_list)
    entity <- synapser::synStore(entity, activity = activity_obj)
    if(ret == "entity") return(entity)
    if(ret == "syn_id") return(entity$properties$id)
}

get_file_df_from_synapse_dir_id <- function(syn_id){
    str_c('select id, name from file where parentId=="', syn_id, '"') %>% 
        synapser::synQuery() %>%
        magrittr::use_series("results") %>% 
        purrr::map(data.frame) %>% 
        dplyr::bind_rows() %>% 
        tibble::as_data_frame() 
}


# data_frame / matrix ---------------------------------------------------------

require(magrittr)
require(tibble)
require(plyr)
require(purrr)

transpose_df <- function(df, id_column, new_col){
    df %>% 
        df_to_matrix(id_column) %>%  
        t() %>% 
        matrix_to_df(new_col)
}

df_to_matrix <- function(df, id_column){
    df %>% 
        data.frame() %>% 
        tibble::column_to_rownames(id_column) %>% 
        as.matrix()
}

matrix_to_df <- function(matrix, new_col){
    matrix %>% 
        data.frame() %>% 
        tibble::rownames_to_column(new_col) %>% 
        tibble::as_data_frame()
}


aggregate_matrix <- function(
    matrix, grouping_df, group_by, group_accross, 
    aggregate_fun = mean,
    split_by_cols = T,
    apply_margin = 1,
    combine_fun = cbind,
    parallel = T){
    
    group_list <- create_group_list(grouping_df, group_by, group_accross)
    matrix %>%
        split_matrix(group_list, split_by_cols, parallel) %>% 
        plyr::llply(apply, apply_margin, aggregate_fun, .parallel = parallel) %>% 
        purrr::invoke(combine_fun, .) 
}

library(plyr)
## Translate/compress genes from one name space (e.g., probe ids) to another (e.g., symbols)
compressGenes <- function(e, mapping, from.col = "from", to.col = "to", fun = mean)
{
  e$to    <- mapping[match(rownames(e), mapping[, from.col]), to.col]
  e           <- e[!is.na(e$to),]
  e           <- ddply(.data = e, .variables = "to", .fun = function(x){apply(x[,-ncol(x)],2,fun)},.parallel = T)
  rownames(e) <- e$to
  e           <- e[,-1]
  return(e)
}

split_matrix <- function(matrix, lst, by_cols = T, parallel = T){
    fun <- ifelse(T, 
                  function(names) matrix[,names, drop = F], 
                  function(names) matrix[names, drop = F])
    plyr::llply(lst, fun, .parallel = parallel)
}


create_group_list <- function(df, group_by, group_accross){
    df %>%
        magrittr::extract2(group_by) %>% 
        split(df, .) %>%
        purrr::map(extract2, group_accross)
}


# misc ------------------------------------------------------------------------

require(magrittr)
require(preprocessCore)

get_summary_by_matrix_cols <- function(columns, matrix, fn){
    m <- matrix[,columns]
    if(length(columns) == 1) return(m)
    apply(m, 1, fn)
}

calculate_cpm <- function(counts){
    1000000 * (counts/sum(counts))
}

zscore_matrix <- function(matrix){
    matrix %>% 
        apply(1, scale) %>% 
        t %>% 
        magrittr::set_colnames(colnames(matrix))
}

quantile_normalize_matrix <- function(matrix){
    matrix %>% 
        preprocessCore::normalize.quantiles() %>% 
        magrittr::set_colnames(colnames(matrix)) %>% 
        magrittr::set_rownames(rownames(matrix))
}

