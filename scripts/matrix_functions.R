require(magrittr)
require(tibble)

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