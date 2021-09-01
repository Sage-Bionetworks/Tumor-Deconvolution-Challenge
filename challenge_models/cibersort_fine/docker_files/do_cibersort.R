
###############################################################################
## This code uses Cibersort to calculate cell type predictions for
## the fine-grained sub-Challenge.
###############################################################################


library(dplyr)
library(readr)
library(purrr)
library(magrittr)
library(tibble)
library(tidyr)
library(e1071)
library(parallel)
library(preprocessCore)

source("CIBERSORT.R")


## Read in the round and sub-Challenge-specific input file 
## listing each of the datasets
print(list.files())
print(getwd())
input_df <- readr::read_csv("input/input.csv")

# convert 

## Extract the names of each datasetit
dataset_names <- input_df$dataset.name

## Extract the names of the expression files that use 
## Hugo symbols as gene identifiers
expression_files  <- input_df$hugo.expr.file

## Form the paths of the expression files
expression_paths <- paste0("input/", expression_files)

## get the scale methods from the input file
scales <- input_df$scale

## get the scale methods from the input file
normalizations <- input_df$normalization

allowed_normalization_methods <- c(
    "CPM", "MAS5", "gcRMA", "RMA", "RMA+quantile normalization+FARMS", 
    "average", "TMM", "RMA+quantile normalization", "normexp", "TPM", "fRMA"
)

# Cibersort translation table
translation_df <- tibble::tribble(
    ~cell.type, ~cibersort.cell.type,
    "memory.B.cells", "B cells memory",
    "naive.B.cells", "B cells naive",
    "memory.CD4.T.cells", "T cells CD4 memory activated",
    "memory.CD4.T.cells", "T cells CD4 memory resting",
    "naive.CD4.T.cells", "T cells CD4 naive",
    "regulatory.T.cells", "T cells regulatory (Tregs)",
    "NK.cells", "NK cells resting",
    "NK.cells", "NK cells activated",
    "neutrophils", "Neutrophils",
    "monocytes", "Monocytes",
    "myeloid.dendritic.cells", "Dendritic cells resting",
    "myeloid.dendritic.cells", "Dendritic cells activated",
    "macrophages", "Macrophages M0",
    "macrophages", "Macrophages M1",
    "macrophages", "Macrophages M2"
)

## Execute Cibersort against a dataset.
## Assumes that expression_path points to a CSV whose gene identifiers
## are HUGO symbols.
do_cibersort <- function(expression_path, dataset_name, scale, normalization){
    
    # normalization must be one of these methods
    if (!normalization %in% allowed_normalization_methods) {
        stop(paste0("Non-accepted normalization method: ", normalization))
    }
    
    # files must be in tsv format with no missing values
    expression_df <- expression_path %>% 
        readr::read_csv() %>% 
        tidyr::drop_na()
    
    if (scale == "Linear") {
        expression_df <- expression_df
    } else if (scale == "Log2") {
        expression_df <- expression_df %>% 
            tibble::column_to_rownames("Gene") %>% 
            as.matrix() %>% 
            2^. %>% 
            as.data.frame() %>% 
            tibble::rownames_to_column("Gene")
    } else if (scale == "Log10") {
        expression_df <- expression_df %>% 
            tibble::column_to_rownames("Gene") %>% 
            as.matrix() %>% 
            10^. %>% 
            as.data.frame() %>% 
            tibble::rownames_to_column("Gene")
    } else {
        stop("non-accepted scale method")
    }
    
    readr::write_tsv(expression_df, "expr.tsv")
    
    result_matrix <- CIBERSORT(
        "LM22.tsv", 
        "expr.tsv", 
        abs_method = "sig.score",
        absmean = TRUE,
        QN = FALSE
    )
    
    # remove expresion file
    file.remove("expr.tsv")
    
    # Convert the result matrix back to a dataframe, remove non cell type, 
    # columns, create non-immune columns
    result_df <- result_matrix %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("sample.id") %>% 
        dplyr::as_tibble() %>% 
        dplyr::select(-c(
            "P-value", 
            "Correlation", 
            "RMSE", 
            dplyr::starts_with("Absolute score")
        )) %>% 
        dplyr::mutate(non.immune = 1 - rowSums(.[-1])) 
        
    
    # Stack the predictions into one column
    result_df <- tidyr::gather(
        result_df,
        key = cibersort.cell.type, 
        value = prediction, 
        -sample.id) 
    
    # Add dataset column
    result_df <- dplyr::mutate(result_df, dataset.name = dataset_name)
}

## Run Cibersort on each of the expression files
result_dfs <- purrr::pmap(
    list(expression_paths, dataset_names, scales, normalizations),
    do_cibersort
) 

## Combine all results into one dataframe
combined_result_df <- dplyr::bind_rows(result_dfs)

tmp <- as.data.frame(translation_df)
col <- "cibersort.cell.type"
flag <- !(tmp[, col] %in% as.data.frame(combined_result_df)[,col])
if(any(flag)) {
    missed.cell.types <- unique(tmp[flag,col])
    stop("Method did not returned a cell type expected by the translation: ",
         paste0(missed.cell.types, collapse = ", "), "\n")
}

## Translate cell type names as output from MCP-Counter to those
## required for the coarse-grained sub-Challenge.
combined_result_df <- combined_result_df %>% 
    dplyr::inner_join(translation_df) %>% 
    dplyr::select(dataset.name, sample.id, cell.type, prediction) %>% 
    dplyr::group_by(dataset.name, sample.id, cell.type) %>% 
    dplyr::summarise_all(mean)

# The last line above has a bug -- 'mean' should be 'sum'.
# However, to avoid having to re-run this through all of the challenge
# infastruture, I will instead manually correct the results downstream.
# Note that the difference between a mean of n cell types output by
# CIBERSORT and aggregated into a single challenge cell type is simply
# 1/n * the sum. Hence, knowing the number of CIBERSORT cell types
# that are aggregated, we can simply scale the results accordingly.
# Rather that determining this logically, I will do it empirically
# by comparing output derived from CIBERSORTx using similarly
# buggy code and that from corrected code.

## Create the directory the output will go into
dir.create("output")

## Write result into output directory
readr::write_csv(combined_result_df, "output/predictions.csv")

print(list.files("output"))
    
