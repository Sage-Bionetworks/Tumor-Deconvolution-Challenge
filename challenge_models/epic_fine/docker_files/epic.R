
###############################################################################
## This code uses Epic to calculate cell type predictions for
## the fine-grained sub-Challenge.
###############################################################################

library(EPIC)
library(readr)
library(tibble)
library(magrittr)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)


## Read in the round and sub-Challenge-specific input file 
## listing each of the datasets
print(list.files())
print(getwd())
input_df <- readr::read_csv("input/input.csv")

## Extract the names of each dataset
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

## Get cancer types 
cancer_types <- input_df$cancer.type

# EPIC to challenge cell type names
translation_df <- tibble::tribble(
    ~cell.type, ~epic.cell.type,
    "NK.cells", "NKcells",
    "macrophages", "Macrophages",
    "fibroblasts", "CAFs",
    "endothelial.cells", "Endothelial",
    "monocytes", "Monocytes",
    "neutrophils", "Neutrophils"
)


## Assumes that expression_path points to a CSV whose gene identifiers
## are HUGO symbols.
do_epic <- function(
    expression_path, 
    dataset_name, 
    scale, 
    normalization, 
    cancer_type
){
    
    # normalization must be one of these methods
    if (!normalization %in% allowed_normalization_methods) {
        stop(paste0("Non-accepted normalization method: ", normalization))
    }
    
    expression_matrix <- expression_path %>% 
        readr::read_csv() %>% 
        as.data.frame() %>% 
        tibble::column_to_rownames(., colnames(.)[[1]]) %>% 
        as.matrix()
    
    if (scale == "Linear") {
        expression_matrix <- expression_matrix
    } else if (scale == "Log2") {
        expression_matrix <- 2^expression_matrix
    } else if (scale == "Log10") {
        expression_matrix <- 10^expression_matrix
    } else {
        stop(paste0("non-accepted scale method: ", scale))
    }
    
    if (is.na(cancer_type)) {
        reference <- "TRef"
    } else if (cancer_type == "") {
        reference <- "TRef"
    } else {
        reference <- "BRef"
    }
    
    expression_matrix %>% 
        EPIC::EPIC(reference = reference, mRNA_cell = FALSE) %>% 
        magrittr::use_series(mRNAProportions) %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("sample.id") %>% 
        dplyr::as_tibble() %>% 
        tidyr::gather(
            key = epic.cell.type, 
            value = prediction, 
            -sample.id) %>% 
        dplyr::mutate(dataset.name = dataset_name)
}

## Run EPIC on each of the expression files
result_dfs <- purrr::pmap(
    list(expression_paths, dataset_names, scales, normalizations, cancer_types),
    do_epic
) 

## Combine all results into one dataframe
combined_result_df <- dplyr::bind_rows(result_dfs)

## This will be dependent on which reference we use
## tmp <- as.data.frame(translation_df)
## col <- "epic.cell.type"
## flag <- !(tmp[, col] %in% as.data.frame(combined_result_df)[,col])
## if(any(flag)) {
##    missed.cell.types <- unique(tmp[flag,col])
##    stop("Method did not returned a cell type expected by the translation: ",
##         paste0(missed.cell.types, collapse = ", "), "\n")
## }

## Translate cell type names as output from MCP-Counter to those
## required for the coarse-grained sub-Challenge.
combined_result_df <- combined_result_df %>% 
    dplyr::inner_join(translation_df) %>% 
    dplyr::select(dataset.name, sample.id, cell.type, prediction) %>% 
    dplyr::group_by(dataset.name, sample.id, cell.type) %>% 
    dplyr::summarise(prediction = sum(prediction)) %>% 
    dplyr::ungroup()

##### MCPcounter example code above ########

## Create the directory the output will go into
dir.create("output")

## Write result into output directory
readr::write_csv(combined_result_df, "output/predictions.csv")

print(list.files("output"))
    
