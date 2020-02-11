
###############################################################################
## This code uses timer to calculate cell type predictions for
## the coarse-grained sub-Challenge.
###############################################################################


library(magrittr)
library(immunedeconv)

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


# timer to challenge cell type names
translation_df <- tibble::tribble(
    ~cell.type, ~quantiseq.cell.type,
    "B.cells", "B_cell",
    "CD4.T.cells", "T_cell.CD4",
    "CD8.T.cells", "T_cell.CD8",
    "neutrophils", "Neutrophil",
    "monocytic.lineage", "Macrophage",
    "monocytic.lineage", "DC"
)

## get the scale methods from the input file
scales <- input_df$scale

## get the scale methods from the input file
normalizations <- input_df$normalization

allowed_normalization_methods <- c(
    "CPM", "MAS5", "gcRMA", "RMA", "RMA+quantile normalization+FARMS", 
    "average", "TMM", "RMA+quantile normalization", "normexp", "TPM"
)

## Get cancer types 
cancer_types <- input_df$cancer.type

## Assumes that expression_path points to a CSV whose gene identifiers
## are HUGO symbols.
do_timer <- function(
    expression_path, 
    dataset_name,
    scale, 
    normalization, 
    cancer_type
){
    
    # normalization must be one of these methods
    if (!normalization %in% allowed_normalization_methods) {
        stop("non-accepted normalization method")
    }
    
    if (cancer_type == "BRCA") {
        indictations <- "brca" 
        reference <- "TRef"
    } else if (cancer_type == "CRC") {
        indications <- "coad"
        reference <- "TRef"
    } else if (cancer_type == "FL") {
        indications <- "dlbc"
        reference <- "TRef"
    } else if (is.na(cancer_type)) {
        indications <- "coad"
        reference <- "BRef"
    } else if (cancer_type == "") {
        indications <- "coad"
        reference <- "BRef"
    } else {
        stop("Unallowed cancer type")
    }
    
    expression_matrix <- expression_path %>%
        readr::read_csv() %>%
        tidyr::drop_na() %>% 
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
        stop("non-accepted scale method")
    }
    
    expression_df <- expression_matrix %>% 
        immunedeconv::deconvolute(
            method = "timer",
            indications = indications,
            reference = reference
        ) %>% 
        dplyr::as_tibble() %>%
        dplyr::rename(quantiseq.cell.type = cell_type) %>% 
        tidyr::gather(
            key = sample.id,
            value = prediction,
            -quantiseq.cell.type
        ) %>% 
        print()
    
    dplyr::tibble(
        "quantiseq.cell.type" = c("Fibroblasts", "Endothelial cells"),
        "prediction" = 0
    ) %>% 
        merge(df$sample.id) %>%
        dplyr::distinct() %>% 
        dplyr::rename(sample.id = y) %>% 
        dplyr::bind_rows(df) %>% 
        dplyr::mutate(dataset.name = dataset_name)
    
    
}

## Run Timer on each of the expression files
result_dfs <- purrr::pmap(
    list(expression_paths, dataset_names, scales, 
         normalizations, cancer_types),
    do_timer
)

## Combine all results into one dataframe
combined_result_df <- dplyr::bind_rows(result_dfs)

## Translate cell type names as output from Quantiseq to those
## required for the coarse-grained sub-Challenge.
combined_result_df2 <- combined_result_df %>% 
    dplyr::inner_join(translation_df) %>% 
    dplyr::select(dataset.name, sample.id, cell.type, prediction) %>% 
    dplyr::group_by(dataset.name, sample.id, cell.type) %>% 
    dplyr::summarise(prediction = sum(prediction)) %>% 
    dplyr::ungroup()


## Create the directory the output will go into
dir.create("output")

## Write result into output directory
readr::write_csv(combined_result_df2, "output/predictions.csv")

print(list.files("output"))
