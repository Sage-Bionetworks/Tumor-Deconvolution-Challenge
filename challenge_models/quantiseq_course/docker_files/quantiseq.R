
###############################################################################
## This code uses quantiseq to calculate cell type predictions for
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


# quantiseq to challenge cell type names
translation_df <- tibble::tribble(
    ~cell.type,          ~quantiseq.cell.type,
    "B.cells",           "B.cells", 
    "CD4.T.cells",       "T.cells.CD4", 
    "CD4.T.cells",       "Tregs", 
    "CD8.T.cells",       "T.cells.CD8", 
    "NK.cells",          "NK.cells", 
    "neutrophils",       "Neutrophils", 
    "monocytic.lineage", "Macrophages.M1", 
    "monocytic.lineage", "Macrophages.M2", 
    "monocytic.lineage", "Monocytes", 
    "monocytic.lineage", "Dendritic.cells"
)
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

## Get platforms
platforms <- input_df$platform

## Assumes that expression_path points to a CSV whose gene identifiers
## are HUGO symbols.
do_quantiseq <- function(
    expression_path, 
    dataset_name,
    scale, 
    normalization, 
    cancer_type,
    platform
){
    if (platform %in% c(
        "Illumina HiSeq 2000", "Illumina NovaSeq", "Illumina", 
        "Illumina HiSeq 4000", "Illumina NextSeq 500")) {
        arrays <- FALSE
    } else if (platform %in% c(
        "Affymetrix HG-U133 Plus 2.0", "Affymetrix Human Gene 1.0 ST",
        "Illumina HumanHT-12 V4.0", "Affymetrix Human Gene 1.1 ST",
        "Affymetrix Human Gene PrimeView", "Illumina HumanHT-12 V4.0",
        "Affymetrix HG-U133A")) {
        arrays <- TRUE
    } else {
        stop(paste0("platform not allowed:", platform))
    }
    
    # normalization must be one of these methods
    if (!normalization %in% allowed_normalization_methods) {
        stop(paste0("Non-accepted normalization method: ", normalization))
    }
    
    if (is.na(cancer_type)) {
        tumor <- FALSE
    } else if (cancer_type == "") {
        tumor <- FALSE
    } else {
        tumor <- TRUE 
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
        stop(paste0("non-accepted scale method: ", scale))
    }
    
    tbl <- 
        immunedeconv::deconvolute_quantiseq(
            expression_matrix,
            tumor = tumor,
            arrays = arrays,
            scale_mrna = FALSE
        ) %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("cell.type") %>% 
        dplyr::as_tibble() %>%
        tidyr::gather(
            key = sample.id, 
            value = prediction,
            -cell.type
        ) %>% 
        dplyr::mutate(dataset.name = dataset_name)
}

## Run Quantiseq on each of the expression files
result_dfs <- purrr::pmap(
    list(expression_paths, dataset_names, scales, normalizations, cancer_types, platforms),
    do_quantiseq
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
