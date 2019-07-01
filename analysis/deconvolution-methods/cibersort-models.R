suppressPackageStartupMessages(p_load(e1071))
source(cibersort.script)
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(purrr))
suppressPackageStartupMessages(p_load(magrittr))
suppressPackageStartupMessages(p_load(tibble))
suppressPackageStartupMessages(p_load(tidyr))
suppressPackageStartupMessages(p_load(reshape2))

## Execute CIBERSORT against a dataset.
## Assumes that expression_matrix_file points to a CSV whose gene identifiers
## are HUGO symbols.
do_cibersort <- function(expression_matrix_file, dataset_name){
    
    # We are using the HUGO version of the expression file, so this needs to
    # indicate that here. probests and genes are the dataframes created above.
    result_matrix <- CIBERSORT(cibersort.lm22.signature.file, expression_matrix_file, absolute = TRUE, abs_method = "no.sumto1")

    cell.type.cols <- c("B cells naive", "B cells memory", "Plasma cells",
                        "T cells CD8", "T cells CD4 naive",
			"T cells CD4 memory resting", "T cells CD4 memory activated",
			"T cells follicular helper", "T cells regulatory (Tregs)",
			"T cells gamma delta", "NK cells resting", "NK cells activated",
			"Monocytes", "Macrophages M0", "Macrophages M1",
			"Macrophages M2", "Dendritic cells resting",
			"Dendritic cells activated", "Mast cells resting",
			"Mast cells activated", "Eosinophils", "Neutrophils")
    result_matrix <- result_matrix[, cell.type.cols]

    # Convert the result matrix back to a dataframe
    result_df <- result_matrix %>% 
        as.matrix() %>%
	melt()

    colnames(result_df) <- c("sample.id", "cibersort.cell.type", "prediction")
    
    # Add dataset column
    result_df <- dplyr::mutate(result_df, dataset.name = dataset_name)
}

run.cibersort.models <- function(expression_matrix_files, dataset_names) {
  ## Create a table that translates the cell types output by
  ## CIBERSORT ('cibersort.cell.type' column) to the cell types
  ## required of the course-grained sub-Challenge ('cell.type' column).
  fine.grained.translations <-
    list("memory.B.cells" = c("B cells memory"),
         "naive.B.cells" = c("B cells naive"),
         "memory.CD4.T.cells" = c("T cells CD4 memory resting", "T cells CD4 memory activated"),
         "naive.CD4.T.cells" = c("T cells CD4 naive"),
         "regulatory.T.cells" = c("T cells regulatory (Tregs)"),
         "memory.CD8.T.cells" = c("T cells CD8"),
         "naive.CD8.T.cells" = c("T cells CD8"),
         "NK.cells" = c("NK cells resting", "NK cells activated"),
         "neutrophils" = c("Neutrophils"),
         "monocytes" = c("Monocytes"),
         "myeloid.dendritic.cells" = c("Dendritic cells resting", "Dendritic cells activated"),
	 "macrophages" = c("Macrophages M0", "Macrophages M1", "Macrophages M2"),
	 "fibroblasts" = c("non.immune"),
	 "endothelial.cells" = c("non.immune")
        )

  coarse.grained.translations <-
    list("B.cells" = c("B cells naive", "B cells memory"),
         "CD4.T.cells" = c("T cells CD4 naive", "T cells CD4 memory resting", "T cells CD4 memory activated"),
         "CD8.T.cells" = c("T cells CD8"),
         "NK.cells" = c("NK cells resting", "NK cells activated"),
         "neutrophils" = c("Neutrophils"),
         "monocytic.lineage" = c("Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2",
	                         "Dendritic cells resting", "Dendritic cells activated"),
	 "fibroblasts" = c("non.immune"),
	 "endothelial.cells" = c("non.immune")
        )

  ## Run CIBERSORT on each of the expression files
  result_dfs <- purrr::map2(expression_matrix_files, dataset_names, do_cibersort) 

  ## Translate each cell type into coarse-grained results
  coarse_result_dfs <-
    llply(result_dfs,
          .fun = function(df) {
		   dataset <- df$dataset[1]
	           mat <- acast(df[, 1:3], sample.id ~ cibersort.cell.type)
		   vec <- unlist(apply(mat, 1, function(row) 1 - sum(row, na.rm=TRUE)))
		   mat <- cbind(mat, "non.immune" = vec)
		   res <- combine.columns(mat, coarse.grained.translations)
		   new.df <- melt(res)
		   colnames(new.df) <- c("sample.id", "cell.type", "prediction")
		   new.df <- cbind(dataset.name = dataset, new.df)
		   new.df
		 })

  ## Combine all results into one dataframe
  coarse_combined_result_df <- dplyr::bind_rows(coarse_result_dfs)

  ## Translate each cell type into fine-grained results
  fine_result_dfs <-
    llply(result_dfs,
          .fun = function(df) {
		   dataset <- df$dataset[1]
	           mat <- acast(df[, 1:3], sample.id ~ cibersort.cell.type)
		   vec <- unlist(apply(mat, 1, function(row) 1 - sum(row, na.rm=TRUE)))
		   mat <- cbind(mat, "non.immune" = vec)
		   res <- combine.columns(mat, fine.grained.translations)
		   new.df <- melt(res)
		   colnames(new.df) <- c("sample.id", "cell.type", "prediction")
		   new.df <- cbind(dataset.name = dataset, new.df)
		   new.df
		 })
  fine_combined_result_df <- dplyr::bind_rows(fine_result_dfs)
  
  return(list("coarse" = coarse_combined_result_df, "fine" = fine_combined_result_df))
}