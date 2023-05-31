
library(pacman)
p_load("tools") # for file_ext

dest.fig.dir <- "manuscript-figs"
dir.create(dest.fig.dir)

figs <- 
  list(
#       "figs1a.svg" = "presentations/challenge-website/Hierarchical\ Relationship\ of\ Immune\ Cells/V1.3-coarse-grained-aggregated-scaled.svg",
#       "figs1b.svg" = "presentations/challenge-website/Hierarchical Relationship of Immune Cells/V1.2-fine-grained-scaled.svg",
       "figs1.png" = "analysis/summary/purified-samples-marker-heatmap-protein-coding-genes-no-title.png",
       "figs1.pdf" = "analysis/summary/purified-samples-marker-heatmap-protein-coding-genes-no-title.pdf",
       "fig1a.png" = "external-figs-and-tables/V4.png",
       "fig1b.png" = "presentations/deconvolution-schematic.png",
       "fig2.png" = "analysis/validation-analysis/figs/fig-validation-round-1-performance.png",
       "fig2.pdf" = "analysis/validation-analysis/figs/fig-validation-round-1-performance.pdf",
       "figs2.png" = "external-figs-and-tables/da_505-supp-fig1.png",
       "figs3.png" = "external-figs-and-tables/Biogem_ComparisonDataset.png",
       "figs4.png" = "analysis/validation-analysis/figs/fig-validation-performance-across-rounds.png",
       "figs4.pdf" = "analysis/validation-analysis/figs/fig-validation-performance-across-rounds.pdf",
       "figs5.png" = "analysis/validation-analysis/figs/fig-validation-all-performance.png",
       "figs5.pdf" = "analysis/validation-analysis/figs/fig-validation-all-performance.pdf",
       #"figs5.png" = "analysis/validation-analysis/figs/fig-validation-round-1-merged-mixture-distribution-effect.png",
       #"figs5.pdf" = "analysis/validation-analysis/figs/fig-validation-round-1-merged-mixture-distribution-effect.pdf",
       "fig3.png" = "analysis/validation-analysis/figs/fig-validation-round-1-strip-and-heatmap-merged-cell-type.png",
       "fig3.pdf" = "analysis/validation-analysis/figs/fig-validation-round-1-strip-and-heatmap-merged-cell-type.pdf",
       "figs6.png" = "analysis/validation-analysis/figs/fig-validation-round-1-merged-strip-cell-type.png",
       "figs6.pdf" = "analysis/validation-analysis/figs/fig-validation-round-1-merged-strip-cell-type.pdf",
       "figs7.png" = "analysis/validation-analysis/figs/fig-validation-round-1-coarse-and-fine-strip-cell-type.png",
       "figs7.pdf" = "analysis/validation-analysis/figs/fig-validation-round-1-coarse-and-fine-strip-cell-type.pdf",
       "figs8.png" = "analysis/validation-analysis/figs/fig-validation-heatmap-round-1-coarse-and-fine-cell-type.png",
       "figs8.pdf" = "analysis/validation-analysis/figs/fig-validation-heatmap-round-1-coarse-and-fine-cell-type.pdf",
       "figs9.png" = "analysis/validation-analysis/figs/fig-validation-round-2-merged-strip-cell-type.png",
       "figs9.pdf" = "analysis/validation-analysis/figs/fig-validation-round-2-merged-strip-cell-type.pdf",
       "figs10.png" = "analysis/validation-analysis/figs/fig-validation-heatmap-rounds-2-and-3-merged-cell-type.png",
       "figs10.pdf" = "analysis/validation-analysis/figs/fig-validation-heatmap-rounds-2-and-3-merged-cell-type.pdf",
       "figs11.png" = "analysis/validation-analysis/figs/fig-validation-round-2-coarse-and-fine-strip-cell-type.png",
       "figs11.pdf" = "analysis/validation-analysis/figs/fig-validation-round-2-coarse-and-fine-strip-cell-type.pdf",
       "figs12.png" = "analysis/validation-analysis/figs/fig-validation-heatmap-round-2-coarse-and-fine-cell-type.png",
       "figs12.pdf" = "analysis/validation-analysis/figs/fig-validation-heatmap-round-2-coarse-and-fine-cell-type.pdf",
       "figs13.png" = "analysis/validation-analysis/figs/fig-validation-round-3-merged-strip-cell-type.png",
       "figs13.pdf" = "analysis/validation-analysis/figs/fig-validation-round-3-merged-strip-cell-type.pdf",
       "figs14.png" = "analysis/validation-analysis/figs/fig-validation-round-3-coarse-and-fine-strip-cell-type.png",
       "figs14.pdf" = "analysis/validation-analysis/figs/fig-validation-round-3-coarse-and-fine-strip-cell-type.pdf",
       "figs15.png" = "analysis/validation-analysis/figs/fig-validation-heatmap-round-3-coarse-and-fine-cell-type.png",
       "figs15.pdf" = "analysis/validation-analysis/figs/fig-validation-heatmap-round-3-coarse-and-fine-cell-type.pdf",
       "fig4.png" = "analysis/sample-level-analysis/figs/sample-level-metric-swarm-round-1.png",
       "fig4.pdf" = "analysis/sample-level-analysis/figs/sample-level-metric-swarm-round-1.pdf",
       "fig5.png" = "analysis/specificity-analysis/figs/spillover-summary.png",
       "fig5.pdf" = "analysis/specificity-analysis/figs/spillover-summary.pdf",
       "figs16.png" = "analysis/specificity-analysis/figs/spillover-all-scores-coarse-grained-round-1.png",
       "figs16.pdf" = "analysis/specificity-analysis/figs/spillover-all-scores-coarse-grained-round-1.pdf",
       "figs17.png" = "analysis/specificity-analysis/figs/spillover-all-scores-fine-grained-round-1.png",
       "figs17.pdf" = "analysis/specificity-analysis/figs/spillover-all-scores-fine-grained-round-1.pdf",
       "fig6.png" = "analysis/in-silico-admixtures/figs/sensitivity-spikein-and-summary.png",
       "fig6.pdf" = "analysis/in-silico-admixtures/figs/sensitivity-spikein-and-summary.pdf",
       "figs18.png" = "analysis/cancer-validation/figs/fig-cancer-validation-heatmap-wu-and-pelka.png",
       "figs18.pdf" = "analysis/cancer-validation/figs/fig-cancer-validation-heatmap-wu-and-pelka.pdf",
       "figs19.png" = "analysis/cancer-validation/figs/fig-cancer-validation-heatmap-all.png",
       "figs19.pdf" = "analysis/cancer-validation/figs/fig-cancer-validation-heatmap-all.pdf",
       "fig7.png" = "analysis/cancer-validation/figs/fig-cancer-validation.png",
       "fig7.pdf" = "analysis/cancer-validation/figs/fig-cancer-validation.pdf",
       "figs20.png" = "analysis/cancer-validation/figs/fig-cancer-validation-per-cell-type.png",
       "figs20.pdf" = "analysis/cancer-validation/figs/fig-cancer-validation-per-cell-type.pdf",
       "fig8.png" = "analysis/summary-table/figs/binned-score-heatmap.png",
       "fig8.pdf" = "analysis/summary-table/figs/binned-score-heatmap.pdf"

)

tbls <-
  list(
       "tables10.csv" = "training-sample-metadata/geo-expression-array-immune-cells.csv",
       "tables11.csv" = "training-sample-metadata/geo-rnaseq-immune-cells.csv",
       "tables12.tsv" = "analysis/timing/method-run-times.tsv",
       "tables19.tsv" = "analysis/cancer-validation/figs/cancer-validation-correlations.tsv",
       "tables20.tsv" = "analysis/cancer-validation/figs/cancer-validation-dataset-pvals.tsv"       
)

for(fig in names(figs)) {
  fig.file <- figs[[fig]]
  if(file.exists(fig.file)) {
    cat(paste0("Found ", fig, ": ", fig.file, "\n"))
  } else {
    stop(paste0("Could not find ", fig, ": ", fig.file, "\n"))
  }
}

for(fig in names(figs)) {
  fig.file <- figs[[fig]]
  #extension <- file_ext(fig.file)
  #dest.file <- paste0(dest.fig.dir, "/", fig, ".", extension)
  dest.file <- paste0(dest.fig.dir, "/", fig)
  file.copy(fig.file, dest.file, overwrite=TRUE)
}

for(tbl in names(tbls)) {
  tbl.file <- tbls[[tbl]]
  if(file.exists(tbl.file)) {
    cat(paste0("Found ", tbl, ": ", tbl.file, "\n"))
  } else {
    stop(paste0("Could not find ", tbl, ": ", tbl.file, "\n"))
  }
}

for(tbl in names(tbls)) {
  tbl.file <- tbls[[tbl]]
  #extension <- file_ext(tbl.file)
  #dest.file <- paste0(dest.fig.dir, "/", tbl, ".", extension)
  dest.file <- paste0(dest.fig.dir, "/", tbl)
  file.copy(tbl.file, dest.file, overwrite=TRUE)
}
