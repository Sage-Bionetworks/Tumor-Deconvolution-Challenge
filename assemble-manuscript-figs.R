
library(pacman)
p_load("tools") # for file_ext

dest.fig.dir <- "manuscript-figs"
dir.create(dest.fig.dir)

figs <- 
  list(
       "figs1a" = "presentations/challenge-website/Hierarchical\ Relationship\ of\ Immune\ Cells/V1.3-coarse-grained-aggregated-scaled.svg",
       "figs1b" = "presentations/challenge-website/Hierarchical Relationship of Immune Cells/V1.2-fine-grained-scaled.svg",
       "figs2" = "analysis/summary/purified-samples-marker-heatmap-protein-coding-genes-no-title.png",
       "fig1" = "presentations/deconvolution-schematic.png",
       "fig2" = "analysis/validation-analysis/figs/fig-validation-round-1-performance.png",
       "figs3" = "analysis/validation-analysis/figs/fig-validation-all-performance.png",
       "figs4" = "analysis/validation-analysis/figs/fig-validation-round-1-merged-mixture-distribution-effect.png",
       "fig3" = "analysis/validation-analysis/figs/fig-validation-round-1-strip-and-heatmap-merged-cell-type.png",
       "figs5" = "analysis/validation-analysis/figs/fig-validation-round-1-merged-strip-cell-type.png",
       "figs6" = "analysis/validation-analysis/figs/fig-validation-round-1-coarse-and-fine-strip-cell-type.png",
       "figs7" = "analysis/validation-analysis/figs/fig-validation-heatmap-round-1-coarse-and-fine-cell-type.png",
       "figs8" = "analysis/validation-analysis/figs/fig-validation-round-2-merged-strip-cell-type.png",
       "figs9" = "analysis/validation-analysis/figs/fig-validation-heatmap-rounds-2-and-3-merged-cell-type.png",
       "figs10" = "analysis/validation-analysis/figs/fig-validation-round-2-coarse-and-fine-strip-cell-type.png",
       "figs11" = "analysis/validation-analysis/figs/fig-validation-heatmap-round-2-coarse-and-fine-cell-type.png",
       "figs12" = "analysis/validation-analysis/figs/fig-validation-round-3-merged-strip-cell-type.png",
       "figs13" = "analysis/validation-analysis/figs/fig-validation-round-3-coarse-and-fine-strip-cell-type.png",
       "figs14" = "analysis/validation-analysis/figs/fig-validation-heatmap-round-3-coarse-and-fine-cell-type.png",
       "fig4" = "analysis/sample-level-analysis/figs/sample-level-metric-swarm-round-1.png",
       "fig5" = "analysis/specificity-analysis/figs/spillover-summary.png",
       "figs15" = "analysis/specificity-analysis/figs/spillover-all-scores-coarse-grained-round-1.png",
       "figs16" = "analysis/specificity-analysis/figs/spillover-all-scores-fine-grained-round-1.png",
       "fig6" = "analysis/in-silico-admixtures/figs/sensitivity-spikein-and-summary.png"
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
  extension <- file_ext(fig.file)
  dest.file <- paste0(dest.fig.dir, "/", fig, ".", extension)
  file.copy(fig.file, dest.file, overwrite=TRUE)
}