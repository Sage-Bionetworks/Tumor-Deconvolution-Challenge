library(plyr)
library(dplyr)
library(ggplot2)

## todo
## - map dataset names
## - fit distributions
## - rank by aurelien and weighting approach
## - how to show rank vs correlations (possibly show top 5 and where they land in distribution)

files <- list("coarse-round1" = "round1-coarse-results.tsv",
              "fine-round1" = "round1-fine-results.tsv")

plot.histograms <- function(tsv, title) {
    tsv <- subset(tsv, is_latest == "true" | grepl(repo_name, pattern = "baseline"))
    tsv <- subset(tsv, metric == "pearson")
    plts <-
        dlply(tsv,
              .variables = c("dataset"),
              .fun = function(df) {
                  g <- ggplot(df)
                  g <- g + geom_histogram(aes(metric_value))
                  g <- g + xlim(c(-1,1))
                  g <- g + facet_wrap(~ celltype, nrow = 2, scale = "free_x")
                  g <- g + xlab("Pearson Correlation (Predicted vs Ground Truth)")
                  g <- g + ggtitle(paste0(title, " ", df$dataset[1], " (n = ", length(unique(df$repo_name)), ")"))
                  g
              })
    plts
}

plot.all.histograms <- function(tsv, title, repo.name.to.highlight) {
    tsv <- subset(tsv, is_latest == "true" | grepl(repo_name, pattern = "baseline"))
    tsv <- subset(tsv, metric == "pearson")
    tsv$facet <- paste(tsv$dataset, tsv$celltype, sep="\n")
    dataInt <- tsv %>%
        group_by(facet) %>%
        filter(repo_name == repo.name.to.highlight)
    g <- ggplot(data = tsv, aes(x = metric_value))
    g <- g + geom_histogram()
    g <- g + xlim(c(-1.2,1.2))
    ##    g <- g + facet_wrap(~ facet, nrow = 3, scale = "free_x")
##    dataInt$metric_value <- as.numeric(dataInt$metric_value)
    g <- g + facet_wrap(~ facet)
    g <- g + geom_vline(data = dataInt, aes(xintercept = metric_value))
    g <- g + xlab("Pearson Correlation (Predicted vs Ground Truth)")
    g <- g + ggtitle(paste0(title, " (n = ", length(unique(tsv$repo_name)), ")"))
    g
}

stop("stop")

nms <- names(files)
for(nm in nms) {
    tsv <- read.table(files[[nm]], sep="\t", header=TRUE, as.is=TRUE)
    pdf(paste0(nm, "-pearson-distributions.pdf"), onefile = TRUE)
    plts <- plot.histograms(tsv, nm)
    for(plt in plts) {
        print(plt)
    }
    d <- dev.off()
}
