suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(cowplot))

suppressPackageStartupMessages(p_load("foreach"))
suppressPackageStartupMessages(p_load("parallel"))
suppressPackageStartupMessages(p_load("reshape2"))

suppressPackageStartupMessages(p_load("patchwork"))
suppressPackageStartupMessages(p_load("data.table"))

# We will apply MCP-Counter to the Challenge data
cat(paste0("Apply MCP-Counter to Challenge validation data\n"))

if(!require(MCPcounter)) {
  suppressPackageStartupMessages(p_load(devtools))
  install_github("ebecht/MCPcounter",ref="master", subdir="Source",upgrade="never")
}
suppressPackageStartupMessages(p_load(MCPcounter))

source("../utils.R")

set.seed(1234)

synLogin()

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(p_load("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

source("validation-analysis-utils.R")

subchallenge.col <- "subchallenge"
measured.col <- "measured"
cell.type.col <- "cell.type"
raw.cell.type.col <- "raw.cell.type"
dataset.name.col <- "dataset.name"
sample.id.col <- "sample.id"
prediction.col <- "prediction"
method.name.col <- "method.name"
round.col <- "submission"

### 1. Load Challenge validation data
cat(paste0("Loading Challenge validataion data\n"))
challenge.expr.mats <- load.challenge.validation.data() 

### 2. Run MCP-Counter against Challenge validation data
cat(paste0("Running MCP-Counter against Challenge validataion data\n"))
# Define a function that encapsulates the deconvolution method.
# It should return a data.frame with columns method.name.col, raw.cell.type.col, sample.id.col, prediction.col 
deconvolve.expr.mat <- 
  function(mat, method.name.col_ = method.name.col, cell.type.col_ = raw.cell.type.col, 
           sample.id.col_ = sample.id.col, prediction.col_ = prediction.col) {
    res <- MCPcounter.estimate(mat, featuresType="HUGO_symbols")
    res <- reshape2::melt(res)
    res <- cbind(method = "my.mcp", res)
    colnames(res) <- c(method.name.col_, cell.type.col_, sample.id.col_, prediction.col_)
    res
  }

my.deconv.res <- ldply(challenge.expr.mats, deconvolve.expr.mat)
colnames(my.deconv.res)[1] <- dataset.name.col

### 3. Translate raw deconvolution cell types into those expected of the 
###    fine- and coarse-grained Challenges
cat(paste0("Translating MCP-Counter cell types\n"))
# Translate raw cell type as output by deconvolution method to
# cell type expected by Challenge
fine_translation_df <- tibble::tribble(
    ~cell.type, ~raw.cell.type,
    "NK.cells", "NK cells",
    "neutrophils", "Neutrophils",
    "myeloid.dendritic.cells", "Myeloid dendritic cells",
    "fibroblasts", "Fibroblasts",
    "endothelial.cells", "Endothelial cells"
)

coarse_translation_df <- tibble::tribble(
    ~cell.type, ~raw.cell.type,
    "B.cells", "B lineage",
    "CD8.T.cells", "CD8 T cells",
    "NK.cells", "NK cells",
    "neutrophils", "Neutrophils",
    "monocytic.lineage", "Monocytic lineage",
    "fibroblasts", "Fibroblasts",
    "endothelial.cells", "Endothelial cells"
)

my.coarse.deconv.res <- aggregate.cell.type.predictions(my.deconv.res, coarse_translation_df)
my.coarse.deconv.res <- cbind(subchallenge = "coarse", my.coarse.deconv.res)

my.fine.deconv.res <- aggregate.cell.type.predictions(my.deconv.res, fine_translation_df)
my.fine.deconv.res <- cbind(subchallenge = "fine", my.fine.deconv.res)

my.deconv.res <- rbind(my.coarse.deconv.res, my.fine.deconv.res)
my.deconv.res <- cbind(my.deconv.res, submission=1)

### 4. Load the Challenge results for the intital submission (=1)
###    (i.e., of comparator and participant deconvolution methods applied to validation data)
cat(paste0("Loading Challenge validation results\n"))
tmp <- load.challenge.validation.results(target.subchallenges = c("coarse", "fine"), target.submissions = c("1"))
challenge.deconv.res <- tmp[[1]]
challenge.dataset.anno <- tmp[[2]]
challenge.ground.truth <- tmp[[3]]

### 5. Compare results computed above for MCP-Counter with those from the Challenge, as an sanity check.
mcp.challenge.deconv.res <- subset(challenge.deconv.res, method.name=="MCP-counter")
cols <- intersect(colnames(mcp.challenge.deconv.res), colnames(my.deconv.res))
cols <- cols[!(cols %in% c("prediction", "method.name"))]
m <- merge(mcp.challenge.deconv.res, my.deconv.res, by = cols)

# Note that, though the results aren't exactly the same, the correlation is very high (>0.99)
# and the _maximum_ relative difference is low (<9%)

rel.diff <- abs(m$prediction.x - m$prediction.y) / pmax(m$prediction.x, m$prediction.y)
cor.p <- cor(m$prediction.x, m$prediction.y)
cat(paste0("Relative difference between re-computed and Challenge MCP-Counter results: ", max(rel.diff), "\n"))
cat(paste0("Pearson correlation between re-computed and Challenge MCP-Counter results: ", cor.p, "\n"))

### 6. Combine MCP-Counter results computed above with Challenge results
cat(paste0("Combinging MCP-Counter and Challenge results\n"))
res.all <- 
  rbind(challenge.deconv.res[, c(dataset.name.col, subchallenge.col, sample.id.col, cell.type.col, prediction.col, method.name.col, round.col)],
        my.deconv.res[, c(dataset.name.col, subchallenge.col, sample.id.col, cell.type.col, prediction.col, method.name.col, round.col)])

### 7. Merge the ground truth into the results, fill in any missing predictions as NA, and rename cell types
cat(paste0("Merging ground truth into results\n"))
res.all <- merge(res.all, challenge.ground.truth)
res.all <- make.missing.predictions.na(res.all)
res.all <- rename.cell.types(res.all, from.col = cell.type.col, to.col = cell.type.col)

### 8. Compute bootstrap statistics over MCP-Counter results computed above and Challenge results.
###    Note that bootstrapped Challenge results are already available at Synapse id syn22951683
cat(paste0("Computing bootstrap statistics\n"))

# Load in the bootstrap samples
synId <- "syn22344963"
obj <- synGet(synId, downloadFile=TRUE)
bootstraps <- readRDS(obj$path)

rds.file <- "boot-res-my-mcp.rds"

# This is long-running ... e.g., ~30minutes on HPC with 32cpus and 500GB mem
if(!file.exists(rds.file)) { 
  boot.res <- do.bootstrap.analysis(res.all, bootstraps, method.name.col, subchallenge.col, measured.col, cell.type.col,
                                    dataset.name.col, sample.id.col, prediction.col, round.col, round = "1")

  saveRDS(boot.res, rds.file)
}

boot.res <- readRDS(rds.file)

### 9. Make plots 
suppressPackageStartupMessages(p_load(xlsx))
method.anno <- get.method.annotations()

plot.spearman.distribution <- TRUE
round <- "1"
postfix <- paste0("-round-", round)

results <- list()
results[[round]] <- boot.res

res.round <- results[[round]][["res.round"]]
res.round <- merge(res.round, challenge.dataset.anno, by = c("dataset.name"))
flag <- res.round[, "distribution.type"] == "Random"
res.round[flag, "distribution.type"] <- "Unconstrained"
bootstrapped.cors <- results[[round]][["bootstrapped.cors"]]
bootstrapped.scores <- results[[round]][["bootstrapped.scores"]]
mean.bootstrapped.scores <- results[[round]][["mean.bootstrapped.scores"]]

## Recalculate bootstrapped score summary using median
# We only use median in the boxplot of pearson (implicitly)
# and in the barplot of spearman (explicitly and for consistency with the pearson boxplot)
summary.fun <- median
median.bootstrapped.scores <-
  llply(bootstrapped.scores, .parallel = TRUE,
        .fun = function(df) {
                 ## Average over bootstraps
                 df <- ddply(df,
                             .variables = c(method.name.col),
                             .fun = function(tmp) {
                                      data.frame(pearson = summary.fun(tmp$pearson), spearman = summary.fun(tmp$spearman), rmse = summary.fun(tmp$rmse), pearson.fc = summary.fun(tmp$pearson.fc))
                              })
                 o <- order(df$pearson, decreasing = TRUE)
                 df <- df[o,]
                 df
             })

means.by.cell.type.method <- results[[round]][["means.by.cell.type.method"]]
means.over.dataset <- results[[round]][["means.over.dataset"]]
top.performers <- results[[round]][["top.performers"]]
bayes.factors <- results[[round]][["bayes.factors"]]

## Need to take annotation for latest round if there isn't one for current round
method.anno.round <- get.round.specific.annotations(method.anno, round)

plots <- list()	   
plots[[round]] <- plot.bootstrap.analysis(res.round, bootstrapped.scores, mean.bootstrapped.scores, median.bootstrapped.scores,
                                means.by.cell.type.method,
                                means.over.dataset, method.anno.round,
                                postfix, plot.spearman.distribution = plot.spearman.distribution)

### 10. Plot heatmap of method x cell type pearson correlation (like Fig 3a)

g.heatmap.merged.round1 <- plots[["1"]][["heatmaps"]][["merged"]]
g.heatmap.merged.round1 <- g.heatmap.merged.round1 + ggtitle("")

png("fig3a-my-mcp.png", width = 2 * 480, height = 1 * 480)
print(g.heatmap.merged.round1)
d <- dev.off()

### 11. Assemble plot like Fig 2 (i.e., distribution over bootstraps of Pearson and Spearman scores for coarse- and fine-grained subchallenges)

suppressPackageStartupMessages(p_load(grid))
suppressPackageStartupMessages(p_load(gridExtra))

g.bootstrap.coarse.pearson.round1 <- plots[["1"]][["barplots"]][["coarse-pearson"]] 
g.bootstrap.coarse.spearman.round1 <- plots[["1"]][["barplots"]][["coarse-spearman"]] 
g.bootstrap.coarse.anno.round1 <- plots[["1"]][["barplots"]][["coarse-anno"]]
g.bootstrap.coarse.anno.legend.round1 <- plots[["1"]][["barplots"]][["coarse-legend"]]

g.bootstrap.coarse.pearson.round1 <- g.bootstrap.coarse.pearson.round1 + ylab("Pearson\nCorrelation") +
    theme(axis.title.y = element_blank())
g.bootstrap.coarse.spearman.round1 <- g.bootstrap.coarse.spearman.round1 + ylab("Spearman\nCorrelation")

title <- "Coarse-Grained (First Submission)"
rel_widths <- c(1,3,3,0.5,0.5)
plot_row <- 
  g.bootstrap.coarse.pearson.round1 + g.bootstrap.coarse.spearman.round1 + g.bootstrap.coarse.anno.round1 + g.bootstrap.coarse.anno.legend.round1 + plot_layout(widths=c(4,4,1,1))
g.bootstrap.coarse.round1 <- plot_grid(textGrob(title, gp = gpar(fontsize = 20)), plot_row, ncol=1, rel_heights = c(0.1, 1))

g.bootstrap.fine.pearson.round1 <- plots[["1"]][["barplots"]][["fine-pearson"]] 
g.bootstrap.fine.spearman.round1 <- plots[["1"]][["barplots"]][["fine-spearman"]] 
g.bootstrap.fine.anno.round1 <- plots[["1"]][["barplots"]][["fine-anno"]]
g.bootstrap.fine.anno.legend.round1 <- plots[["1"]][["barplots"]][["fine-legend"]]

g.bootstrap.fine.pearson.round1 <- g.bootstrap.fine.pearson.round1 + ylab("Pearson\nCorrelation") +
    theme(axis.title.y = element_blank())
g.bootstrap.fine.spearman.round1 <- g.bootstrap.fine.spearman.round1 + ylab("Spearman\nCorrelation")

title <- "Fine-Grained (First Submission)"
rel_widths <- c(1,3,3,0.5,0.5)
plot_row <- 
  g.bootstrap.fine.pearson.round1 + g.bootstrap.fine.spearman.round1 + g.bootstrap.fine.anno.round1 + g.bootstrap.fine.anno.legend.round1 + plot_layout(widths=c(4,4,1,1))
g.bootstrap.fine.round1 <- plot_grid(textGrob(title, gp = gpar(fontsize = 20)), plot_row, ncol=1, rel_heights = c(0.1, 1))

g <- plot_grid(g.bootstrap.coarse.round1, g.bootstrap.fine.round1, labels = c("A", "B"))

png("fig2-my-mcp.png", width = 2.5 * 480, height = 1 * 480)
print(g)
d <- dev.off()

cat("Exiting successfully\n")
q(status=0)
