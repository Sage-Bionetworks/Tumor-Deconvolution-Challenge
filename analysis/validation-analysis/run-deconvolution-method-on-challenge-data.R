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

if(!require(xCell)) {
  suppressPackageStartupMessages(p_load(devtools))
  devtools::install_github('dviraran/xCell', upgrade = "never")
}
cat(paste0("Apply xCell to Challenge validation data\n"))
suppressPackageStartupMessages(p_load(xCell))

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

### 2. Run xCell against Challenge validation data
cat(paste0("Running xCell against Challenge validataion data\n"))
# Define a function that encapsulates the deconvolution method.
# It should return a data.frame with columns method.name.col, raw.cell.type.col, sample.id.col, prediction.col 
deconvolve.expr.mat <- 
  function(mat, cell.types, method.name.col_ = method.name.col, cell.type.col_ = raw.cell.type.col, 
           sample.id.col_ = sample.id.col, prediction.col_ = prediction.col) {
    res <- xCell::xCellAnalysis(mat, rnaseq = TRUE, cell.types.use = cell.types)
    res <- reshape2::melt(res)
    res <- cbind(method = "my.xcell", res)
    colnames(res) <- c(method.name.col_, cell.type.col_, sample.id.col_, prediction.col_)
    res
  }

### 3. Translate raw deconvolution cell types into those expected of the 
###    fine- and coarse-grained Challenges
cat(paste0("Translating MCP-Counter cell types\n"))
# Translate raw cell type as output by deconvolution method to
# cell type expected by Challenge
fine_translation_df <- tibble::tribble(
    ~cell.type, ~raw.cell.type,
    "memory.B.cells", "Memory B-cells",
    "naive.B.cells", "naive B-cells",
    "memory.CD4.T.cells", "CD4+ memory T-cells",
    "naive.CD4.T.cells", "CD4+ naive T-cells",
    "regulatory.T.cells", "Tregs",
    "memory.CD8.T.cells", "CD8+ Tem",
    "naive.CD8.T.cells", "CD8+ naive T-cells",
    "NK.cells", "NK cells",
    "neutrophils", "Neutrophils",
    "monocytes", "Monocytes",
    "myeloid.dendritic.cells", "DC",
    "macrophages", "Macrophages",
    "fibroblasts", "Fibroblasts",
    "endothelial.cells", "Endothelial cells"
)

coarse_translation_df <- tibble::tribble(
    ~cell.type, ~raw.cell.type,
    "B.cells", "B-cells",
    "CD4.T.cells", "CD4+ T-cells",
    "CD8.T.cells", "CD8+ T-cells",
    "NK.cells", "NK cells",
    "neutrophils", "Neutrophils",
    "monocytic.lineage", "Monocytes",
    "fibroblasts", "Fibroblasts",
    "endothelial.cells","Endothelial cells"
)

my.coarse.deconv.res <- ldply(challenge.expr.mats, .fun = function(mat) deconvolve.expr.mat(mat, coarse_translation_df$raw.cell.type))
colnames(my.coarse.deconv.res)[1] <- dataset.name.col

my.fine.deconv.res <- ldply(challenge.expr.mats, .fun = function(mat) deconvolve.expr.mat(mat, fine_translation_df$raw.cell.type))
colnames(my.fine.deconv.res)[1] <- dataset.name.col

my.coarse.deconv.res <- aggregate.cell.type.predictions(my.coarse.deconv.res, coarse_translation_df)
my.coarse.deconv.res <- cbind(subchallenge = "coarse", my.coarse.deconv.res)

my.fine.deconv.res <- aggregate.cell.type.predictions(my.fine.deconv.res, fine_translation_df)
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

### 5. Compare results computed above for xCell with those from the Challenge, as an sanity check.
xcell.challenge.deconv.res <- subset(challenge.deconv.res, method.name=="xCell")
cols <- intersect(colnames(xcell.challenge.deconv.res), colnames(my.deconv.res))
cols <- cols[!(cols %in% c("prediction", "method.name"))]
m <- merge(xcell.challenge.deconv.res, my.deconv.res, by = cols)

# Note that, though the results aren't exactly the same, the correlation is very high (>0.99)
# and the _maximum_ relative difference is low (<9%)

eps <- 10^-6
rel.diff <- abs(m$prediction.x - m$prediction.y) / pmax(abs(m$prediction.x), abs(m$prediction.y), eps)
cor.p <- cor(m$prediction.x, m$prediction.y)
cat(paste0("Relative difference between re-computed and Challenge xCell results: ", max(rel.diff), "\n"))
cat(paste0("Pearson correlation between re-computed and Challenge xCell results: ", cor.p, "\n"))

### 6. Combine xCell results computed above with Challenge results
cat(paste0("Combinging xCell and Challenge results\n"))
res.all <- 
  rbind(challenge.deconv.res[, c(dataset.name.col, subchallenge.col, sample.id.col, cell.type.col, prediction.col, method.name.col, round.col)],
        my.deconv.res[, c(dataset.name.col, subchallenge.col, sample.id.col, cell.type.col, prediction.col, method.name.col, round.col)])

### 7. Merge the ground truth into the results, fill in any missing predictions as NA, and rename cell types
cat(paste0("Merging ground truth into results\n"))
res.all <- merge(res.all, challenge.ground.truth)
res.all <- make.missing.predictions.na(res.all)
res.all <- rename.cell.types(res.all, from.col = cell.type.col, to.col = cell.type.col)

### 8. Compute bootstrap statistics over xCell results computed above and Challenge results.
###    Note that bootstrapped Challenge results are already available at Synapse id syn22951683
cat(paste0("Computing bootstrap statistics\n"))

# Load in the bootstrap samples
synId <- "syn22344963"
obj <- synGet(synId, downloadFile=TRUE)
bootstraps <- readRDS(obj$path)

rds.file <- "boot-res-my-xcell.rds"

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

png("fig3a-my-xcell.png", width = 2 * 480, height = 1 * 480)
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

png("fig2-my-xcell.png", width = 2.5 * 480, height = 1 * 480)
print(g)
d <- dev.off()

# Plot the pearson fold changes as well

g.bootstrap.coarse.pearson.fc.round1 <- plots[["1"]][["barplots"]][["coarse-pearson.fc"]] 
g.bootstrap.coarse.pearson.fc.round1 <- g.bootstrap.coarse.pearson.fc.round1 + ylab("Fold Change\n(Pearson)") +
    theme(axis.title.y = element_blank())
g.bootstrap.fine.pearson.fc.round1 <- plots[["1"]][["barplots"]][["fine-pearson.fc"]] 
g.bootstrap.fine.pearson.fc.round1 <- g.bootstrap.fine.pearson.fc.round1 + ylab("Fold Change\n(Pearson)") +
    theme(axis.title.y = element_blank())

g.bootstrap.coarse.pearson.fc.round1 <- g.bootstrap.coarse.pearson.fc.round1 + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
g.bootstrap.fine.pearson.fc.round1 <- g.bootstrap.fine.pearson.fc.round1 + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())

title <- "Coarse-Grained (First Submission)"
plot_row <- 
  g.bootstrap.coarse.pearson.round1 + g.bootstrap.coarse.spearman.round1 + g.bootstrap.coarse.pearson.fc.round1 + g.bootstrap.coarse.anno.round1 + g.bootstrap.coarse.anno.legend.round1 + plot_layout(widths=c(4,4,4,1,1))
g.bootstrap.coarse.round1 <- plot_grid(textGrob(title, gp = gpar(fontsize = 20)), plot_row, ncol=1, rel_heights = c(0.1, 1))

title <- "Fine-Grained (First Submission)"
plot_row <- 
  g.bootstrap.fine.pearson.round1 + g.bootstrap.fine.spearman.round1 + g.bootstrap.fine.pearson.fc.round1 + g.bootstrap.fine.anno.round1 + g.bootstrap.fine.anno.legend.round1 + plot_layout(widths=c(4,4,4,1,1))
g.bootstrap.fine.round1 <- plot_grid(textGrob(title, gp = gpar(fontsize = 20)), plot_row, ncol=1, rel_heights = c(0.1, 1))

g <- plot_grid(g.bootstrap.coarse.round1, g.bootstrap.fine.round1, labels = c("A", "B"))

png("fig2-my-xcell-with-fold-change.png", width = 3 * 480, height = 1 * 480)
print(g)
d <- dev.off()


cat("Exiting successfully\n")
q(status=0)

head(subset(boot.res[["bootstrapped.cors"]][["coarse"]], method.name=="Aginome-XMU" & cell.type == "B"))

# Aginome has some negative values -- what is the fold change across positive and negative?
ds.name <- "AA"
boot <- 1
ct <- "B"
mn <- "Aginome-XMU"
sb <- subset(boot.res[["res.round"]], dataset.name == ds.name & method.name == mn & cell.type == ct & submission == 1)

mn <- "CIBERSORTx"
sb <- subset(boot.res[["res.round"]], dataset.name == ds.name & method.name == mn & cell.type == ct & submission == 1)

# CIBERSORTx has all positive values, but the correlation is still poor.
# Note that the measured values aren't spaced very far apart -- as we would have (e.g., by decades)
# to see fold change.

df <- sb


pred <- sb$prediction
measured <- sb$measured
o <- order(measured, decreasing=FALSE)
eps <- 10^-5
gt.ordered <- measured[o] + eps
# NB: the ordering is established by the _ground truth_ values and applied to those
# values and the predicted values
pred.ordered <- pred[o] + eps
gt.fc <- unlist(lapply(2:length(gt.ordered), function(indx) gt.ordered[indx]/gt.ordered[indx-1]))
pred.fc <- unlist(lapply(2:length(pred.ordered), function(indx) pred.ordered[indx]/pred.ordered[indx-1]))
val <- cor(pred.fc, gt.fc, method = "pearson")

df <- data.frame(prediction = pred.fc, measured = gt.fc)

library(ggpubr) 

# Draw the scatterplot of measured vs predicted and of measured fold change vs predicted fold change
g1 <- ggplot(sb, aes(x=measured, y=prediction)) + geom_point() + stat_cor(method = "pearson", label.x = min(sb$measured)*1.1, label.y = max(sb$prediction)*0.9, size=10)
g1 <- g1 + theme(text = element_text(size=18))
g1 <- g1 + xlab("Measured") + ylab("Predicted")

g2 <- ggplot(df, aes(x=measured, y=prediction)) + geom_point() + stat_cor(method = "pearson", label.x = min(df$measured)*1.1, label.y = max(df$prediction)*0.9, size=10)
g2 <- g2 + theme(text = element_text(size=18))
g2 <- g2 + xlab("Measured Fold Change") + ylab("Predicted Fold Change")

g <- plot_grid(g1, g2, labels=c("A","B"))

png("correlation-vs-fold-change.png", width = 2 * 480)
print(g)
d <- dev.off()