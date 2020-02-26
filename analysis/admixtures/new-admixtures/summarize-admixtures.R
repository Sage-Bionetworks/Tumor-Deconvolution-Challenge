
tumor.col <- "tumor.fraction"
## drop the tumor column before creating correlations?
drop.tumor.col <- TRUE

##mask <- matrix(data=FALSE, nrow=nrow(broken.stick.cor), ncol=ncol(broken.stick.cor))
##mask[1:10,1:10] <- TRUE
##mask[11:20,11:20] <- TRUE
##mask[21:30,21:20] <- TRUE

##same.cors <- broken.stick.cor[mask & lower.tri(broken.stick.cor)]
##diff.cors <- broken.stick.cor[!mask & lower.tri(broken.stick.cor)]

broken.stick.cols <- colnames(broken.stick.admixtures)
broken.stick.unconstrained.cols <- colnames(broken.stick.admixtures.unconstrained)
admixture.cols <- colnames(admixture.tbl)

if(drop.tumor.col) {
  broken.stick.cols <- broken.stick.cols[!(broken.stick.cols == tumor.col)]
  broken.stick.unconstrained.cols <- broken.stick.unconstrained.cols[!(broken.stick.unconstrained.cols == tumor.col)]
  admixture.cols <- admixture.cols[!(admixture.cols == tumor.col)]
}

broken.stick.cor <- cor(t(broken.stick.admixtures[, broken.stick.cols]), method="spearman")
broken.stick.unconstrained.cor <- cor(t(broken.stick.admixtures.unconstrained[, broken.stick.unconstrained.cols]), method="spearman")
bio.cor <- cor(t(admixture.tbl[, admixture.cols]), method="spearman")

broken.cors <- broken.stick.cor[lower.tri(broken.stick.cor)]
broken.unconstrained.cors <- broken.stick.unconstrained.cor[lower.tri(broken.stick.unconstrained.cor)]
bio.cors <- bio.cor[lower.tri(bio.cor)]

h1 <- hist(bio.cors, plot=FALSE)
h2 <- hist(broken.cors, plot=FALSE)
h3 <- hist(broken.unconstrained.cors, plot=FALSE)

g1 <- ggplot(data = data.frame(x = bio.cors),aes(x=x)) + geom_histogram() + xlab("Spearman Correlation")
g1 <- g1 + ggtitle("Admixture-wise biological correlations")

g2 <- ggplot(data = data.frame(x = broken.cors),aes(x=x)) + geom_histogram() + xlab("Spearman Correlation")
g2 <- g2 + ggtitle("Admixture-wise random correlations")

g3 <- ggplot(data = data.frame(x = broken.unconstrained.cors),aes(x=x)) + geom_histogram() + xlab("Spearman Correlation")
g3 <- g3 + ggtitle("Admixture-wise random (unconstrained) correlations")

xmin <- min(ggplot_build(g1)$layout$panel_params[[1]]$x.range, ggplot_build(g2)$layout$panel_params[[1]]$x.range, ggplot_build(g3)$layout$panel_params[[1]]$x.range)
xmax <- max(ggplot_build(g1)$layout$panel_params[[1]]$x.range, ggplot_build(g2)$layout$panel_params[[1]]$x.range, ggplot_build(g3)$layout$panel_params[[1]]$x.range) 
ymin <- min(ggplot_build(g1)$layout$panel_params[[1]]$y.range, ggplot_build(g2)$layout$panel_params[[1]]$y.range, ggplot_build(g3)$layout$panel_params[[1]]$y.range)
ymax <- max(ggplot_build(g1)$layout$panel_params[[1]]$y.range, ggplot_build(g2)$layout$panel_params[[1]]$y.range, ggplot_build(g3)$layout$panel_params[[1]]$y.range)

g1 <- g1 + xlim(c(xmin, xmax)) + ylim(c(ymin, ymax))
g2 <- g2 + xlim(c(xmin, xmax)) + ylim(c(ymin, ymax))
g3 <- g3 + xlim(c(xmin, xmax)) + ylim(c(ymin, ymax))

library(gridExtra)
file <- paste0(file.prefix, "admixture-correlation-histograms.png")
png(file)
grid.arrange(g1, g2, g3)
d <- dev.off()

broken.stick.cell.type.cor <- cor(broken.stick.admixtures[, broken.stick.cols], method="spearman")
broken.stick.cell.type.unconstrained.cor <- cor(broken.stick.admixtures.unconstrained[, broken.stick.unconstrained.cols], method="spearman")
bio.cell.type.cor <- cor(admixture.tbl[, admixture.cols], method="spearman")

broken.cell.type.cors <- broken.stick.cell.type.cor[lower.tri(broken.stick.cell.type.cor)]
broken.cell.type.unconstrained.cors <- broken.stick.cell.type.unconstrained.cor[lower.tri(broken.stick.cell.type.unconstrained.cor)]
bio.cell.type.cors <- bio.cell.type.cor[lower.tri(bio.cell.type.cor)]

g1 <- ggplot(data = data.frame(x = bio.cell.type.cors),aes(x=x)) + geom_histogram() + xlab("Spearman Correlation")
g1 <- g1 + ggtitle("Cell-wise biological correlations")
g2 <- ggplot(data = data.frame(x = broken.cell.type.cors),aes(x=x)) + geom_histogram() + xlab("Spearman Correlation")
g2 <- g2 + ggtitle("Cell-wise random correlations")
g3 <- ggplot(data = data.frame(x = broken.cell.type.unconstrained.cors),aes(x=x)) + geom_histogram() + xlab("Spearman Correlation")
g3 <- g3 + ggtitle("Cell-wise random (unconstrained) correlations")

xmin <- min(ggplot_build(g1)$layout$panel_params[[1]]$x.range, ggplot_build(g2)$layout$panel_params[[1]]$x.range, ggplot_build(g3)$layout$panel_params[[1]]$x.range)
xmax <- max(ggplot_build(g1)$layout$panel_params[[1]]$x.range, ggplot_build(g2)$layout$panel_params[[1]]$x.range, ggplot_build(g3)$layout$panel_params[[1]]$x.range) 
ymin <- min(ggplot_build(g1)$layout$panel_params[[1]]$y.range, ggplot_build(g2)$layout$panel_params[[1]]$y.range, ggplot_build(g3)$layout$panel_params[[1]]$y.range)
ymax <- max(ggplot_build(g1)$layout$panel_params[[1]]$y.range, ggplot_build(g2)$layout$panel_params[[1]]$y.range, ggplot_build(g3)$layout$panel_params[[1]]$y.range)
g1 <- g1 + xlim(c(xmin, xmax)) + ylim(c(ymin, ymax))
g2 <- g2 + xlim(c(xmin, xmax)) + ylim(c(ymin, ymax))
g3 <- g3 + xlim(c(xmin, xmax)) + ylim(c(ymin, ymax))

library(gridExtra)
file <- paste0(file.prefix, "cell-correlation-histograms.png")
png(file)
grid.arrange(g1, g2, g3)
d <- dev.off()

## broken.stick.admixtures <- generate.random.uniform.admixtures(colnames(flat.model), 30, tumor.type = tumor.col, min.prop = 0.01, lbs = lbs, ubs = ubs)

cat("Done zero'ing out populations\n")
save.image(".Rdata")

library(openxlsx)
wb <- createWorkbook("Admixtures")

addWorksheet(wb, "model")
writeData(wb, sheet = 1, flat.model, colNames = TRUE, rowNames = TRUE)

addWorksheet(wb, "bio")
o <- order(admixture.tbl[, tumor.col])
writeData(wb, sheet = 2, admixture.tbl, colNames = TRUE, rowNames = FALSE)

addWorksheet(wb, "rand")
writeData(wb, sheet = 3, broken.stick.admixtures, colNames = TRUE, rowNames = FALSE)

addWorksheet(wb, "rand-unconstrained")
writeData(wb, sheet = 4, broken.stick.admixtures.unconstrained, colNames = TRUE, rowNames = FALSE)

addWorksheet(wb, "rand-zero")
writeData(wb, sheet = 5, broken.stick.admixtures.zero, colNames = TRUE, rowNames = FALSE)

addWorksheet(wb, "rand-unconstrained-zero")
writeData(wb, sheet = 6, broken.stick.admixtures.unconstrained.zero, colNames = TRUE, rowNames = FALSE)

saveWorkbook(wb, paste0(file.prefix, "admixtures.xlsx"), overwrite = TRUE)


## Now output the same admixtures, but in a simlified XLS sheet that
## only has the admixtures we are interested in.

library(openxlsx)
wb <- createWorkbook("Admixtures")

addWorksheet(wb, "bio")
writeData(wb, sheet = 1, admixture.tbl, colNames = TRUE, rowNames = FALSE)

addWorksheet(wb, "rand-unconstrained-zero")
writeData(wb, sheet = 2, broken.stick.admixtures.unconstrained.zero, colNames = TRUE, rowNames = FALSE)

saveWorkbook(wb, paste0(file.prefix, "simplified-admixtures.xlsx"), overwrite = TRUE)


prefix <- "biological"
file <- paste0(file.prefix, prefix, "-cell-type-and-sample-correlations.png")
png(file)
par(mfrow=c(1,2))
mar <- c(0, 0, 0, 0)
plot.admixture.correlations(admixture.tbl[, admixture.cols], mar = mar)
title(main="Biological\ncell-wise correlations", line=-3)
plot.admixture.correlations(t(admixture.tbl[, admixture.cols]), mar = mar)
title(main="Biological\nadmixture-wise correlations", line=-3)
d <- dev.off()

prefix <- "biological"
file <- paste0(file.prefix, prefix, "-cell-type-correlations.png")
png(file)
plot.admixture.correlations(admixture.tbl[, admixture.cols], main="Biological cell-wise correlations")
d <- dev.off()

prefix <- "biological"
file <- paste0(file.prefix, prefix, "-sample-correlations.png")
png(file)
plot.admixture.correlations(t(admixture.tbl[, admixture.cols]), main="Biological admixture-wise correlations")
d <- dev.off()

prefix <- "random"
file <- paste0(file.prefix, prefix, "-cell-type-and-sample-correlations.png")
png(file)
par(mfrow=c(1,2))
mar <- c(0, 0, 0, 0)
plot.admixture.correlations(broken.stick.admixtures[, broken.stick.cols], mar=mar)
title(main="Random\ncell-wise correlations", line=-3)
plot.admixture.correlations(t(broken.stick.admixtures[, broken.stick.cols]), mar=mar)
title(main="Random\nadmixture-wise correlations", line=-3)
d <- dev.off()


prefix <- "random"
file <- paste0(file.prefix, prefix, "-cell-type-correlations.png")
png(file)
plot.admixture.correlations(broken.stick.admixtures[, broken.stick.cols], main="Random cell-wise correlations")
d <- dev.off()

prefix <- "random"
file <- paste0(file.prefix, prefix, "-sample-correlations.png")
png(file)
plot.admixture.correlations(t(broken.stick.admixtures[, broken.stick.cols]), main="Random admixture-wise correlations")
d <- dev.off()

prefix <- "random-unconstrained"
file <- paste0(file.prefix, prefix, "-cell-type-and-sample-correlations.png")
png(file)
par(mfrow=c(1,2))
mar <- c(0, 0, 0, 0)
plot.admixture.correlations(broken.stick.admixtures.unconstrained[, broken.stick.unconstrained.cols], mar=mar)
title(main="Random (unconstrained)\ncell-wise correlations", line=-3)
plot.admixture.correlations(t(broken.stick.admixtures.unconstrained[, broken.stick.unconstrained.cols]), mar=mar)
title(main="Random (unconstrained)\nadmixture-wise correlations", line=-3)
d <- dev.off()

prefix <- "random-unconstrained"
file <- paste0(file.prefix, prefix, "-cell-type-correlations.png")
png(file)
plot.admixture.correlations(broken.stick.admixtures.unconstrained[, broken.stick.unconstrained.cols], main="Random (unconstrained) cell-wise correlations")
d <- dev.off()

prefix <- "random-unconstrained"
file <- paste0(file.prefix, prefix, "-sample-correlations.png")
png(file)
plot.admixture.correlations(t(broken.stick.admixtures.unconstrained[, broken.stick.unconstrained.cols]), main="Random (unconstrained) admixture-wise correlations")
d <- dev.off()


save.image(".Rdata")
cat("Exiting successfully\n")

q(status=0)
