suppressPackageStartupMessages(library("pacman"))
suppressPackageStartupMessages(p_load("synapser"))
suppressPackageStartupMessages(p_load("ggrepel"))
suppressPackageStartupMessages(p_load("ggcorrplot"))
suppressPackageStartupMessages(p_load("ggdendro"))
suppressPackageStartupMessages(p_load("plyr"))
suppressPackageStartupMessages(p_load("dplyr"))
suppressPackageStartupMessages(p_load("reshape2"))

suppressPackageStartupMessages(p_load("MCPcounter"))
suppressPackageStartupMessages(p_load("corrplot"))
suppressPackageStartupMessages(p_load("openxlsx"))

suppressPackageStartupMessages(p_load("sva"))
suppressPackageStartupMessages(p_load("gridExtra"))
suppressPackageStartupMessages(p_load("grid"))
suppressPackageStartupMessages(p_load("gridGraphics"))
suppressPackageStartupMessages(p_load("ComplexHeatmap"))
suppressPackageStartupMessages(p_load("gplots"))
suppressPackageStartupMessages(p_load("scales"))

suppressPackageStartupMessages(p_load("preprocessCore")) # for normalize.quantiles
synLogin()

source("../../utils.R")

## Load the TPM validation data
##synId <- "syn21574299"
##obj <- synGet(synId, downloadFile = TRUE)
##cpm.expr <- read.table(obj$path, sep = "\t", header = TRUE)

synId <- "syn21576632"
obj <- synGet(synId, downloadFile = TRUE)
cpm.expr <- read.table(obj$path, sep = ",", header = TRUE)

## Load the counts
synId <- "syn21576630"
obj <- synGet(synId, downloadFile = TRUE)
cnts.expr <- read.table(obj$path, sep = ",", header = TRUE)

tmp <- cnts.expr[,-1]
round.cnts.expr <- round(tmp)
rownames(round.cnts.expr) <- cnts.expr[,1]

suppressPackageStartupMessages(p_load(edgeR))
log.cpm.expr <- cpm(round.cnts.expr, log=TRUE)

rownames(cpm.expr) <- as.character(cpm.expr$Gene)
cpm.expr <- cpm.expr[, !(colnames(cpm.expr) %in% c("Gene"))]

purified.samples <- c(
    "Naive_B_cells_1",
    "Macrophages_2",
    "Dendritic_cells_1",
    "CRC",
    "Macrophages_1",
    "Fibroblasts",
    "Endothelial_cells",
    "Memory_CD8_T_cells_1",
    "Monocytes_1",
    "Neutrophils_2",
    "Breast",
    "Dendritic_cells_2",
    "NK_cells_2",
    "Monocytes_2",
    "NK_cells_1",
    "Memory_CD8_T_cells_2",
    "Memory_CD4_T_cells_2",
    "Naive_CD4_T_cells_1",
    "Tregs",
    "Memory_CD4_T_cells_1",
    "Naive_CD8_T_cells_2",
    "Naive_CD4_T_cells_2")

rename.samples <- function(mat) {
    lst <- list(
        "Naive_B_cells" = "B_naive",
        "Macrophages" = "Macro",
        "Dendritic_cells" = "DC",
        "Endothelial_cells" = "Endo",
        "Memory_CD8_T_cells" = "CD8T_mem",
        "Naive_CD8_T_cells" = "CD8T_naive",
        "Memory_CD4_T_cells" = "CD4T_mem",
        "Naive_CD4_T_cells" = "CD4T_naive",
        "Monocytes" = "Mono",
        "Neutrophils" = "Neutro",
        "NK_cells" = "NK")
    for(entry in names(lst)) {
      colnames(mat) <- gsub(colnames(mat), pattern = entry, replacement = lst[[entry]])
    }
    mat
}

purified.cpm.expr <- cpm.expr[, purified.samples]
purified.cpm.expr <- rename.samples(purified.cpm.expr)

purified.log.cpm.expr <- log.cpm.expr[, purified.samples]
purified.log.cpm.expr <- rename.samples(purified.log.cpm.expr)

all.zero <- unlist(apply(purified.cpm.expr, 1, function(row) all(row == 0)))
zero.genes <- rownames(purified.cpm.expr)[all.zero]

purified.cpm.expr <- purified.cpm.expr[!(rownames(purified.cpm.expr) %in% zero.genes),]
purified.log.cpm.expr <- purified.log.cpm.expr[!(rownames(purified.log.cpm.expr) %in% zero.genes),]

subset.expressed.and.highly.variable.trend.line <- function(log.expr.mat) {

    ## Filter genes whose max expression (across samples) is less than mean expression
    ## minus 2 standard deviations
    max.log.expr <- unlist(apply(log.expr.mat, 1, function(row) max(row, na.rm=TRUE)))
    
    ret <- calc.density.of.expressed.genes(log.expr.mat)
    cutoff <- ret$z$mean - 2 * ret$z$sd

    g1 <- plot.density.of.expressed.genes(log.expr.mat)
    g1 <- g1 + geom_vline(xintercept = ret$z$mean)
    g1 <- g1 + geom_vline(xintercept = cutoff, linetype = "dashed")
    
    flag <- max.log.expr > cutoff
    expressed.log.expr.mat <- log.expr.mat[flag,]
    
    ## Filter out expressed genes whose variance is below the mean-variance trend line
    cat("Calculating mean / variance trend\n")
    mean.var.df <- calc.mean.variance(expressed.log.expr.mat)
    
    g2 <- ggplot_smooth_scatter(mean.var.df, aes(x = mean, y = sqrt.std))
    g2 <- g2 + geom_line(data = mean.var.df, aes(x = mean, y = sqrt.std.pred))
    g2 <- g2 + theme(legend.position = "none") 
    g2 <- g2 + xlab(expression("Mean"~log[2]~"(CPM)")) + ylab(expression(atop("Sqrt Std Dev",~log[2]~"(CPM)")))
    

    highly.variably.expressed.genes <- rownames(subset(mean.var.df, sqrt.std > sqrt.std.pred))
    highly.variably.expressed.log.expr.mat <- expressed.log.expr.mat[highly.variably.expressed.genes, ]

    return(list("filtered.mat" = highly.variably.expressed.log.expr.mat, "g1" = g1, "g2" = g2,
                "highly.variably.expressed.genes" = highly.variably.expressed.genes))
}

res <- subset.expressed.and.highly.variable.trend.line(purified.log.cpm.expr)


highly.variably.expressed.genes <- res$highly.variably.expressed.genes
highly.variably.expressed.purified.log.cpm.expr <- purified.log.cpm.expr[highly.variably.expressed.genes, ]
highly.variably.expressed.purified.cpm.expr <- purified.cpm.expr[highly.variably.expressed.genes, ]

if(FALSE) {
plot.pca <- function(mat) {
    pc <- prcomp(mat, scale = TRUE, center = TRUE)
    df <- data.frame(PC1 = pc$rotation[,1], PC2 =  pc$rotation[,2], label = rownames(pc$rotation))    
    g <- ggplot(data = df, aes(x = PC1, y = PC2, label = label))
    g <- g + geom_point()
    g <- g + geom_text_repel()
    g
}
}


plot.pca <- function(mat, colour = NULL, shape = NULL, pcs = c(1,2), show.labels = TRUE) {
    pc <- mat
    if(class(mat) != "prcomp") {
        constant <- unlist(apply(mat, 1, function(row) all(row == row[1])))
        mat <- mat[!constant,]
        pc <- prcomp(t(mat), scale = TRUE, center = TRUE)
    }
    id <- rownames(pc$x)
    df <- data.frame(id = id, x = pc$x[,pcs[1]], y = pc$x[,pcs[2]])
    if(!is.null(colour)) {
        colour.df <- data.frame(id = names(colour), colour = unname(colour))
        df <- merge(df, colour.df)
    }
    if(!is.null(shape)) {
        shape.df <- data.frame(id = names(shape), shape = unname(shape))
        df <- merge(df, shape.df)
    }
    df$id <- gsub(df$id, pattern = "_1", replacement = "")
    df$id <- gsub(df$id, pattern = "_2", replacement = "")    
    g <- ggplot(data = df, aes(x = x, y = y))
    if(!is.null(colour) && !is.null(shape)) {
        g <- g + geom_point(aes(colour = colour, shape = shape))
    } else if(!is.null(colour)) {
        g <- g + geom_point(aes(colour = colour))
    } else if(!is.null(shape)) {
        g <- g + geom_point(aes(shape = shape))
    } else {
        g <- g + geom_point()
    }
    if(show.labels) {
        g <- g + geom_text_repel(aes(label = id))
    }
    pc1.var <- round(100 * (summary(pc)$importance)["Proportion of Variance",paste0("PC", pcs[1])])
    pc2.var <- round(100 * (summary(pc)$importance)["Proportion of Variance",paste0("PC", pcs[2])])
    g <- g + xlab(paste0("PC", pcs[1], " (", pc1.var, "%)"))
    g <- g + ylab(paste0("PC", pcs[2], " (", pc2.var, "%)"))
    g
}


plot.correlation <- function(mat) {
    corr <- round(cor(mat), 1)
    ## g <- ggcorrplot(corr, type = "lower", method = "circle", hc.order = TRUE)
    g <- ggcorrplot(corr, method = "circle", hc.order = TRUE)
    ## The circles are too big. Scale them down so they don't overlap
    g <- g + scale_size(range = c(4, 7.5))
    ## g <- ggcorrplot(corr, type = "lower", hc.order = TRUE)
    g
}

## pairing[s] gives the sample pair of s. e.g., s might be an in vitro sample
## and pairing[s] the corresponding in silico sample
plot.neighbor.correlation <- function(mat, pairing) {
    cr <- cor(mat)
    ## Subtract off the diagonal
    crp <- cr - diag(nrow = nrow(cr), ncol = ncol(cr))

    nms <- names(pairing)
    names(nms) <- nms
    ## Plot correlation as a function of distance/rank
    top.ns <- 1:10
    tbl <-
        ldply(nms,
              .fun = function(col) {
                  vec <- crp[col,]
                  vec <- vec[order(vec, decreasing = TRUE)]
                  paired.cor <- crp[col, pairing[col]]
                  data.frame(rank = c("paired", top.ns), c(paired.cor, as.numeric(vec)[top.ns]))
              })

    colnames(tbl) <- c("id", "rank", "val")
    tbl$rank <- factor(tbl$rank, c("paired", top.ns))
    g <- ggplot(data = tbl, aes(x = rank, y = val))
    g <- g + geom_boxplot()
    g <- g + geom_point()
    g <- g + xlab("Rank") + ylab("Pearson Correlation")
    g <- g + theme(text = element_text(size = 20))
    g <- g + geom_hline(yintercept = median(subset(tbl, rank == "paired")$val), linetype = "dashed")
    
    g
}

plot.dendrogram <- function(mat) {
    dst <- dist(t(as.matrix(mat)))
    hc <- hclust(dst^2, "cen")
    hcdata <- dendro_data(hc)
    
    g <- ggdendrogram(hcdata, rotate = TRUE)
    g <- g + theme(text = element_text(size= 20))
    return(g)
    
    g <- ggplot(segment(hcdata)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
    g <- g + geom_text(aes(x = x, y = y, label = label, angle = 45, hjust = 1), data= label(hcdata), size = 5)
    g <- g + scale_x_continuous(expand = c(0.1, 0.1)) + scale_y_continuous(expand = c(0.2,0.1))
    g
}

    
make.all.plots <- function(mat, log.mat, title, out.png) {
    g.pca <- plot.pca(log.mat)
    ## g.pca <- g.pca + ggtitle(title)


    g.corr <- plot.correlation(mat)

    g.dendro <- plot.dendrogram(log.mat)
    
    res <- list("g.pca" = g.pca, "g.corr" = g.corr, "g.dendro" = g.dendro)

    ## pdf(out.pdf, width = 14, height = 14)
    png(out.png, width = 2 * 480, height = 2 * 480)
    title.grob <- textGrob(title, gp = gpar(fontsize = 20))
    grid.arrange(res$g.pca, res$g.corr, res$g.dendro, layout_matrix = rbind(c(1,2),c(3,3)), top = title.grob)
    d <- dev.off()

    res
}

deconv.genes <- get.deconvolution.genes()

deconv.purified.log.cpm.expr <- purified.log.cpm.expr[rownames(purified.log.cpm.expr) %in% deconv.genes,]
deconv.purified.cpm.expr <- purified.cpm.expr[rownames(purified.cpm.expr) %in% deconv.genes,]

title <- "Highly Variable/Expressed Genes"
out.png <- "highly-variable-expressed-plots.png"
make.all.plots(highly.variably.expressed.purified.cpm.expr,
               highly.variably.expressed.purified.log.cpm.expr, title = title, out.png)


title <- "All Genes"
out.png <- "all-genes-plots.png"
make.all.plots(purified.cpm.expr, purified.log.cpm.expr, title = title, out.png) 

title <- "Deconvolution Genes"
out.png <- "deconvolution-genes-plots.png"
make.all.plots(deconv.purified.cpm.expr, deconv.purified.log.cpm.expr, title = title, out.png) 

make.pca.plots.with.factors <- function(mat,
                                        insilico.vs.invitro.factor,
                                        bm.vs.rm.factor,
                                        tumor.factor,
                                        title, out.png) {
    
    g1 <- plot.pca(mat, colour = insilico.vs.invitro.factor, shape = bm.vs.rm.factor, show.labels = FALSE)
    g1 <- g1 + ggtitle("All Admxitures")
    
    g2 <- plot.pca(mat[,bm.samples], colour = insilico.vs.invitro.factor[bm.samples],
                   shape = tumor.factor[bm.samples], show.labels = FALSE)
    g2 <- g2 + ggtitle("Biological Admxitures")
    
    g3 <- plot.pca(mat[,rm.samples], colour = insilico.vs.invitro.factor[rm.samples],
                   shape = tumor.factor[rm.samples], show.labels = FALSE)
    g3 <- g3 + ggtitle("Random Admxitures")
    
    res <- list("g1" = g1, "g2" = g2, "g3" = g3)

    ## pdf(out.pdf, width = 14, height = 14)
    png(out.png, width = 2 * 480, height = 2 * 480)
    title.grob <- textGrob(title, gp = gpar(fontsize = 20))
    grid.arrange(res$g1, res$g2, res$g3, top = title.grob)
    d <- dev.off()

    res
}

## pairing[s] gives the sample pair of s. e.g., s might be an in vitro sample
## and pairing[s] the corresponding in silico sample
make.paired.in.vitro.vs.in.silico.plots <- function(mat, pairing, title, out.png) {
    
    g.dendro <- plot.dendrogram(mat)

    ## g.corr <- plot.correlation(mat)

    corr <- cor(mat)
    corr <- corr[pairing, names(pairing)]
    g.corr <- ggcorrplot(corr, method = "circle", hc.order = FALSE)
    g.corr <- g.corr + scale_size(range = c(4, 7.5))
    
    g.neighbor <- plot.neighbor.correlation(mat, pairing = pairing)
    
    res <- list("g.corr" = g.corr, "g.dendro" = g.dendro, "g.neighbor" = g.neighbor)

    ## pdf(out.pdf, width = 14, height = 14)
    png(out.png, width = 2 * 480, height = 2 * 480)
    title.grob <- textGrob(title, gp = gpar(fontsize = 20))
    grid.arrange(res$g.dendro, res$g.corr, res$g.neighbor, layout_matrix = rbind(c(1,1),c(2,3)), top = title.grob)
    d <- dev.off()

    res
}

make.all.paired.plots <- function(mat, pairing, bm.samples, rm.samples, crc.samples, brca.samples,
                                  title.prefix, png.prefix) {

    samples <- colnames(mat)[colnames(mat) %in% intersect(bm.samples, crc.samples)]
    title <- paste0(title.prefix, " (Biological; CRC)")
    out.png <- paste0(png.prefix, "-biological-crc-paired.png")
    tmp <- mat[, samples]
    make.paired.in.vitro.vs.in.silico.plots(tmp, pairing = pairing[pairing %in% samples], title, out.png)

    samples <- colnames(mat)[colnames(mat) %in% intersect(bm.samples, brca.samples)]
    title <- paste0(title.prefix, " (Biological; BRCA)")
    out.png <- paste0(png.prefix, "-biological-brca-paired.png")
    tmp <- mat[, samples]
    make.paired.in.vitro.vs.in.silico.plots(tmp, pairing = pairing[pairing %in% samples], title, out.png)

    samples <- colnames(mat)[colnames(mat) %in% intersect(rm.samples, crc.samples)]
    title <- paste0(title.prefix, " (Random; CRC)")
    out.png <- paste0(png.prefix, "-random-crc-paired.png")
    tmp <- mat[, samples]
    make.paired.in.vitro.vs.in.silico.plots(tmp, pairing = pairing[pairing %in% samples], title, out.png)

    samples <- colnames(mat)[colnames(mat) %in% intersect(rm.samples, brca.samples)]
    title <- paste0(title.prefix, " (Random; BRCA)")
    out.png <- paste0(png.prefix, "-random-brca-paired.png")
    tmp <- mat[, samples]
    make.paired.in.vitro.vs.in.silico.plots(tmp, pairing = pairing[pairing %in% samples], title, out.png)
}

## Create in silico admixtures mimicking the in vitro admixtures

## Get the ground truth that we will use to create the in silico admixtures.
## NB: these exclude admixtures for which we don't have all corresponding purified profiles ...
sample.ratios <- get.actual.validation.ground.truth.sample.names()

invitro.cpm.admixtures <- cpm.expr[, !(colnames(cpm.expr) %in% purified.samples)]

## ... hence the need to take this intersection of samples
admix.names <- intersect(colnames(invitro.cpm.admixtures), unique(sample.ratios$sample))
invitro.cpm.admixtures <- invitro.cpm.admixtures[, admix.names]
sample.ratios <- subset(sample.ratios, sample %in% admix.names)

insilico.cpm.admixtures <-
    create.in.silico.admixtures(cpm.expr, sample.ratios,
                                sample.col = "sample", cell.type.col = "cell.type",
                                fraction.col = "actual")

## Combine in silico and in vitro admixtures after appending "s" to the former and "v" to the latter
colnames(insilico.cpm.admixtures) <- paste0(colnames(insilico.cpm.admixtures), "s")
colnames(invitro.cpm.admixtures) <- paste0(colnames(invitro.cpm.admixtures), "v")
combined.cpm.admixtures <- cbind(insilico.cpm.admixtures, invitro.cpm.admixtures)

insilico.samples <- colnames(insilico.cpm.admixtures)
invitro.samples <- colnames(invitro.cpm.admixtures)
insilico.vs.invitro.factor <- c(rep("in silico", length(insilico.samples)),
                                rep("in vitro", length(invitro.samples)))
names(insilico.vs.invitro.factor) <- c(insilico.samples, invitro.samples)

crc.samples <- subset(sample.ratios, cell.type == "CRC" & actual > 0)$sample
crc.samples <- c(paste0(crc.samples, "s"), paste0(crc.samples, "v"))
brca.samples <- subset(sample.ratios, cell.type == "Breast" & actual > 0)$sample
brca.samples <- c(paste0(brca.samples, "s"), paste0(brca.samples, "v"))
tumor.factor <- c(rep("CRC", length(crc.samples)), rep("BRCA", length(brca.samples)))
names(tumor.factor) <- c(crc.samples, brca.samples)

## Colour by random vs biological admixture
bm.samples <- colnames(combined.cpm.admixtures)[grepl(colnames(combined.cpm.admixtures), pattern="BM")]
rm.samples <- colnames(combined.cpm.admixtures)[grepl(colnames(combined.cpm.admixtures), pattern="RM")]
bm.vs.rm.factor <- c(rep("Biological", length(bm.samples)), rep("Random", length(rm.samples)))
names(bm.vs.rm.factor) <- c(bm.samples, rm.samples)

title <- "All Genes"
out.png <- "all-genes-admixture-pca-plots.png"
tmp <- combined.cpm.admixtures
make.pca.plots.with.factors(tmp, insilico.vs.invitro.factor, bm.vs.rm.factor, tumor.factor, title, out.png)

title <- "Highly Variable/Expressed Genes"
out.png <- "highly-variable-expressed-admixture-pca-plots.png"
flag <- rownames(combined.cpm.admixtures) %in% highly.variably.expressed.genes
tmp <- combined.cpm.admixtures[flag, ]
make.pca.plots.with.factors(tmp, insilico.vs.invitro.factor, bm.vs.rm.factor, tumor.factor, title, out.png)

title <- "Deconvolution Genes"
out.png <- "deconvolution-genes-admixture-pca-plots.png"
flag <- rownames(combined.cpm.admixtures) %in% deconv.genes
tmp <- combined.cpm.admixtures[flag, ]
make.pca.plots.with.factors(tmp, insilico.vs.invitro.factor, bm.vs.rm.factor, tumor.factor, title, out.png)

pairing <- colnames(insilico.cpm.admixtures)
names(pairing) <- pairing
pairing <- gsub(pairing, pattern="s$", replacement="v")

title.prefix <- "All Genes"
png.prefix <- "all-genes"
tmp <- combined.cpm.admixtures
make.all.paired.plots(tmp, pairing, bm.samples, rm.samples, crc.samples, brca.samples, title.prefix, png.prefix)

title.prefix <- "Highly Variable/Expressed Genes"
png.prefix <- "highly-variable-expressed-genes"
flag <- rownames(combined.cpm.admixtures) %in% highly.variably.expressed.genes
tmp <- combined.cpm.admixtures[flag, ]
make.all.paired.plots(tmp, pairing, bm.samples, rm.samples, crc.samples, brca.samples, title.prefix, png.prefix)


title.prefix <- "Deconvolution Genes"
png.prefix <- "deconvolution-genes"
flag <- rownames(combined.cpm.admixtures) %in% deconv.genes
tmp <- combined.cpm.admixtures[flag, ]
make.all.paired.plots(tmp, pairing, bm.samples, rm.samples, crc.samples, brca.samples, title.prefix, png.prefix)


## Repeat above for count-based admixtures
## invitro.cnt.admixtures <- cnts.expr[, !(colnames(cnts.expr) %in% purified.samples)]
invitro.cnt.admixtures <- cnts.expr[, admix.names]
rownames(invitro.cnt.admixtures) <- cnts.expr$Gene

purified.cnts.expr <- cnts.expr[, purified.samples]
rownames(purified.cnts.expr) <- cnts.expr$Gene
tot.purified.cnts <- colSums(purified.cnts.expr)
median.tot.purified.cnts <- median(tot.purified.cnts)

norm.purified.cnts.expr <- sweep(purified.cnts.expr, 2, colSums(purified.cnts.expr),`/`) * median.tot.purified.cnts

insilico.cnt.admixtures <-
    create.in.silico.admixtures(norm.purified.cnts.expr, sample.ratios,
                                sample.col = "sample", cell.type.col = "cell.type",
                                fraction.col = "actual")

## Combine in silico and in vitro admixtures after appending "s" to the former and "v" to the latter
colnames(insilico.cnt.admixtures) <- paste0(colnames(insilico.cnt.admixtures), "s")
colnames(invitro.cnt.admixtures) <- paste0(colnames(invitro.cnt.admixtures), "v")
combined.cnt.admixtures <- cbind(insilico.cnt.admixtures, invitro.cnt.admixtures)
combined.cpm.cnt.admixtures <- cpm(combined.cnt.admixtures, log = FALSE)

do.all <- function(insilico.mat, invitro.mat, combined.mat, sample.ratios, png.suffix, title.suffix) {

    insilico.samples <- colnames(insilico.mat)
    invitro.samples <- colnames(invitro.mat)
    insilico.vs.invitro.factor <- c(rep("in silico", length(insilico.samples)),
                                    rep("in vitro", length(invitro.samples)))
    names(insilico.vs.invitro.factor) <- c(insilico.samples, invitro.samples)

    crc.samples <- subset(sample.ratios, cell.type == "CRC" & actual > 0)$sample
    crc.samples <- c(paste0(crc.samples, "s"), paste0(crc.samples, "v"))
    brca.samples <- subset(sample.ratios, cell.type == "Breast" & actual > 0)$sample
    brca.samples <- c(paste0(brca.samples, "s"), paste0(brca.samples, "v"))
    tumor.factor <- c(rep("CRC", length(crc.samples)), rep("BRCA", length(brca.samples)))
    names(tumor.factor) <- c(crc.samples, brca.samples)

    ## Colour by random vs biological admixture
    bm.samples <- colnames(combined.mat)[grepl(colnames(combined.mat), pattern="BM")]
    rm.samples <- colnames(combined.mat)[grepl(colnames(combined.mat), pattern="RM")]
    bm.vs.rm.factor <- c(rep("Biological", length(bm.samples)), rep("Random", length(rm.samples)))
    names(bm.vs.rm.factor) <- c(bm.samples, rm.samples)

    title <- paste0("All Genes", title.suffix)
    out.png <- paste0("all-genes-admixture-pca-plots", png.suffix, ".png")
    tmp <- combined.mat
    make.pca.plots.with.factors(tmp, insilico.vs.invitro.factor, bm.vs.rm.factor, tumor.factor, title, out.png)

    title <- paste0("Highly Variable/Expressed Genes", title.suffix)
    out.png <- paste0("highly-variable-expressed-admixture-pca-plots", png.suffix, ".png")
    flag <- rownames(combined.mat) %in% highly.variably.expressed.genes
    tmp <- combined.mat[flag, ]
    make.pca.plots.with.factors(tmp, insilico.vs.invitro.factor, bm.vs.rm.factor, tumor.factor, title, out.png)

    title <- paste0("Deconvolution Genes", title.suffix)
    out.png <- paste0("deconvolution-genes-admixture-pca-plots", png.suffix, ".png")
    flag <- rownames(combined.mat) %in% deconv.genes
    tmp <- combined.mat[flag, ]
    make.pca.plots.with.factors(tmp, insilico.vs.invitro.factor, bm.vs.rm.factor, tumor.factor, title, out.png)

    pairing <- colnames(insilico.mat)
    names(pairing) <- pairing
    pairing <- gsub(pairing, pattern="s$", replacement="v")

    title.prefix <- paste0("All Genes", title.suffix)
    png.prefix <- paste0("all-genes", png.suffix)
    tmp <- combined.mat
    make.all.paired.plots(tmp, pairing, bm.samples, rm.samples, crc.samples, brca.samples, title.prefix, png.prefix)
    
    title.prefix <- paste0("Highly Variable/Expressed Genes", title.suffix)
    png.prefix <- paste0("highly-variable-expressed-genes", png.suffix)
    flag <- rownames(combined.mat) %in% highly.variably.expressed.genes
    tmp <- combined.mat[flag, ]
    make.all.paired.plots(tmp, pairing, bm.samples, rm.samples, crc.samples, brca.samples, title.prefix, png.prefix)
    
    title.prefix <- paste0("Deconvolution Genes", title.suffix)
    png.prefix <- paste0("deconvolution-genes", png.suffix)
    flag <- rownames(combined.mat) %in% deconv.genes
    tmp <- combined.mat[flag, ]
    make.all.paired.plots(tmp, pairing, bm.samples, rm.samples, crc.samples, brca.samples, title.prefix, png.prefix)
}

png.suffix <- "-cpm-cnt-mixed"
title.suffix <- "; CPM of mixed counts"
do.all(insilico.cnt.admixtures, invitro.cnt.admixtures, combined.cpm.cnt.admixtures, sample.ratios, png.suffix, title.suffix)

png.suffix <- "-cpm-mixed"
title.suffix <- "; mixtures of CPMs"
do.all(insilico.cpm.admixtures, invitro.cpm.admixtures, combined.cpm.admixtures, sample.ratios, png.suffix, title.suffix)


cat("Exiting\n")
q(status = 0)



suppressPackageStartupMessages(p_load(AnnotationHub))
suppressPackageStartupMessages(p_load(ensembldb))
ah <- AnnotationHub()
flag <- (ah$species == "Homo sapiens") & (ah$genome == "GRCh38") & (ah$dataprovider == "Ensembl") & (ah$rdataclass == "EnsDb")
ah2 <- ah[flag, ]
## as.data.frame(mcols(ah2))[1:10,c("title"),drop=FALSE]
nm <- "AH73881"
if(!(nm %in% names(ah2))) { nm <- names(ah2)[1] }
edb <- ah2[[nm]]

## keytypes(edb)
## columns(edb)
keys <- keys(edb, "GENENAME")
columns <- c("GENEID", "ENTREZID", "GENEBIOTYPE")
tbl <- ensembldb::select(edb, keys, columns, keytype = "GENENAME")
pc.tbl <- subset(tbl, GENEBIOTYPE == "protein_coding")

cpm.pc.expr <- cpm.expr[rownames(cpm.expr) %in% pc.tbl$GENENAME, ]
## Renormalize TPMs to protein coding genes
cpm.pc.expr <- apply(cpm.pc.expr, 2, function(col) col * 10^6 / sum(col))
purified.pc.expr <- cpm.pc.expr[, purified.samples]
purified.pc.expr <- rename.samples(purified.pc.expr)

vendors <- get.vendor.sample.assignments()



admix.names <- intersect(colnames(insilico.admixtures), colnames(admixture.expr))
insilico.admixtures <- insilico.admixtures[, admix.names]
admixture.expr <- admixture.expr[, admix.names]

insilico.mat <- cbind(Gene = rownames(purified.mat), insilico.admixtures)

rand.cols <- c("Gene", colnames(insilico.mat)[grepl(colnames(insilico.mat), pattern="^R")])
bio.cols <- c("Gene", colnames(insilico.mat)[grepl(colnames(insilico.mat), pattern="^B")])
rand.mat <- insilico.mat[, rand.cols]
bio.mat <- insilico.mat[, bio.cols]


batch <- c(rep(0, ncol(insilico.admixtures)), rep(1, ncol(admixture.expr)))
flag <- apply(mat, 1, function(row) all(row == row[1]))
mat.c <- ComBat(as.matrix(mat)[!flag,], batch = batch)

## NB: correlations are all very high.
## Limit to deconv genes



all.zero <- unlist(apply(purified.expr, 1, function(row) all(row == 0)))
purified.expr <- purified.expr[!all.zero,]

means <- unlist(apply(purified.expr, 1, mean))
vars <- unlist(apply(purified.expr, 1, var))
cov <- vars / means

cov.genes <- names(tail(cov[order(cov)],n=100))




dst <- dist(t(as.matrix(mat.c[rownames(mat.c) %in% deconv.genes,])))
hc <- hclust(dst^2, "cen")
str(as.dendrogram(hc))

mt <- mat.c

## Plot in vitro vs in silico PCA (deconvolution genes)
pc <- prcomp(mt[rownames(mt) %in% deconv.genes,], scale = TRUE, center = TRUE)
df <- data.frame(PC1 = pc$rotation[,1], PC2 =  pc$rotation[,2])
batch <- rep("in vitro", nrow(df))
flag <- grepl(rownames(df), pattern="i")
batch[flag] <- "in silico"
pdf("gp.pdf")
g <- ggplot(data = df, aes(x = PC1, y = PC2, colour = batch))
g <- g + geom_point()
print(g)
d <- dev.off()



## Plot in vitro vs in silico correlation plot (deconvolution genes)
head(cor(mat.c[rownames(mat.c) %in% deconv.genes,])[,1:5])
pdf("cr.pdf", width = 14)
corrplot(cr, order="hclust")
dst <- dist(t(as.matrix(mat[rownames(mat) %in% deconv.genes,])))
hc <- hclust(dst^2, "cen")
## plot(hc)
d <- dev.off()

dst <- dist(t(as.matrix(mat[rownames(mat) %in% deconv.genes,])))
hc <- hclust(dst^2, "cen")
str(as.dendrogram(hc))

cov.purified.expr <- purified.expr[rownames(purified.expr) %in% cov.genes,]
corrplot(cor(cov.purified.expr), order="hclust")


deconv.purified.expr <- purified.expr[rownames(purified.expr) %in% deconv.genes,]
corrplot(cor(deconv.purified.expr), order="hclust")

foo <- normalize.quantiles(as.matrix(deconv.purified.expr))
colnames(foo) <- colnames(deconv.purified.expr)
m <- as.matrix(foo)
heatmap(as.matrix(t(scale(t(m)))), scale="none")

m <- as.matrix(deconv.purified.expr)
heatmap(as.matrix(t(scale(t(m)))), scale="none")

pc <- prcomp(t(purified.expr), scale = TRUE, center = TRUE)

pc <- prcomp(t(foo), scale = TRUE, center = TRUE)
plot(pc$x[,1], pc$x[,2])
text(pc$x[,1], pc$x[,2], rownames(pc$x))
plot(pc$x[,1], pc$x[,2])
text(pc$x[,1], pc$x[,2], vendors[rownames(pc$x), "vendor"])
g <- ggplot(data = data.frame(x = pc$x[,1], y = pc$x[,2], colour = vendors[rownames(pc$x), "vendor"],
                              label = rownames(pc$x)),
            aes(x = x, y = y))
g <- g + geom_point(aes(colour = colour))
g <- g + geom_text_repel(aes(label = label))
pc1.var <- round(100 * (summary(pc)$importance)["Proportion of Variance","PC1"])
pc2.var <- round(100 * (summary(pc)$importance)["Proportion of Variance","PC2"])
g <- g + xlab(paste0("PC1 (", pc1.var, "%)"))
g <- g + ylab(paste0("PC2 (", pc2.var, "%)"))

colour <- vendors$vendor
names(colour) <- vendors$sample

purified.expr.qn <- normalize.quantiles(as.matrix(purified.expr))
rownames(purified.expr.qn) <- rownames(purified.expr)
colnames(purified.expr.qn) <- colnames(purified.expr)
all.zero <- unlist(apply(purified.expr.qn, 1, function(row) all(row == row[1])))
purified.expr.qn <- purified.expr.qn[!all.zero,]

deconv.purified.expr.qn <- purified.expr.qn[rownames(purified.expr.qn) %in% deconv.genes,]

cols <- c("uid", "accession", "gds", "title", "summary", "gpl", "gse", "taxon", "entrytype", "gdstype",
          "ptechtype", "valtype", "ssinfo", "subsetinfo", "pdat", "suppfile", "samples", "relations",
          "extrelations", "n_samples", "seriestitle", "platformtitle", "platformtaxa", "samplestaxa",
          "pubmedids", "projects", "ftplink", "geo2r")
