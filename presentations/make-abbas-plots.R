library(synapser)
library(plyr)
library(ggplot2)
library(grid)


synLogin()

source("../scripts/utils.R")

## Get CIBERSORT genes
cibersort.genes.synId <- "syn12184137"
cibersort.genes <- read.table(synGet(cibersort.genes.synId)$path, sep="\t",
                              header=TRUE, as.is=TRUE)

cibersort.genes <- subset(cibersort.genes, Method == "cibersort")

cd8.t.cells = c("T.cells.CD8")
cd4.t.cells = c("T.cells.CD4.naive",
                "T.cells.CD4.memory.resting",
                "T.cells.CD4.memory.activated",
                "T.cells.regulatory.(Tregs)",
                "T.cells.follicular.helper")

cell.types <-
    list("B cells" = c("B.cells.naive", "B.cells.memory"),
         "CD8 T cells" = c(cd8.t.cells, cd4.t.cells),
         "NK cells" = c("NK.cells.activated", "NK.cells.resting"))
for(cell.type in names(cell.types)) {
    flag <- cibersort.genes$cell_type %in% cell.types[[cell.type]]
    cibersort.genes$cell_type[flag] <- cell.type
}

cibersort.genes <- subset(cibersort.genes,
                          cell_type %in% names(cell.types))

## Load Abbas data
abbas.expr.file <- "GSE22886-GPL96_series_matrix.txt-expr.tsv"
abbas.anno.file <- "GSE22886-GPL96_series_matrix.txt-anno.tsv"

orig.expr <- read.table(abbas.expr.file, sep="\t", header=TRUE)
orig.anno <- read.table(abbas.anno.file, sep="\t", header=TRUE, as.is=TRUE)


expr <- orig.expr
anno <- orig.anno

common.samples <- intersect(orig.anno$geo_accession, colnames(orig.expr))
anno <- subset(orig.anno, geo_accession %in% common.samples)
common.samples <- anno$geo_accession

expr <- orig.expr[, c("Hugo", common.samples)]

cell.type.means <- 
    dlply(anno, .variables = c("cell.type.ch1"),
          .fun = function(anno.cell) {
              expr.sub <- expr[, colnames(expr) %in% anno.cell$geo_accession]
              log2(rowMeans(expr.sub))
          })
tmp <- ldply(cell.type.means)
rownames(tmp) <- tmp$cell.type.ch1
tmp <- tmp[, !(colnames(tmp) == "cell.type.ch1")]
cell.type.means <- t(tmp)

old.cols <- c("B cells", "IgM memory B cells", "CD8+ T cells",
              "CD4+ T cells", "CD4+ CD45RO+ CD45RA- T cells",
              "NK cells")
new.cols <- c("B cells", "Memory B cells", "CD8 T cells",
              "CD4 T cells", "Memory CD4 cells", "NK cells")
cell.type.means <- cell.type.means[, old.cols]
colnames(cell.type.means) <- new.cols

mat <- cell.type.means[, c("B cells", "CD8 T cells", "CD4 T cells", "Memory CD4 cells", "NK cells")]
rownames(mat) <- expr$Hugo


png("abbas-correlation.png")
library(corrplot)
## plot.admixture.correlations(cell.type.means)
fc <- as.matrix(cor(cell.type.means, method = "spearman"))
mar <- c(1, 1, 0, 1) + 0.1
corrplot(fc, method = "ellipse", type = "upper", order = "original", tl.cex = 0.8, diag = FALSE, mar = mar)
d <- dev.off()

cell.types <-
    list("B cells" = c("B cells"),
         ##         "T cells" = c("CD4+ T cells", "CD8+ T cells"),
         "CD8 T cells" = c("CD8+ T cells"),
         "NK cells" = c("NK cells"))
for(cell.type in names(cell.types)) {
    flag <- anno$cell.type.ch1 %in% cell.types[[cell.type]]
    anno$cell.type.ch1[flag] <- cell.type
}

tmp <- as.data.frame(table(anno$cell.type.ch1))
tot <- sum(tmp$Freq)
tmp$proportions = tmp$Freq / tot
colnames(tmp) <- c("cell.type", "frequency", "proportions")
proportions <- tmp


anno <- subset(anno, cell.type.ch1 %in% unique(cibersort.genes$cell_type))

common.samples <- intersect(anno$geo_accession, colnames(expr))

anno <- subset(anno, geo_accession %in% common.samples)
anno <- anno[order(anno$cell.type.ch1),]
common.samples <- anno$geo_accession

expr <- subset(expr, Hugo %in% cibersort.genes$Hugo)



rownames(expr) <- expr[, "Hugo"]
expr <- expr[, !(colnames(expr) == "Hugo")]
common.genes <- cibersort.genes$Hugo
common.genes <- common.genes[common.genes %in% rownames(expr)]
expr <- expr[common.genes,]

mat <- mat[common.genes,]

all.expr <- expr

all.expr.log.scale <- t(scale(t(log2(as.matrix(expr))), center=TRUE, scale=TRUE))
expr <- expr[, common.samples]

expr <- as.matrix(expr)

##heatmap(log2(expr), Rowv=NA, Colv=NA, scale="row")

expr.log.scale <- t(scale(t(log2(expr)), center=TRUE, scale=TRUE))

gene.levels <- rownames(expr.log.scale)
sample.levels <- colnames(expr.log.scale)
cell.type.levels <- unique(anno$cell.type.ch1)

library("tidyr")
expr.long <- expr.log.scale %>%
    as.data.frame() %>%
    rownames_to_column('Var1') %>%
    gather(Var2, value, -Var1) %>%
    mutate(Var2 = factor(Var2, levels = sample.levels)) %>%
    mutate(Var1 = factor(Var1, levels = gene.levels))

## g <- ggplot(data=expr.long, aes(x=Var2, y=Var1, fill=value)) + geom_tile()

remove.axes <- function(g, remove.title = TRUE, remove.axis.labels = TRUE) {
    g <- g + theme(axis.line=element_blank(),
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks=element_blank(),
                   legend.position="none",
                   panel.background=element_blank(),
                   panel.border=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   plot.background=element_blank())
    if(remove.title) { g <- g + labs(title = NULL) }
    if(remove.axis.labels) {
        g <- g + theme(axis.title.x=element_blank(),
                       axis.title.y=element_blank())
        g <- g + labs(x=NULL, y=NULL)
    }
##    g <- g + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
    g
}

plot.col <- function(vals, low.col = "darkblue", mid.col = "white",
                     high.col = "darkgreen", ...) {

    df <- data.frame(fill = vals, y = 1:length(vals))

    g <- ggplot(df, aes(x = 1, y = y, fill = fill))
    g <- g + geom_tile(show.legend=FALSE)
    g <- g + scale_fill_gradient2(low = low.col, mid = mid.col, high = high.col)
    g <- remove.axes(g, ...)
    g
}

low.col <- "darkblue"
high.col <- "darkgreen"
mid.col <- "white"


nsamples <- ncol(expr.log.scale)
##gt <- ggplotGrob(g)
prefix <- "training"
##png(paste0(prefix, ".png"),nsamples*10/4,50)
png(paste0(prefix, ".png"))
g <- ggplot(data=expr.long, aes(x=Var2, y=Var1, fill=value)) 
g <- g + geom_tile()
g <- g + scale_fill_gradient2(low = low.col, mid = mid.col, high = high.col)
g <- remove.axes(g, remove.title=FALSE, remove.axis.labels=FALSE)
g <- g + xlab("Samples") + ylab("Genes")
g <- g + scale_x_discrete(position="top")
g <- g + theme(plot.title = element_text(size = 30),
               axis.title.x = element_text(size = 25),
               axis.title.y = element_text(size = 25))
g <- g + ggtitle("Training Data") + theme(plot.title = element_text(hjust = 0.5))
print(g)
##grid.draw(gtable::gtable_filter(gt, "panel"))
d <- dev.off()





noise.cols <- colnames(all.expr.log.scale)[!(colnames(all.expr.log.scale) %in% anno$geo_accession)]
noise.vec <- rowMeans(all.expr.log.scale[, noise.cols])
admixture.vec <- rowMeans(all.expr.log.scale)

set.seed(1234)
sample.indices1 <- sample.int(ncol(all.expr.log.scale),
                              size = ncol(all.expr.log.scale), replace = TRUE)
sample.indices2 <- sample.int(ncol(all.expr.log.scale),
                              size = ncol(all.expr.log.scale), replace = TRUE)
sample.indices3 <- sample.int(ncol(all.expr.log.scale),
                              size = ncol(all.expr.log.scale), replace = TRUE)
admixture1.vec <- rowMeans(all.expr.log.scale[, sample.indices1])
admixture2.vec <- rowMeans(all.expr.log.scale[, sample.indices2])
admixture3.vec <- rowMeans(all.expr.log.scale[, sample.indices3])


admixture1.orig.vec <- rowMeans(all.expr)
admixture2.orig.vec <- rowMeans(all.expr[, sample.indices2])
admixture3.orig.vec <- rowMeans(all.expr[, sample.indices3])

mean.vecs <-
    dlply(anno, .variables = c("cell.type.ch1"),
          .fun = function(anno.cell) {
              expr.sub <- expr.log.scale[, colnames(expr.log.scale) %in% anno.cell$geo_accession]
              rowMeans(expr.sub)
          })
tmp <- ldply(mean.vecs)
rownames(tmp) <- tmp$cell.type.ch1
tmp <- tmp[, !(colnames(tmp) == "cell.type.ch1")]
cell.type.mean.expr <- t(tmp)

for(i in 1:ncol(cell.type.mean.expr)) {
    nm <- colnames(cell.type.mean.expr)[i]
    vals <- cell.type.mean.expr[,i]
    g <- plot.col(vals)

    gt <- ggplotGrob(g)
    prefix <- make.names(nm)
    prefix <- gsub(prefix, pattern="\\.", replacement="")
##    png(paste0(prefix, ".png"),10,50)
    grid.draw(gtable::gtable_filter(gt, "panel"))
    ggsave(paste0(prefix, ".png"),width=10,height=50,units="cm")
##    d <- dev.off()
}

g <- plot.col(admixture.vec)
gt <- ggplotGrob(g)
prefix <- "admixture"
##png(paste0(prefix, ".png"),10,50)
##png(paste0(prefix, ".png"))
grid.draw(gtable::gtable_filter(gt, "panel"))
ggsave(paste0(prefix, ".png"),width=10,height=50,units="cm")
##d <- dev.off()

g <- plot.col(admixture2.vec)
gt <- ggplotGrob(g)
prefix <- "admixture2"
##png(paste0(prefix, ".png"),10,50)
##png(paste0(prefix, ".png"))
grid.draw(gtable::gtable_filter(gt, "panel"))
ggsave(paste0(prefix, ".png"),width=10,height=50,units="cm")
##d <- dev.off()

g <- plot.col(admixture3.vec)
gt <- ggplotGrob(g)
prefix <- "admixture3"
##png(paste0(prefix, ".png"),10,50)
##png(paste0(prefix, ".png"))
grid.draw(gtable::gtable_filter(gt, "panel"))
ggsave(paste0(prefix, ".png"),width=10,height=50,units="cm")
##d <- dev.off()

g <- plot.col(admixture.vec, remove.title=TRUE, remove.axis.labels=FALSE)
g <- g + xlab("Admixture") + labs(y=NULL)
g <- g + theme(axis.title.x = element_text(size = 25))
prefix <- "admixture-title"
##png(paste0(prefix, ".png"))
g
ggsave(paste0(prefix, ".png"),width=10,height=50,units="cm")
## d <- dev.off()


g <- plot.col(noise.vec)
gt <- ggplotGrob(g)
##png(paste0(prefix, ".png"),10,50)
grid.draw(gtable::gtable_filter(gt, "panel"))
prefix <- "noise"
ggsave(paste0(prefix, ".png"),width=10,height=50,units="cm")
##d <- dev.off()

tmp <- subset(proportions, cell.type %in% cell.type.levels)
tmp$cell.type <- factor(tmp$cell.type, levels = cell.type.levels)
g <- ggplot()
g <- g + geom_col(data = tmp,aes(x=cell.type, y=proportions), width=0.5, fill="blue")
g <- g + coord_flip()
g <- g + theme_classic()
g <- g + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
g <- g + scale_y_continuous(limits=c(0,max(tmp$proportions)),expand=c(0,0))
##g <- g + theme(text = element_text(size = 25))
g <- g + theme(text = element_text(size = 35))
g <- g + ylab("Proportion") + xlab("")
g
prefix <- "proportions"
ggsave(paste0(prefix, ".png"), width=10,height=5)

b.cell.genes <- subset(cibersort.genes, cell_type == "B cells")$Hugo
## Cheat here to get some diversity in the scores
admixtures <- list("A" = admixture1.orig.vec,
                   "C" = all.expr[, 41],
                   ## "B" = admixture2.orig.vec,
                   "B" = all.expr[, 109]
                   ## "C" = admixture3.orig.vec
                   )
names(admixtures[["B"]]) <- rownames(all.expr)
names(admixtures[["C"]]) <- rownames(all.expr)
b.cell.scores <-
    llply(admixtures,
          .fun = function(vec) {
              cols <- b.cell.genes[b.cell.genes %in% names(vec)]
              mean(log2(vec[cols]))
          })
llply(1:ncol(all.expr),
      .fun = function(i) {
          rows <- b.cell.genes[b.cell.genes %in% rownames(all.expr)]
          mean(log2(all.expr[rows,i]))
      })

tmp <- ldply(b.cell.scores)
colnames(tmp) <- c("Admixture", "b.cell.score")
tmp$Admixture <- factor(tmp$Admixture, levels=rev(c("A","B","C")))

g <- ggplot()
g <- g + geom_col(data = tmp,aes(x=Admixture, y=b.cell.score), width=0.5, fill="darkgreen")
g <- g + coord_flip()
g <- g + theme_classic()
## g <- g + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
g <- g + scale_y_continuous(limits=c(0,max(tmp$b.cell.score)),expand=c(0,0))
g <- g + theme(text = element_text(size = 35))
g <- g + ylab("B cell score") + xlab("Admixture")
g
prefix <- "b-cell-score"
ggsave(paste0(prefix, ".png"), width=10,height=5)

create.admixture <- function(all.expr.log.scale, seed = 1234) {
  set.seed(seed)
  n.col <- 22
  tmp <- llply(1:n.col,
               .fun = function(i) {
                   sample.indices <- sample.int(ncol(all.expr.log.scale),
                                                size = ncol(all.expr.log.scale), replace = TRUE)
                   admixture.vec <- rowMeans(all.expr.log.scale[, sample.indices])
               })
  tmp <- do.call("cbind", tmp)
  tmp.long <- reshape2::melt(as.matrix(tmp))
  tmp.long
}

prefix <- "admixture-matrix-1"
png(paste0(prefix, ".png"))
tmp.long <- create.admixture(all.expr.log.scale, seed = 1234)
g <- ggplot(data=tmp.long, aes(x=Var2, y=Var1, fill=value)) 
g <- g + geom_tile()
g <- g + scale_fill_gradient2(low = low.col, mid = mid.col, high = high.col)
g <- remove.axes(g, remove.title=FALSE, remove.axis.labels=FALSE)
g <- g + ylab("Genes") + xlab("Admixtures")
## g <- g + scale_x_discrete(position="top")
g <- g + theme(text = element_text(size = 35),
               axis.title.x = element_text(size = 35),
               axis.title.y = element_text(size = 35))
## g <- g + ggtitle("Admixtures") + theme(plot.title = element_text(hjust = 0.5))
g
print(g)
##grid.draw(gtable::gtable_filter(gt, "panel"))
d <- dev.off()

prefix <- "admixture-matrix-2"
png(paste0(prefix, ".png"))
tmp.long <- create.admixture(all.expr.log.scale, seed = 4321)
g <- ggplot(data=tmp.long, aes(x=Var2, y=Var1, fill=value)) 
g <- g + geom_tile()
g <- g + scale_fill_gradient2(low = low.col, mid = mid.col, high = high.col)
g <- remove.axes(g, remove.title=FALSE, remove.axis.labels=FALSE)
g <- g + ylab("Genes") + xlab("Admixtures")
## g <- g + scale_x_discrete(position="top")
g <- g + theme(text = element_text(size = 35),
               axis.title.x = element_text(size = 35),
               axis.title.y = element_text(size = 35))
## g <- g + ggtitle("Admixtures") + theme(plot.title = element_text(hjust = 0.5))
g
print(g)
##grid.draw(gtable::gtable_filter(gt, "panel"))
d <- dev.off()


tmp <- mat
tmp[tmp < 9] <- 0
non.zero <- unlist(apply(tmp, 1, function(row) length(which(row > 0))))
tmp <- tmp[non.zero == 1,]
hc <- hclust(dist(tmp))
tmp <- tmp[hc$order,]
hc <- hclust(dist(t(tmp)))
tmp <- tmp[, hc$order]
tmp <- tmp[, c(3,4,2,5,1)]
n.zero <- length(which(tmp == 0))
tmp[tmp == 0] <- rnorm(n.zero, mean = 5)
tmp.long <- reshape2::melt(as.matrix(tmp))
prefix <- "abbas-purified"
png(paste0(prefix, ".png"))
g <- ggplot(data=tmp.long, aes(x=Var2, y=Var1, fill=value)) 
g <- g + geom_tile()
print(c(low.col, mid.col, high.col))
g <- g + scale_fill_gradient2(low = low.col, mid = mid.col, high = high.col, midpoint = median(tmp))
g <- remove.axes(g, remove.title=FALSE, remove.axis.labels=FALSE)
g <- g + ylab("Genes") + xlab("Purified Samples")
## g <- g + scale_x_discrete(position="top")
g <- g + theme(text = element_text(size = 35),
               axis.title.x = element_text(size = 35),
               axis.title.y = element_text(size = 35))
## g <- g + ggtitle("Admixtures") + theme(plot.title = element_text(hjust = 0.5))
g
print(g)
##grid.draw(gtable::gtable_filter(gt, "panel"))
d <- dev.off()

## TODO

## write text


## email group about: time and poster

## email group about admixtures
