suppressPackageStartupMessages(library("pacman"))
suppressPackageStartupMessages(p_load("synapser"))
suppressPackageStartupMessages(p_load("ggrepel"))
suppressPackageStartupMessages(p_load("plyr"))
suppressPackageStartupMessages(p_load("dplyr"))
suppressPackageStartupMessages(p_load("reshape2"))

suppressPackageStartupMessages(p_load("immunedeconv"))
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

synLogin()

## Load the TPM validation data
##synId <- "syn21574299"
##obj <- synGet(synId, downloadFile = TRUE)
##cpm.expr <- read.table(obj$path, sep = "\t", header = TRUE)

synId <- "syn21576632"
obj <- synGet(synId, downloadFile = TRUE)
cpm.expr <- read.table(obj$path, sep = ",", header = TRUE)


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

purified.expr <- cpm.expr[, purified.samples]
purified.expr <- rename.samples(purified.expr)

admixture.expr <- cpm.expr[, !(colnames(cpm.expr) %in% purified.samples)]

plot.deconv.scores <- function(res, order.by = NULL) {
  long.df <- melt(res)
  colnames(long.df) <- c("cell.type", "sample", "score")
  new_order <- NULL
  g <- NULL
  if(!is.null(order.by)) {
    vec <- res[order.by, ]
    o <- order(as.numeric(vec))
    vec <- vec[o]
    levels <- names(vec)
    long.df$sample <- factor(long.df$sample, levels = levels)
    g <- ggplot(data = long.df,
                aes(x = sample, y = score, cell.type = cell.type))
  } else {
    new_order <- long.df %>% 
      group_by(cell.type) %>% 
      do(data_frame(al=levels(reorder(interaction(.$cell.type, .$sample, drop=TRUE), .$score)))) %>% 
        pull(al)
    long.df <- long.df %>% 
        mutate(al=factor(interaction(.$cell.type, .$sample), levels=new_order))
    g <- ggplot(data = long.df,
                aes(x = al, y = score, cell.type = cell.type))
  }
  g <- g + facet_wrap(~ cell.type, scales = "free")
  g <- g + theme(axis.text.x = element_text(angle = 90))
  g <- g + geom_point()
  if(is.null(order.by)) {
      g <- g + scale_x_discrete(breaks= new_order, labels=gsub("^.*\\.", "", new_order))
  }
  g <- g + xlab("cell type")
  g
}

## expr.mat should be genes x samples, with genes indicated as symbols
run.deconv.methods <- function(expr.mat) {

    files <- c("/home/bwhite/LM22.txt", "/Users/Brian/Downloads/new-cibersort-code/LM22.txt")
    for(file in files) {
        if(file.exists(file)) {
            cat(paste0("Setting CIBERSORT mat to ", file, "\n"))
            immunedeconv::set_cibersort_mat(file)
        }
    }
    files <- c("/home/bwhite/CIBERSORT.R", "/Users/Brian/Downloads/new-cibersort-code/CIBERSORT.R")
    for(file in files) {
        if(file.exists(file)) {
            cat(paste0("Setting CIBERSORT script to ", file, "\n"))
            immunedeconv::set_cibersort_binary(file)
        }
    }
  xcell.fine.grained.cell.types <-
      c("CD4+ memory T-cells", "CD4+ naive T-cells",
        "CD8+ T-cells", "DC", "Endothelial cells",
        "Fibroblasts", "Macrophages", "Memory B-cells", "Monocytes", "naive B-cells",
        "Neutrophils", "NK cells", "Tregs")

  xcell.fine.grained.cell.types <-
      c("T cell CD4+ memory", "T cell CD4+ naive",
        "T cell CD8+", "Myeloid dendritic cell", "Endothelial cell",
        "Cancer associated fibroblast", "Macrophage", "B cell memory",
        "Monocyte", "B cell naive", "Neutrophil", "NK cell", "T cell regulatory (Tregs)")

  expected_cell_types = xcell.fine.grained.cell.types

  arrays <- FALSE
  tumor <- FALSE
  scale_mrna <- FALSE
  ret.list <- list()
  ret.list[["mcp"]] <- MCPcounter.estimate(expr.mat, featuresType = "HUGO_symbols")
  ret.list[["cibersort"]] <- deconvolute_cibersort(expr.mat, arrays = arrays, absmean = TRUE)
  ret.list[["xcell"]] <-
      deconvolute_xcell(expr.mat, arrays = arrays, expected_cell_types = expected_cell_types)
  ret.list[["epic"]] <- deconvolute_epic(expr.mat, tumor = tumor, scale_mrna = scale_mrna)
  indications <- rep("coad", ncol(expr.mat))
  ret.list[["timer"]] <- deconvolute_timer(expr.mat, indications = indications)
  ret.list[["quantiseq"]] <-
      deconvolute_quantiseq(expr.mat, tumor = tumor, arrays = arrays, scale_mrna = scale_mrna)
  return(ret.list)
}

plot.all.purified.results <- function(res.list, suffix) {
  nms <- names(res.list)
  names(nms) <- nms
  deconv.plots <-
      llply(nms,
            .fun = function(nm) {
                g <- plot.deconv.scores(res.list[[nm]])
                g <- g + ggtitle(nm)
                g
            })
  
  nm <- "cibersort"
  mid.pt <- floor(nrow(res.list[[nm]])/2)
  rows <- 1:mid.pt
  g <- plot.deconv.scores(res.list[[nm]][rows, ])
  g <- g + ggtitle(nm)
  deconv.plots[["cibersort-1"]] <- g
  
  rows <- (mid.pt+1):nrow(res.list[[nm]])
  g <- plot.deconv.scores(res.list[[nm]][rows, ])
  g <- g + ggtitle(nm)
  deconv.plots[["cibersort-2"]] <- g
  
  nms <- names(deconv.plots)
  names(nms) <- nms
  
  ## Save the plots
  l_ply(nms,
        .fun = function(nm) {
            g <- deconv.plots[[nm]]
            pdf(paste0(nm, suffix, ".pdf"), width = 14)
            print(g)
            d <- dev.off()
        })
}

purified.deconv <- run.deconv.methods(purified.expr)
plot.all.purified.results(purified.deconv, suffix = "-purified-deconv")

suppressPackageStartupMessages(p_load(AnnotationHub))
suppressPackageStartupMessages(p_load(ensembldb))
ah <- AnnotationHub()
flag <- (ah$species == "Homo sapiens") & (ah$genome == "GRCh38") & (ah$dataprovider == "Ensembl") & (ah$rdataclass == "EnsDb")
ah2 <- ah[flag, ]
## as.data.frame(mcols(ah2))[1:10,c("title"),drop=FALSE]
edb <- ah2[["AH73881"]]

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

purified.pc.deconv <- run.deconv.methods(purified.pc.expr)
plot.all.purified.results(purified.pc.deconv, suffix = "-purified-pc-deconv")



## Our original admixture specification includes the vendor for each sample
synId <- "syn21577258"
obj <- synGet(synId, downloadFile = TRUE)
vendors1 <- read.xlsx(obj$path, sheet = 1)[3,,drop=F]
vendors2 <- read.xlsx(obj$path, sheet = 2)[3,,drop=F]
vendors1 <- vendors1[,-1,drop=F]
vendors2 <- vendors2[,-1,drop=F]
exclude.cols <- c("Tregs", "Endothelial_cells", "breast", "CRC", "Fibroblasts")
colnames(vendors1) <-
    unlist(lapply(colnames(vendors1),
                  function(str) ifelse(str %in% exclude.cols,
                                       str,
                                       paste0(str, "_1"))))
colnames(vendors2) <-
    unlist(lapply(colnames(vendors2),
                  function(str) ifelse(str %in% exclude.cols,
                                       str,
                                       paste0(str, "_2"))))
vendors <- cbind(vendors1, vendors2[, !(colnames(vendors2) %in% exclude.cols), drop=F])
vendors <- t(rename.samples(vendors))
rownames(vendors)[grepl(rownames(vendors), pattern="breast")] <- "Breast"
colnames(vendors) <- c("vendor")
vendors <- data.frame(sample = rownames(vendors), vendors)

## Get the final admixture ratios
## Note these rows to not sum to 1, so need to renormalize
## This is the old/first "ratio" file that John Coller sent us "Final_Mixture_Ratios.xlsx"
## Instead use the subsequent file he sent "Final_Mixture_Fractions.xlsx" in which
## he appears to have just decided the rows by the row sum (which is what I do for the ratio file below).
if(FALSE) {
  synId <- "syn21577248"
  obj <- synGet(synId, downloadFile = TRUE)
  stop("Need to read _both_ sheets -- this code is not doing that. See below for fractions file")
  ratios <- read.xlsx(obj$path, sheet = 1, startRow = 2)
  ratios <- ratios[-1,]
  nms <- ratios$X1
  ratios <- ratios[,-1]
  ratios <-
      ldply(1:nrow(ratios), .fun = function(i) ratios[i,] / sum(ratios[i,]))
  rownames(ratios) <- nms
}

synId <- "syn21598638"
obj <- synGet(synId, downloadFile = TRUE)
old.col.names <- c("breast", "CRC", "Fibroblasts", "Endothelial_cells", "Dendritic_cells",
                   "Monocytes", "Macrophages", "NK_cells", "Tregs", "Naive_CD4_T_cells",
		   "Memory_CD4_T_cells", "Memory_CD8_T_cells", "Naive_B_cells")
new.col.names <- c("breast", "CRC", "fibroblasts", "endothelial.cells", "DC", "monocytes",
                   "macrophages", "NK.cells", "regulatory.T.cells", "naive.CD4.T.cells",
                   "memory.CD4.T.cells", "memory.CD8.T.cells", "naive.B.cells")

exclude.cols <- c("Tregs", "Endothelial_cells", "breast", "CRC", "Fibroblasts")

## Leave off the second sheet (note 1:1 below) since Naive_B_cells_2
## appear to have been used in admixtures, but not sequenced.
sample.ratios <-
  ldply(1:1,
        .fun = function(sheetIndex) {
                 df <- read.xlsx(obj$path, sheet = sheetIndex, startRow = 2)
                 df <- df[-1,]
                 nms <- df$X1
                 df <- df[,-1]
		 df <- df[, old.col.names]
                 colnames(df) <-
                   unlist(lapply(colnames(df),
                                 function(str) ifelse(str %in% exclude.cols,
                                                      str,
                                                      paste0(str, "_", sheetIndex))))
                 rownames(df) <- nms
                 for(col in 1:ncol(df)) { df[,col] <- as.numeric(df[,col]) }		 
                 m <- melt(as.matrix(df))
		 colnames(m) <- c("sample", "cell.type", "actual")
		 m
		 
	})
sample.ratios$cell.type <- as.character(sample.ratios$cell.type)
sample.ratios$sample <- as.character(sample.ratios$sample)

sample.ratios[sample.ratios$cell.type == "breast", "cell.type"] <- "Breast"
eps <- 10^-4
chk <- ddply(sample.ratios, .variables = c("sample"), .fun = function(df) sum(df$actual))
if(any(abs(chk$V1 - 1) > eps)) { stop("Something doesn't sum to one") }

table(sample.ratios$cell.type %in% colnames(cpm.expr))

write.table(file = "admixture-sample-fractions.tsv", sample.ratios, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

purified.mat <- cpm.expr[, purified.samples]
insilico.admixtures <-
  dlply(sample.ratios,
        .variables = c("sample"),
	.fun = function(df) {
                 cols <- as.character(df$cell.type)
		 fracs <- as.numeric(df$actual)
		 mat <- purified.mat[, cols]
		 ret <- as.matrix(mat) %*% fracs
		 colnames(ret) <- df$sample[1]
		 ret
               })
insilico.admixtures <- do.call(cbind, insilico.admixtures)
admix.names <- intersect(colnames(insilico.admixtures), colnames(admixture.expr))
insilico.admixtures <- insilico.admixtures[, admix.names]
admixture.expr <- admixture.expr[, admix.names]
colnames(insilico.admixtures) <- paste0(colnames(insilico.admixtures), "i")
mat <- cbind(insilico.admixtures, admixture.expr)

dst <- dist(t(as.matrix(mat)))
hc <- hclust(dst^2, "cen")

## stop("stop")

ratios <-
  ldply(1:2,
        .fun = function(sheetIndex) {
                 df <- read.xlsx(obj$path, sheet = sheetIndex, startRow = 2)
                 df <- df[-1,]
                 nms <- df$X1
                 df <- df[,-1]
                 df <- df[, old.col.names]
                 colnames(df) <- new.col.names
		 df$sample <- nms
		 df
	})

rownames(ratios) <- ratios$sample
ratios <- ratios[, !(colnames(ratios) == "sample")]
for(col in 1:ncol(ratios)) { ratios[,col] <- as.numeric(ratios[,col]) }

l_ply(1:nrow(ratios),
      .fun = function(i) if(abs(sum(as.numeric(ratios[i,])) - 1) > eps) { stop("Rows don't sum to one") })

melted.ratios <- melt(as.matrix(ratios))
colnames(melted.ratios) <- c("sample", "cell.type", "actual")
melted.ratios$cell.type <- as.character(melted.ratios$cell.type)

admixture.deconv <- list()
admixture.deconv[["mcp"]] <- MCPcounter.estimate(admixture.expr, featuresType = "HUGO_symbols")

m <- melt(admixture.deconv[["mcp"]])
colnames(m) <- c("cell.type", "sample", "predicted")
m$cell.type <- as.character(m$cell.type)

m <- m[grepl(m$sample, pattern="BM"),]

gt.ct <- "Naive_B_cells"
pred.ct <- "B lineage"
gt.ct <- "Memory_CD8_T_cells"
pred.ct <- "CD8 T cells"
gt.ct <- "fibroblasts"
pred.ct <- "Fibroblasts"
mer <- merge(subset(melted.ratios, cell.type == gt.ct)[, c("sample", "actual")],
             subset(m, cell.type == pred.ct)[, c("sample", "predicted")])

plot(mer$actual, mer$predicted)
cor(mer$actual, mer$predicted)


prediction.files <-
    list("mcp" = list("coarse" = "syn21576653", "fine" = "syn21576652"),
         "epic" = list("coarse" = "syn21576657", "fine" = "syn21576656"),
         "cibersort" = list("coarse" = "syn21576643", "fine" = "syn21576642"),
         "xcell" = list("coarse" = "syn21576865", "fine" = "syn21576864"))

## epic fine
synId <- "syn21576656"

## xcell fine
synId <- "syn21576864"
meth <- "cibersort"
meth <- "epic"
meth <- "xcell"
sc <- "fine"
synId <- prediction.files[[meth]][[sc]]
obj <- synGet(synId, downloadFile = TRUE)
tbl <- read.table(obj$path, sep = ",", header = TRUE)
colnames(tbl) <- c("dataset.name", "sample", "cell.type", "prediction")

mer <- merge(melted.ratios, tbl)
g <- ggplot(data = mer, aes(x = actual, y = prediction))
g <- g + geom_point()
g <- g + facet_wrap(~ cell.type, scale = "free_y")

xcell.deconv.genes <-
    sort(unique(unlist(lapply(xCell.data$signatures, function(x) x@geneIds))))

lm22 <- NULL
files <- c("/home/bwhite/LM22.txt", "/Users/Brian/Downloads/new-cibersort-code/LM22.txt")
for(file in files) {
    if(file.exists(file)) {
        lm22 <- read.table(file, sep="\t", header=TRUE)
    }
}

cibersort.deconv.genes <- as.character(lm22$Gene.symbol)
mcp.deconv.genes <-
    sort(unique(
        read.table(curl:::curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"), 
                   sep = "\t", stringsAsFactors = FALSE, header = TRUE, 
                   colClasses = "character", check.names = FALSE)$`HUGO symbols`))

epic.deconv.genes <- sort(unique(c(TRef$sigGenes, BRef$sigGenes)))

signame <- "TIL10"
sig.mat.file <- system.file("extdata", "quantiseq",
                            paste0(signame, 
                                   "_signature.txt"), package = "immunedeconv", mustWork = TRUE)
df <- read.table(sig.mat.file, sep="\t", header=TRUE)
quantiseq.deconv.genes <- as.character(unique(df$ID))

all.zero <- unlist(apply(purified.expr, 1, function(row) all(row == 0)))
purified.expr <- purified.expr[!all.zero,]

means <- unlist(apply(purified.expr, 1, mean))
vars <- unlist(apply(purified.expr, 1, var))
cov <- vars / means

cov.genes <- names(tail(cov[order(cov)],n=100))

deconv.genes <-
    sort(unique(c(cibersort.deconv.genes, mcp.deconv.genes, epic.deconv.genes, quantiseq.deconv.genes)))
deconv.genes <- deconv.genes


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

plot.pca <- function(mat, colour = NULL, pcs = c(1,2), show.labels = TRUE) {
    pc <- mat
    if(class(mat) != "prcomp") {
        pc <- prcomp(t(mat), scale = TRUE, center = TRUE)
    }
    id <- rownames(pc$x)
    df <- data.frame(id = id, x = pc$x[,pcs[1]], y = pc$x[,pcs[2]])
    if(!is.null(colour)) {
        colour.df <- data.frame(id = names(colour), colour = unname(colour))
        df <- merge(df, colour.df)
    }
    df$id <- gsub(df$id, pattern = "_1", replacement = "")
    df$id <- gsub(df$id, pattern = "_2", replacement = "")    
    g <- ggplot(data = df, aes(x = x, y = y))
    if(!is.null(colour)) {
        g <- g + geom_point(aes(colour = colour))
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

colour <- vendors$vendor
names(colour) <- vendors$sample

purified.expr.qn <- normalize.quantiles(as.matrix(purified.expr))
rownames(purified.expr.qn) <- rownames(purified.expr)
colnames(purified.expr.qn) <- colnames(purified.expr)
all.zero <- unlist(apply(purified.expr.qn, 1, function(row) all(row == row[1])))
purified.expr.qn <- purified.expr.qn[!all.zero,]

deconv.purified.expr.qn <- purified.expr.qn[rownames(purified.expr.qn) %in% deconv.genes,]

pca.plots <- list()
g1 <- plot.pca(as.matrix(purified.expr), colour = colour, pcs = c(1,2))
g1 <- g1 + ggtitle("All genes; Uncorrected")
pca.plots[["all-uncorrected-g1"]] <- g1
g2 <- plot.pca(as.matrix(purified.expr), colour = colour, pcs = c(3,4))
g2 <- g2 + ggtitle("All genes; Uncorrected")
pca.plots[["all-uncorrected-g2"]] <- g2

g1 <- plot.pca(as.matrix(deconv.purified.expr), colour = colour, pcs = c(1,2))
g1 <- g1 + ggtitle("Deconvolution genes; Uncorrected")
pca.plots[["deconv-uncorrected-g1"]] <- g1
g2 <- plot.pca(as.matrix(deconv.purified.expr), colour = colour, pcs = c(3,4))
g2 <- g2 + ggtitle("Deconvolution genes; Uncorrected")
pca.plots[["deconv-uncorrected-g2"]] <- g2

g1 <- plot.pca(as.matrix(purified.expr.qn), colour = colour, pcs = c(1,2))
g1 <- g1 + ggtitle("All genes; Quantile normalized")
pca.plots[["all-qn-g1"]] <- g1
g2 <- plot.pca(as.matrix(purified.expr.qn), colour = colour, pcs = c(3,4))
g2 <- g2 + ggtitle("All genes; Quantile normalized")
pca.plots[["all-qn-g2"]] <- g2

g1 <- plot.pca(as.matrix(deconv.purified.expr.qn), colour = colour, pcs = c(1,2))
g1 <- g1 + ggtitle("Deconvolution genes; Quantile normalized")
pca.plots[["deconv-qn-g1"]] <- g1
g2 <- plot.pca(as.matrix(deconv.purified.expr.qn), colour = colour, pcs = c(3,4))
g2 <- g2 + ggtitle("Deconvolution genes; Quantile normalized")
pca.plots[["deconv-qn-g2"]] <- g2

batch <- unname(vendors[colnames(purified.expr),"vendor"])
batch <- as.numeric(unname(vendors[colnames(purified.expr),"vendor"]))
purified.expr.combat <- ComBat(as.matrix(purified.expr), batch = batch)
deconv.purified.expr.combat <- purified.expr.combat[rownames(purified.expr.combat) %in% deconv.genes,]

g1 <- plot.pca(as.matrix(purified.expr.combat), colour = colour, pcs = c(1,2))
g1 <- g1 + ggtitle("All genes; ComBat corrected")
pca.plots[["all-combat-g1"]] <- g1
g2 <- plot.pca(as.matrix(purified.expr.combat), colour = colour, pcs = c(3,4))
g2 <- g2 + ggtitle("All genes; ComBat corrected")
pca.plots[["all-combat-g2"]] <- g2

g1 <- plot.pca(as.matrix(deconv.purified.expr.combat), colour = colour, pcs = c(1,2))
g1 <- g1 + ggtitle("Deconvolution genes; ComBat corrected")
pca.plots[["deconv-combat-g1"]] <- g1
g2 <- plot.pca(as.matrix(deconv.purified.expr.combat), colour = colour, pcs = c(3,4))
g2 <- g2 + ggtitle("Deconvolution genes; ComBat corrected")
pca.plots[["deconv-combat-g2"]] <- g2

pdf("uncorrected-pca.pdf", width = 14)
grid.arrange(pca.plots[["all-uncorrected-g1"]], pca.plots[["deconv-uncorrected-g1"]], nrow = 1)
d <- dev.off()

pdf("qn-pca.pdf", width = 14)
grid.arrange(pca.plots[["all-qn-g1"]], pca.plots[["deconv-qn-g1"]], nrow = 1)
d <- dev.off()

pdf("combat-pca.pdf", width = 14)
grid.arrange(pca.plots[["all-combat-g1"]], pca.plots[["deconv-combat-g1"]], nrow = 1)
d <- dev.off()

grab_grob <- function(){
  grid.echo()
  grid.grab()
}

##main <- "Blah genes; uncorrected"
####par(mfrow=c(1,2))
##corrplot(cor(m), order="hclust", main = main)
##g1 <- grab_grob()
##heatmap(as.matrix(t(scale(t(m)))), scale="none", main = main)
##g2 <- grab_grob()
####heatmap.2(as.matrix(t(scale(t(m)))), scale="none", main = main)
####Heatmap(as.matrix(t(scale(t(m)))))
##pdf("uncorrected-heatmap-corrplot.pdf", width = 14)
##grid.arrange(g1, g2, ncol = 2)
##d <- dev.off()

m <- deconv.purified.expr
main <- "Deconvolution genes; uncorrected"
pdf("uncorrected-corrplot.pdf")
corrplot(cor(m), order="hclust", main = main)
d <- dev.off()

batch <- unname(vendors[colnames(m),"vendor"])
n <- length(unique(batch))
dat.col <- data.frame(batch = unique(batch),
                      color = brewer_pal()(n))

rownames(dat.col) <- dat.col$batch
col.side.colors <- as.character(dat.col[vendors[colnames(m), "vendor"], "color"])
pdf("uncorrected-heatmap.pdf")
hc <- hclust(dist(t(scale(as.matrix(t(m))))), "ave")
null.row.labels <- rep("", nrow(m))
heatmap(as.matrix(m), scale="row", ColSideColors = col.side.colors, main = main, Rowv = as.dendrogram(hc),
        labRow = null.row.labels)
d <- dev.off()


m <- deconv.purified.expr.qn
main <- "Deconvolution genes; quantile normalized"
pdf("qn-corrplot.pdf")
corrplot(cor(m), order="hclust", main = main)
d <- dev.off()

batch <- unname(vendors[colnames(m),"vendor"])
n <- length(unique(batch))
dat.col <- data.frame(batch = unique(batch),
                      color = brewer_pal()(n))

rownames(dat.col) <- dat.col$batch
col.side.colors <- as.character(dat.col[vendors[colnames(m), "vendor"], "color"])
pdf("qn-heatmap.pdf")
hc <- hclust(dist(t(scale(as.matrix(t(m))))), "ave")
null.row.labels <- rep("", nrow(m))
heatmap(as.matrix(m), scale="row", ColSideColors = col.side.colors, main = main, Rowv = as.dendrogram(hc),
        labRow = null.row.labels)
d <- dev.off()

batch.colour = colnames(admixture.expr)
names(batch.colour) <- batch.colour
batch.colour <- gsub(batch.colour, pattern = "\\d", replacement="")

all.zero <- unlist(apply(admixture.expr, 1, function(row) all(row == row[1])))
admixture.expr <- admixture.expr[!all.zero,]

g1 <- plot.pca(admixture.expr, colour = batch.colour, pcs = c(1,2), show.labels = FALSE)

immune.markers <-
  c("CD4", "CD8A", "CD8B", "CD14", "GZMB", "NCAM1", "CD2", "CD3D", "CD3E", "CD3G", "CD5","PTPRC","FCGR3A",
    "CD19","MS4A1","CR2","CD28","FCGR2A","FCGR2B", "FCGR2C", "CDR1","CD40","KLRB1", "FOXP3","IL2RA","ICAM1",
    "CD79A","CD79B","CD86","FAS", "TNFSF13", "TNFRSF13C", "TNFSF13")

## Our original admixture specification includes the vendor for each sample
synId <- "syn21577258"
obj <- synGet(synId, downloadFile = TRUE)
bm1 <- read.xlsx(obj$path, sheet = 1)
bm1 <- bm1[-c(1:3),]
bm2 <- read.xlsx(obj$path, sheet = 2)
bm2 <- bm2[-c(1:3),]
rm1 <- read.xlsx(obj$path, sheet = 3)
rm1 <- rm1[-c(1:3),]
rm2 <- read.xlsx(obj$path, sheet = 4)
rm2 <- rm2[-c(1:3),]

admx.batch.colour <- colnames(admixture.expr)
names(admx.batch.colour) <- admx.batch.colour
flag <- admx.batch.colour %in% bm1[,1]
admx.batch.colour[flag] <- "BM1"
flag <- admx.batch.colour %in% bm2[,1]
admx.batch.colour[flag] <- "BM2"
flag <- admx.batch.colour %in% rm1[,1]
admx.batch.colour[flag] <- "RM1"
flag <- admx.batch.colour %in% rm2[,1]
admx.batch.colour[flag] <- "RM2"

pc <- prcomp(t(admixture.expr), scale = TRUE, center = TRUE)
df <- data.frame(id = rownames(pc$x), x = pc$x[, 1], y = pc$x[, 2])
lst <- list(bm1, bm2, rm1, rm2)
df2 <-
    do.call(rbind, lapply(lst, function(entry) data.frame(id = entry[,1], label = ifelse(entry$CRC > 0, "CRC", ""))))
df <- merge(df, df2)
g2 <- plot.pca(admixture.expr, colour = admx.batch.colour, pcs = c(1,2), show.labels = FALSE)
## g2 <- g2 + geom_text_repel(aes(x = df$x, y = df$y, label = df$label))
g2 <- g2 + geom_text_repel(aes(x = df$x, y = df$y, label = df$label))

pdf("admixture-pca.pdf")
print(g2)
d <- dev.off()

df <-
    do.call(rbind, lapply(lst,
                          function(entry) {
                              label = ifelse(entry$CRC > 0, "CRC",
                                             ifelse(entry$breast > 0, "BRCA", ""))                              
                              data.frame(id = entry[,1], tumor.type = label)
                          }))
df$batch <- "empty"
flag <- df$id %in% bm1[,1]
df$batch[flag] <- "BM1"
flag <- df$id %in% bm2[,1]
df$batch[flag] <- "BM2"
flag <- df$id %in% rm1[,1]
df$batch[flag] <- "RM1"
flag <- df$id %in% rm2[,1]
df$batch[flag] <- "RM2"

df$mixture.type <- "empty"
flag <- ( df$id %in% bm1[,1] ) | ( df$id %in% bm2[,1] )
df$mixture.type[flag] <- "BM"
flag <- ( df$id %in% rm1[,1] ) | ( df$id %in% rm2[,1] )
df$mixture.type[flag] <- "RM"

df$dataset <- "empty"
flag <- ( df$tumor.type == "BRCA" ) & ( df$mixture.type == "BM" )
df$dataset[flag] <- "DS1"
flag <- ( df$tumor.type == "CRC" ) & ( df$mixture.type == "BM" )
df$dataset[flag] <- "DS2"
flag <- ( df$tumor.type == "BRCA" ) & ( df$mixture.type == "RM" )
df$dataset[flag] <- "DS3"
flag <- ( df$tumor.type == "CRC" ) & ( df$mixture.type == "RM" )
df$dataset[flag] <- "DS4"

synId <- "syn21576632"
obj <- synGet(synId, downloadFile = TRUE)
cpm.expr <- read.table(obj$path, sep = ",", header = TRUE)

samples <- colnames(cpm.expr)
## Exclude the admixtures samples
flag <- grepl(samples, pattern="RM")
samples <- samples[!flag]
flag <- grepl(samples, pattern="BM")
samples <- samples[!flag]
samples <- samples[!(samples == "Gene")]
names(samples) <- samples

for(col in colnames(df)) { df[, col] <- as.character(df[, col]) }
df <- rbind(data.frame(id = unname(samples), tumor.type = "CRC", batch = "P", mixture.type = "P", dataset = "DS5", stringsAsFactors = FALSE), df)

write.table(file = "in-vitro-admixture-sample-datasets.tsv", df, sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE)


## Map populations to samples
populations <- samples
populations <- gsub(populations, pattern="_1", replacement="")
populations <- gsub(populations, pattern="_2", replacement="")
pop.df <- data.frame(population = unname(populations), sample = names(populations))

## Translate the population names to those we use in the challenge
fine.grained.map <- list(
    "Breast" = "Breast",
    "CRC" = "CRC",
    "Dendritic_cells" = "myeloid.dendritic.cells",
    "Endothelial_cells" = "endothelial.cells",
    "Fibroblasts" = "fibroblasts",
    "Macrophages" = "macrophages",
    "Memory_CD4_T_cells" = "memory.CD4.T.cells",
    "Memory_CD8_T_cells" = "memory.CD8.T.cells",
    "Monocytes" = "monocytes",
    "NK_cells" = "NK.cells",
    "Naive_B_cells" = "naive.B.cells",
    "Naive_CD4_T_cells" = "naive.CD4.T.cells",
    "Naive_CD8_T_cells" = "naive.CD8.T.cells",
    "Neutrophils" = "neutrophils",
    "Tregs" = "regulatory.T.cells")

fine.grained.map.df <- data.frame(population = names(fine.grained.map), challenge.population = unname(unlist(fine.grained.map)))
fine.grained.pop.df <- merge(pop.df, fine.grained.map.df)

ddply(melted.ratios, .variables = c("sample"), .fun = function(df) sum(df$actual))
gold.standard.tbl <- melted.ratios
colnames(gold.standard.tbl) <- c("sample.id", "cell.type", "measured")
gold.standard.tbl <- merge(gold.standard.tbl, unique(df[, c("id", "dataset")]), by.x = c("sample.id"), by.y = c("id"))

gold.standard.tbl <- gold.standard.tbl[, c("dataset", "sample.id", "cell.type", "measured")]
colnames(gold.standard.tbl) <- c("dataset.name", "sample.id", "cell.type", "measured")
gold.standard.tbl <- subset(gold.standard.tbl, !(cell.type %in% c("breast", "CRC")))
flag <- gold.standard.tbl$cell.type == "DC"
gold.standard.tbl$cell.type[flag] <- c("myeloid.dendritic.cells")

## Add purified samples to the gold standard
tmp <- fine.grained.pop.df
tmp$val <- 1
tmp <- tmp[, c("sample", "challenge.population", "val")]
tmp2 <- melt(acast(tmp, sample ~ challenge.population, fill = 0))
colnames(tmp2) <- c("sample.id", "cell.type", "measured")
tmp2 <- subset(tmp2, (cell.type != "Breast") & (cell.type != "CRC"))
tmp2 <- cbind(dataset.name = "DS5", tmp2)
for(col in c("dataset.name", "sample.id", "cell.type")) {
    gold.standard.tbl[, col] <- as.character(gold.standard.tbl[, col])
    tmp2[, col] <- as.character(tmp2[, col])
}
print(dim(gold.standard.tbl))
print(dim(tmp2))
gold.standard.tbl <- rbind(gold.standard.tbl, tmp2)
print(dim(gold.standard.tbl))


## fine-grained
# memory.B.cells (missing)
# naive.B.cells
# memory.CD4.T.cells
# naive.CD4.T.cells
# regulatory.T.cells
# memory.CD8.T.cells
# naive.CD8.T.cells
# NK.cells
# neutrophils (missing)
# monocytes
# myeloid.dendritic.cells
# macrophages
# fibroblasts
# endothelial.cells
fine.grained.gold.standard <- gold.standard.tbl



## coarse-grained
# B.cells
# CD4.T.cells
# CD8.T.cells
# NK.cells
# neutrophils (missing)
# monocytic.lineage
# fibroblasts
# endothelial.cells
coarse.cell.types <-
  c("B.cells", "CD4.T.cells", "CD8.T.cells", "NK.cells", "neutrophils", "monocytic.lineage",
    "fibroblasts", "endothelial.cells")

tmp <- acast(gold.standard.tbl , sample.id ~ cell.type, value.var = "measured", fill = 0)

## B.cells = memory.B.cells + naive.B.cells
flag <- colnames(tmp) %in% c("memory.B.cells", "naive.B.cells")
tmp <- cbind(tmp, B.cells = as.numeric(unlist(apply(tmp[, flag, drop=F], 1, function(row) sum(row)))))

## CD4.T.cells = memory.CD4.T.cells + naive.CD4.T.cells + regulatory.T.cells
flag <- colnames(tmp) %in% c("memory.CD4.T.cells", "naive.CD4.T.cells", "regulatory.T.cells")
tmp <- cbind(tmp, CD4.T.cells = as.numeric(unlist(apply(tmp[, flag, drop=F], 1, function(row) sum(row)))))

## CD8.T.cells = memory.CD8.T.cells + naive.CD8.T.cells
flag <- colnames(tmp) %in% c("memory.CD8.T.cells", "naive.CD8.T.cells")
tmp <- cbind(tmp, CD8.T.cells = as.numeric(unlist(apply(tmp[, flag, drop=F], 1, function(row) sum(row)))))

## monocytic.lineage = 
flag <- colnames(tmp) %in% c("myeloid.dendritic.cells", "macrophages", "monocytes")
tmp <- cbind(tmp, monocytic.lineage = as.numeric(unlist(apply(tmp[, flag, drop=F], 1, function(row) sum(row)))))

tmp <- tmp[, colnames(tmp) %in% coarse.cell.types]

coarse.grained.gold.standard <- melt(tmp)
colnames(coarse.grained.gold.standard) <- c("sample.id", "cell.type", "measured")
coarse.grained.gold.standard <-
  merge(coarse.grained.gold.standard, unique(df[, c("id", "dataset")]), by.x = c("sample.id"), by.y = c("id"))
colnames(coarse.grained.gold.standard) <- c("sample.id", "cell.type", "measured", "dataset.name")
coarse.grained.gold.standard <- coarse.grained.gold.standard[, c("dataset.name", "sample.id", "cell.type", "measured")]



write.table(file = "val_coarse.csv", coarse.grained.gold.standard, sep=",", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(file = "val_fine.csv", fine.grained.gold.standard, sep=",", row.names = FALSE, col.names = TRUE, quote = FALSE)
