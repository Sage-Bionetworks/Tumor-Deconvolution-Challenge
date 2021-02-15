suppressPackageStartupMessages(library("pacman"))
suppressPackageStartupMessages(p_load("synapser"))
suppressPackageStartupMessages(p_load("plyr"))
suppressPackageStartupMessages(p_load("dplyr"))
suppressPackageStartupMessages(p_load("reshape2"))
suppressPackageStartupMessages(p_load("ggplot2"))

synLogin()

synIds <-
    list("csx" = "syn22780960",
         "gt" = "syn22907846")

read.file.from.synapse <- function(synId, sep=",") {
    obj <- synGet(synId, downloadFile = TRUE)
    file <- obj$path
    read.table(file, sep=sep, header=TRUE, stringsAsFactors = FALSE)
}

tbls <- list()
tbls[["csx"]] <- read.file.from.synapse(synIds[["csx"]], sep="\t")
tbls[["gt"]] <- read.file.from.synapse(synIds[["gt"]], sep=",")

## Read in sample annotations
anno <- read.file.from.synapse("syn22999327", sep=",")

extract.id <- function(str) {
    ret <- gsub(str, pattern="X", replacement="")
    ret
}

tbls[["csx"]]$sample.id <- extract.id(tbls[["csx"]]$sample.id)
anno$sample.id <- extract.id(anno$id)

m <- merge(tbls[["gt"]], tbls[["csx"]], by.x = c("sample.id", "cell.type", "dataset.name"), by.y = c("sample.id", "cell.type", "dataset.name"))
m <- merge(m, anno, by.x = c("sample.id"), by.y = c("sample.id"))

create.correlation.plots <- function(m) {

    cors <- 
        dlply(m, .variables = c("cell.type", "subchallenge"),
              .fun = function(df) {
                  cell.type <- df[1, "cell.type"]
                  sc <- df[1, "subchallenge"]
                  pred <- as.numeric(df$prediction)
                  obs <- as.numeric(df$measured)
                  ct.pearson <- cor.test(pred, obs, method = "pearson")
                  ct.spearman <- cor.test(pred, obs, method = "spearman")
                  rmse <- sqrt(mean((pred-obs)^2))
                  rel.rmse <- rmse / mean(obs)
                  g <- ggplot(data = df)
                  g <- g + geom_point(aes(y = prediction, x = measured, colour = Timepoint))
                  g <- g + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
                  g <- g + geom_smooth(aes(y = prediction, x = measured), method = "lm", formula = "y ~ x")
                  g <- g + ylab("Prediction") + xlab("Measured")
                  r.str <- round(ct.pearson$estimate, digits=2)
                  r.pval <- ct.pearson$p.value
                  digits <- 2
                  r.pval <- format(r.pval, digits=digits, scientific=0)
                  rho.str <- round(ct.spearman$estimate, digits=2)
                  rmse.str <- round(rmse, digits=2)
                  rel.rmse.str <- round(rel.rmse, digits=2)              
                  rho.pval <- ct.spearman$p.value
                  rho.pval <- format(rho.pval, digits=digits, scientific=0)
                  g <- g + ggtitle(paste0(cell.type, " ", sc, ": r = ", r.str, " (p = ", r.pval, ")",
                                          " rho = ", rho.str, " (p = ", rho.pval, ")",
                                          "\nrmse = ", rmse.str, " (rel.rmse = ", rel.rmse.str, ")"))
                  ## file <- paste0("csx-ball-", make.names(cell.type), "-", sc, ".png")
                  ## png(file)
                  ## print(g)
                  ## d <- dev.off()
                  ret.list <- list(g = g,
                                   r = ct.pearson$estimate, p.pearson = ct.pearson$p.value,
                                   rho = ct.spearman$estimate, p.spearman = ct.spearman$p.value,
                                   rmse = rmse, rel.rmse = rel.rmse, range = max(obs) - min(obs))
                  return(ret.list)
                  data.frame(r = ct.pearson$estimate, p.pearson = ct.pearson$p.value,
                             rho = ct.spearman$estimate, p.spearman = ct.spearman$p.value,
                             rmse = rmse, rel.rmse = rel.rmse)
              })
    cors
}
## d <- dev.off()

cors <- create.correlation.plots(m)

all.plots <- lapply(cors, function(x) x$g)
all.cors <- ldply(cors, .fun = function(lst) {
    row <- as.data.frame(lst[names(lst) != "g"])
    row
})
cors <- all.cors

dx.cors <- create.correlation.plots(subset(m, Timepoint == "Dx"))
dx.plots <- lapply(dx.cors, function(x) x$g)
dx.cors <- ldply(dx.cors, .fun = function(lst) as.data.frame(lst[names(lst) != "g"]))

eoi.cors <- create.correlation.plots(subset(m, Timepoint == "EOI"))
eoi.plots <- lapply(eoi.cors, function(x) x$g)
eoi.cors <- ldply(eoi.cors, .fun = function(lst) as.data.frame(lst[names(lst) != "g"]))

format.tbl <- function(cors) {
    o <- order(cors$rho, decreasing = TRUE)
    cors <- cors[o,] 
    
    cols <- colnames(cors)
    cols <- cols[cols != "subchallenge"]
    cors <- unique(cors[, cols])
    
    summary.tbl <- cors
    estimate.cols <- c("r", "rho", "rmse", "rel.rmse", "range")
    for(col in estimate.cols) {
        summary.tbl[, col] <- round(summary.tbl[, col], digits=2)
    }
    
    pval.cols <- c("p.pearson", "p.spearman")
    for(col in pval.cols) {
        summary.tbl[, col] <- format(summary.tbl[, col], digits=2, scientific=0)
    }
    summary.tbl
}
    
cors.summary <- format.tbl(cors)
dx.cors.summary <- format.tbl(dx.cors)
eoi.cors.summary <- format.tbl(eoi.cors)


p_load(gridExtra)
p_load(grid)

pdf("b-all-dx-and-eoi-plots.pdf")
grid.table(cors.summary, rows = NULL)
for(g in all.plots) { print(g) }
## png("csx-ball-summary.png")
d <- dev.off()

pdf("b-all-dx-plots.pdf")
grid.table(dx.cors.summary, rows = NULL)
for(g in dx.plots) { print(g) }
## png("csx-ball-summary.png")
d <- dev.off()

pdf("b-all-eoi-plots.pdf")
grid.table(eoi.cors.summary, rows = NULL)
for(g in eoi.plots) { print(g) }
## png("csx-ball-summary.png")
d <- dev.off()
