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

extract.id <- function(str) {
    ret <- gsub(str, pattern="X", replacement="")
    ret
}

tbls[["csx"]]$id <- extract.id(tbls[["csx"]]$sample.id)

m <- merge(tbls[["gt"]], tbls[["csx"]], by = c("id", "cell.type"))

## pdf("plots.pdf", onefile=TRUE)
cors <- 
    ddply(m, .variables = c("cell.type", "subchallenge"),
          .fun = function(df) {
              cell.type <- df[1, "cell.type"]
              sc <- df[1, "subchallenge"]
              ct.pearson <- cor.test(df$prediction, df$measured, method = "pearson")
              ct.spearman <- cor.test(df$prediction, df$measured, method = "spearman")              
              g <- ggplot(data = df)
              g <- g + geom_point(aes(x = prediction, y = measured))
              g <- g + xlab("Prediction") + ylab("Measured")
              r.str <- round(ct.pearson$estimate, digits=2)
              r.pval <- ct.pearson$p.value
              digits <- 2
              r.pval <- format(r.pval, digits=digits, scientific=0)
              rho.str <- round(ct.spearman$estimate, digits=2)
              rho.pval <- ct.spearman$p.value
              rho.pval <- format(rho.pval, digits=digits, scientific=0)              
              g <- g + ggtitle(paste0(cell.type, " ", sc, ": r = ", r.str, " (p = ", r.pval, ")",
                                      " rho = ", rho.str, " (p = ", rho.pval, ")"))
              file <- paste0("csx-ball-", make.names(cell.type), "-", sc, ".png")
              png(file)
              print(g)
              d <- dev.off()
              data.frame(r = ct.pearson$estimate, p.pearson = ct.pearson$p.value,
                         rho = ct.spearman$estimate, p.spearman = ct.spearman$p.value)
          })
## d <- dev.off()

o <- order(cors$rho, decreasing = TRUE)
cors <- cors[o,] 

cols <- colnames(cors)
cols <- cols[cols != "subchallenge"]
cors <- unique(cors[, cols])

summary.tbl <- cors
estimate.cols <- c("r", "rho")
for(col in estimate.cols) {
    summary.tbl[, col] <- round(summary.tbl[, col], digits=2)
}

pval.cols <- c("p.pearson", "p.spearman")
for(col in pval.cols) {
    summary.tbl[, col] <- format(summary.tbl[, col], digits=2, scientific=0)
}

p_load(gridExtra)
p_load(grid)
png("csx-ball-summary.png")
grid.table(summary.tbl, rows = NULL)
d <- dev.off()
