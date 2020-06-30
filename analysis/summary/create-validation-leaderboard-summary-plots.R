suppressPackageStartupMessages(library(pacman))

suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(tidyr))
suppressPackageStartupMessages(p_load(grid))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(synapser))

suppressPackageStartupMessages(p_load("openxlsx"))
suppressPackageStartupMessages(p_load("reshape2"))


## Get the validation results ("all_predictions.csv")
synLogin()

file <- "validation-csx-all-gene-predictions.tsv"
csx.res <- read.table(file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
source("../utils.R")

synId <- "syn21590364"
path <- synGet(synId, downloadFile = TRUE)$path
val.coarse <- read.table(path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
val.coarse$subchallenge <- "coarse"

synId <- "syn21590365"
path <- synGet(synId, downloadFile = TRUE)$path
val.fine <- read.table(path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
val.fine$subchallenge <- "fine"

val.all <- rbind(val.coarse, val.fine)
flag <- val.all$sample.id == "Breast"
val.all[flag, "sample.id"] <- "BRCA"

flag <- csx.res$sample.id == "Breast"
csx.res[flag, "sample.id"] <- "BRCA"

csx.res <- merge(csx.res[, !(colnames(csx.res) == "measured")], val.all)

csx.res <- csx.res[, c("method", "dataset.name", "cell.type", "subchallenge", "measured", "prediction")]
colnames(csx.res) <- c("repo_name", "dataset", "celltype", "subchallenge", "measured", "prediction")

metrics <- list("pearson" = "pearson", "spearman" = "spearman")
cors <-
    ldply(metrics,
          .fun = function(metric) {
              ddply(csx.res,
                    .variables = c("repo_name", "dataset", "celltype", "subchallenge"),
                    .fun = function(df) {
                        val <- cor(df$measured, df$prediction, method = metric)
                        data.frame(metric_value = val)
                    })
          })
colnames(cors) <- c("metric", "repo_name", "dataset", "celltype", "subchallenge", "metric_value")
cors <- cors[, c("repo_name", "dataset", "celltype", "metric", "metric_value", "subchallenge")]

fine.grained.synId <- "syn21822682"
coarse.grained.synId <- "syn21822681"
synIds <- list("fine" = fine.grained.synId, "coarse" = coarse.grained.synId)

nms <- names(synIds)
val.res <- llply(nms,
                 .fun = function(nm) {
                     synId <- synIds[[nm]]
                     df <- synTableQuery(paste0("SELECT * FROM ", synId))$asDataFrame()
                     df$subchallenge <- nm
                     df
                 })

val.res <- do.call("rbind", val.res)
val.res <- val.res[, c("repo_name", "dataset", "celltype", "metric", "metric_value", "subchallenge")]

val.res <- rbind(val.res, cors)

## Exclude the means
val.res <- subset(val.res, !(dataset %in% c("Celltype mean", "Grand mean")))

## Exclude neutrophils
val.res <- subset(val.res, !(celltype == "neutrophils"))

## Translate the baseline methods
trans <-
    list("baseline_method1" = "CIBERSORT",
         "baseline_method2" = "MCP-counter",
         "baseline_method3" = "quanTIseq",
         "baseline_method4" = "xCell",
         "baseline_method5" = "EPIC (Sage)",
         "baseline_method6" = "TIMER",
         "CelEst" = "EPIC")
for(nm in names(trans)) {
    flag <- grepl(val.res$repo_name, pattern = nm)
    val.res[flag, "repo_name"] <- trans[[nm]]
}

val.metadata <- get.validation.metadata()
val.metadata <- unique(val.metadata[, c("tumor.type", "batch", "mixture.type", "dataset")])

val.res <- merge(val.res, val.metadata)

means.by.cell.type.method.mixture.type <-
    ddply(val.res,
          .variables = c("subchallenge", "repo_name", "celltype", "mixture.type", "metric"),
          .fun = function(df) {
              ## average over dataset
              data.frame(metric_value = mean(df$metric_value))
          })


cor.types <- list("Pearson" = "Pearson", "Spearman" = "Spearman")

mt.plts <-
    dlply(means.by.cell.type.method.mixture.type,
          .variables = c("subchallenge", "mixture.type", "metric"),
          .fun = function(df) {
              cor.type <- unique(df$metric)
              g <- plot.cell.type.correlation.heatmap(df, show.corr.text = TRUE, id.var = "repo_name",
                                                      cell.type.var = "celltype",
                                                      cor.var = "metric_value",
                                                      cor.type.label = paste0(firstup(cor.type), "\nCorrelation"))
              sc <- as.character(df$subchallenge[1])
              mt <- as.character(df$mixture.type[1])
              g <- g + ggtitle(paste0(firstup(sc), "-Grained Sub-Challenge\n(Validation; ", mt, ")"))
              g
          })


for(cor.type in c("pearson", "spearman")) {
    for(sc in c("coarse", "fine")) {
        g1 <- mt.plts[[paste0(sc, ".Random.", cor.type)]]
        g1 <- g1 + ggtitle("Random Admixtures")
        g2 <- mt.plts[[paste0(sc, ".Biological.", cor.type)]]        
        g2 <- g2 + ggtitle("Biological Admixtures")
        title <- paste0(firstup(sc), "-Grained Sub-Challenge (Validation", ")")
        png(paste0("validation-leaderboard-synapse-results-", sc, "-bio-and-rand-", cor.type, ".png"), width = 2 * 480)
        g <- grid.arrange(g1, g2, ncol = 2, top = textGrob(title, gp = gpar(fontsize = 25)))
        grid.draw(g)
        d <- dev.off()
    }
}
