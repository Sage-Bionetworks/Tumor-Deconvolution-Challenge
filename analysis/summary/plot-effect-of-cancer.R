
suppressPackageStartupMessages(library(pacman))

suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(tidyr))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(ggpubr))
suppressPackageStartupMessages(p_load(grid))
suppressPackageStartupMessages(p_load(ggbeeswarm))

synLogin()
synId <- "syn21739521"
in.silico.results.file <- synGet(synId, downloadFile = TRUE)$path
res <- read.table(in.silico.results.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)

synId <- "syn21752552"
in.silico.coarse.admixture.file <- synGet(synId, downloadFile = TRUE)$path
in.silico.coarse.admixtures <- read.table(in.silico.coarse.admixture.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
in.silico.coarse.admixtures <- subset(in.silico.coarse.admixtures, spike.in.pop %in% c("Breast", "CRC"))
in.silico.coarse.admixtures <- in.silico.coarse.admixtures[, c("dataset.name", "sample.id", "spike.in.pop", "spike.in.prop")]
in.silico.coarse.admixtures <- na.omit(in.silico.coarse.admixtures)

synId <- "syn21752551"
in.silico.fine.admixture.file <- synGet(synId, downloadFile = TRUE)$path
in.silico.fine.admixtures <- read.table(in.silico.fine.admixture.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
in.silico.fine.admixtures <- subset(in.silico.fine.admixtures, spike.in.pop %in% c("Breast", "CRC"))
in.silico.fine.admixtures <- in.silico.fine.admixtures[, c("dataset.name", "sample.id", "spike.in.pop", "spike.in.prop")]
in.silico.fine.admixtures <- na.omit(in.silico.fine.admixtures)

res.coarse <- merge(subset(res, subchallenge == "coarse"), in.silico.coarse.admixtures)
res.fine <- merge(subset(res, subchallenge == "fine"), in.silico.fine.admixtures)

if(FALSE) {
all.gt <- rbind(in.silico.coarse.admixtures, in.silico.fine.admixtures)
all.gt <- subset(all.gt, spike.in.pop %in% c("Breast", "CRC"))
all.gt <- unique(all.gt[, c("dataset.name", "sample.id", "spike.in.pop", "spike.in.prop")])
all.gt <- na.omit(all.gt)
res <- merge(res, all.gt)
}

trans <-
    list("baseline_method1" = "CIBERSORT",
         "baseline_method2" = "MCP-counter",
         "baseline_method3" = "quanTIseq",
         "baseline_method4" = "xCell",
         "baseline_method5" = "EPIC",
         "baseline_method6" = "TIMER",
         "CelEst" = "EPIC (applicant)")
for(nm in names(trans)) {
    flag <- grepl(res$repo_name, pattern = nm)
    res[flag, "repo_name"] <- trans[[nm]]
}

res.coarse <- subset(res, dataset.name %in% c("G", "H"))

sm <-
    ddply(res.coarse,
          .variables = c("round", "subchallenge", "repo_name", "cell.type"),
          .fun = function(df) {
              ret <-
                  ddply(df,
                        .variables = c("dataset.name"),
                        .fun = function(df.ds) {
                            print(df.ds)
                            cr <- cor(df.ds$prediction, df.ds$measured)
                            data.frame(cor = cr)
                        })
              data.frame(cor = mean(ret$cor))
          })

sm <-
    ddply(res.coarse,
          .variables = c("round", "subchallenge", "repo_name", "cell.type",
                         "spike.in.pop", "spike.in.prop"),
          .fun = function(df) {
              ret <-
                  ddply(df,
                        .variables = c("dataset.name"),
                        .fun = function(df.ds) {
                            print(df.ds)
                            cr <- cor(df.ds$prediction, df.ds$measured)
                            data.frame(cor = cr)
                        })
              data.frame(cor = mean(ret$cor))
          })

sm <-
    ddply(res.coarse,
          .variables = c("dataset.name", "round", "subchallenge", "repo_name", "cell.type",
                         "spike.in.pop", "spike.in.prop"),
          .fun = function(df) {
              cr <- cor(df$prediction, df$measured)
              data.frame(cor = cr)
          })

sm <-
    ddply(res.fine,
          .variables = c("round", "subchallenge", "repo_name", "cell.type",
                         "spike.in.pop", "spike.in.prop"),
          .fun = function(df) {
              ret <-
                  ddply(df,
                        .variables = c("dataset.name"),
                        .fun = function(df.ds) {
                            print(df.ds)
                            stop("stop")
                            cr <- cor(df.ds$prediction, df.ds$measured)
                            data.frame(cor = cr)
                        })
              data.frame(cor = mean(ret$cor))
          })

                                        
plts <-
    dlply(sm,
          .variables = c("repo_name"),
          .fun = function(df) {
              g <- ggplot(data = df, aes(x = spike.in.prop, y = cor))
              g <- g + geom_beeswarm()
              g <- g + facet_wrap(~ cell.type)
              g <- g + ggtitle(df$repo_name[1])
          })

cat("Exiting successfully\n")
q(status=0)

