
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
##synId <- "syn21739521"
##in.silico.results.file <- synGet(synId, downloadFile = TRUE)$path
##res <- read.table(in.silico.results.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
in.silico.results.file <- "all-predictions.tsv"
res <- read.table(in.silico.results.file, sep="\t", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
source("../utils.R")
res$repo_name <- res$method

val.metadata <- get.in.silico.metadata()
## res <- merge(res, val.metadata, all.x = TRUE, by = c("dataset.name", "subchallenge"))
res <- merge(res, val.metadata, all.x = FALSE, by = c("dataset.name", "subchallenge"))

## Let's exclude memory.B.cells (which always have measured == 0, which causes problems with correlation)
res <- subset(res, !(cell.type == "memory.B.cells"))

## synId <- "syn21752552"
synId <- "syn21763908"
in.silico.coarse.admixture.file <- synGet(synId, downloadFile = TRUE)$path
in.silico.coarse.admixtures <- read.table(in.silico.coarse.admixture.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
## in.silico.coarse.admixtures <- subset(in.silico.coarse.admixtures, spike.in.pop %in% c("Breast", "CRC"))
in.silico.coarse.admixtures <- in.silico.coarse.admixtures[, c("dataset.name", "sample", "sample.id", "spike.in.pop", "spike.in.prop", "measured")]
colnames(in.silico.coarse.admixtures) <- c("dataset.name", "cell.type", "sample.id", "spike.in.pop", "spike.in.prop", "measured")
in.silico.coarse.admixtures <-
    ddply(in.silico.coarse.admixtures,
          .variables = c("dataset.name", "cell.type", "sample.id", "spike.in.pop", "spike.in.prop"),
          .fun = function(df) {
              df$measured <- sum(df$measured)
              df[1,,drop=F]
          })
tmp <- na.omit(unique(in.silico.coarse.admixtures[, c("dataset.name", "sample.id", "spike.in.pop", "spike.in.prop")]))
tmp <- subset(tmp, spike.in.pop %in% c("CRC", "BRCA", "Breast", "breast"))
in.silico.coarse.admixtures <- merge(in.silico.coarse.admixtures[, !(colnames(in.silico.coarse.admixtures) %in% c("spike.in.pop", "spike.in.prop"))],
                                     tmp)


## synId <- "syn21752551"
synId <- "syn21763907"
in.silico.fine.admixture.file <- synGet(synId, downloadFile = TRUE)$path
in.silico.fine.admixtures <- read.table(in.silico.fine.admixture.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
## in.silico.fine.admixtures <- subset(in.silico.fine.admixtures, spike.in.pop %in% c("Breast", "CRC"))
in.silico.fine.admixtures <- in.silico.fine.admixtures[, c("dataset.name", "sample", "sample.id", "spike.in.pop", "spike.in.prop", "measured")]
colnames(in.silico.fine.admixtures) <- c("dataset.name", "cell.type", "sample.id", "spike.in.pop", "spike.in.prop", "measured")

tmp <- na.omit(unique(in.silico.fine.admixtures[, c("dataset.name", "sample.id", "spike.in.pop", "spike.in.prop")]))
tmp <- subset(tmp, spike.in.pop %in% c("CRC", "BRCA", "Breast", "breast"))
in.silico.fine.admixtures <- merge(in.silico.fine.admixtures[, !(colnames(in.silico.fine.admixtures) %in% c("spike.in.pop", "spike.in.prop"))],
                                     tmp)



res.coarse <- merge(subset(res, subchallenge == "coarse"), in.silico.coarse.admixtures)
res.fine <- merge(subset(res, subchallenge == "fine"), in.silico.fine.admixtures)

all.gt <- rbind(in.silico.coarse.admixtures, in.silico.fine.admixtures)
all.gt <- subset(all.gt, spike.in.pop %in% c("Breast", "CRC"))
all.gt <- unique(all.gt[, c("dataset.name", "cell.type", "sample.id", "spike.in.pop", "spike.in.prop", "measured")])


## all.gt <- na.omit(all.gt)
flag <- !is.na(all.gt$spike.in.pop) & (all.gt$spike.in.pop == "Breast")
all.gt$spike.in.pop[flag] <- "BRCA"
res <- merge(res, all.gt)

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

res.cor <-
    ddply(res,
          .variables = c("dataset.name", "mixture.type", "cell.type", "tumor.type", "spike.in.pop",
                         "spike.in.prop", "repo_name", "subchallenge"),
          .fun = function(df) {
              data.frame(cor.p = cor(df$measured, df$prediction, method = "pearson"),
                         cor.s = cor(df$measured, df$prediction, method = "spearman"))
          })

res.cor$spike.in.percent <- 100 * res.cor$spike.in.prop
res.cor$spike.in.percent <- factor(res.cor$spike.in.percent, levels = sort(unique(res.cor$spike.in.percent), decreasing = FALSE))
## res.cor$spike.in.prop <- factor(res.cor$spike.in.prop, levels = sort(unique(res.cor$spike.in.prop), decreasing = FALSE))

cor.types <- list("Pearson" = "Pearson", "Spearman" = "Spearman")
plts <-
    llply(cor.types,
          .fun = function(cor.type) {
              dlply(res.cor,
                    .variables = c("subchallenge", "repo_name", "mixture.type"),
                    .fun = function(df) {
                        cor.var <- "null"
                        ylab <- "null"
                        if(cor.type == "Pearson") {
                            cor.var <- "cor.p"
                            ylab <- "Pearson Correlation"
                        }
                        if(cor.type == "Spearman") {
                            cor.var <- "cor.s"
                            ylab <- "Spearman Correlation"
                        }
                        g <- ggplot(data = df, aes_string(x = "spike.in.percent", y = cor.var, colour = "spike.in.pop"))
                        g <- g + geom_point()
                        ## g <- g + geom_smooth(data = df, aes(x = spike.in.prop, y = cor.var), method = "lm")
                        g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
                        g <- g + facet_wrap(~ cell.type, scales = "free")
                        meth <- df$repo_name[1]
                        sc <- df$subchallenge[1]
                        mt <- df$mixture.type[1]
                        g <- g + xlab("Cancer Spike-in (%)") + ylab(ylab)
                        g <- g + ggtitle(paste0(meth, " (", firstup(sc), "-Grained Sub-Challenge; ", mt, " Admixtures)"))
                        g <- g + labs(colour = "Cancer Type")
                        g
                    })
          })

all.plts <-
    llply(cor.types,
          .fun = function(cor.type) {
              dlply(res.cor,
                    .variables = c("subchallenge", "repo_name", "mixture.type", "cell.type"),
                    .fun = function(df) {
                        cor.var <- "null"
                        ylab <- "null"
                        if(cor.type == "Pearson") {
                            cor.var <- "cor.p"
                            ylab <- "Pearson Correlation"
                        }
                        if(cor.type == "Spearman") {
                            cor.var <- "cor.s"
                            ylab <- "Spearman Correlation"
                        }
                        g <- ggplot(data = df, aes_string(x = "spike.in.percent", y = cor.var, colour = "spike.in.pop"))
                        g <- g + geom_point()
                        ## g <- g + geom_smooth(data = df, aes(x = spike.in.prop, y = cor.var), method = "lm")
                        g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
                        meth <- df$repo_name[1]
                        sc <- df$subchallenge[1]
                        mt <- df$mixture.type[1]
                        ct <- df$cell.type[1]
                        g <- g + xlab("Cancer Spike-in (%)") + ylab(ylab)
                        g <- g + ggtitle(paste0(meth, " ", ct, "\n(", firstup(sc), "-Grained Sub-Challenge; ", mt, " Admixtures)"))
                        g <- g + labs(colour = "Cancer Type")
                        g
                    })
          })

cancer.cor <-
    ddply(res.cor,
          .variables = c("subchallenge", "repo_name", "mixture.type", "cell.type"),
          .fun = function(df) {
              cor.of.cor.p <- cor(df$cor.p, df$spike.in.prop, method = "pearson")
              cor.of.cor.s <- cor(df$cor.p, df$spike.in.prop, method = "spearman")
              lm.p <- lm(data = df, cor.p ~ spike.in.prop)
              coeffs.p <- coefficients(summary(lm.p))
              lm.p.effect <- coeffs.p["spike.in.prop","Estimate"]
              lm.p.pval <- coeffs.p["spike.in.prop","Pr(>|t|)"]
              lm.s <- lm(data = df, cor.s ~ spike.in.prop)
              coeffs.s <- coefficients(summary(lm.s))
              lm.s.effect <- coeffs.s["spike.in.prop","Estimate"]
              lm.s.pval <- coeffs.s["spike.in.prop","Pr(>|t|)"]
              data.frame(cor.of.cor.p = cor.of.cor.p, cor.of.cor.s = cor.of.cor.s,
                         lm.p.effect = lm.p.effect, lm.p.pval = lm.p.pval,
                         lm.s.effect = lm.s.effect, lm.s.pval = lm.s.pval)
          })


if(FALSE) {
cancer.cor.plots <-
    llply(cor.types,
          .fun = function(cor.type) {
              cor.var <- "null"
              if(cor.type == "Pearson") { cor.var <- "cor.of.cor.p" }
              if(cor.type == "Spearman") { cor.var <- "cor.of.cor.s" }
              
              dlply(cancer.cor,
                    .variables = c("subchallenge", "mixture.type"),
                    .fun = function(df) {
                        g <- plot.cell.type.correlation.heatmap(df, show.corr.text = TRUE,
                                                                id.var = "repo_name", cor.var = cor.var,
                                                                cor.type.label = cor.type, digits = 1)
                        g
                    })
              
          })
}
cancer.cor.plots <-
    llply(cor.types,
          .fun = function(cor.type) {
              cor.var <- "null"
              pval.var <- "null"
              if(cor.type == "Pearson") {
                  cor.var <- "lm.p.effect"
                  pval.var <- "lm.p.pval"
              }
              if(cor.type == "Spearman") {
                  cor.var <- "lm.s.effect"
                  pval.var <- "lm.s.pval"                  
              }
              dlply(cancer.cor,
                    .variables = c("subchallenge", "mixture.type"),
                    .fun = function(df) {
                        label <- paste0("Effect on\n", cor.type, "\nCorrelation")
                        g <- plot.cell.type.correlation.heatmap(df, show.corr.text = TRUE,
                                                                id.var = "repo_name", cor.var = cor.var,
                                                                pval.var = pval.var,
                                                                cor.type.label = label, digits = 1,
                                                                limits = c(-1,1))
                        g
                    })
              
          })

for(cor.type in c("pearson", "spearman")) {
    for(sc in c("coarse", "fine")) {
        g1 <- cancer.cor.plots[[firstup(cor.type)]][[paste0(sc, ".Random")]]
        g1 <- g1 + ggtitle("Random Admixtures")
        g2 <- cancer.cor.plots[[firstup(cor.type)]][[paste0(sc, ".Biological")]]
        g2 <- g2 + ggtitle("Biological Admixtures")
        title <- paste0(firstup(sc), "-Grained Sub-Challenge (Validation)")
        png(paste0("validation-cancer-", sc, "-bio-and-rand-", cor.type, ".png"), width = 2 * 480)
        g <- grid.arrange(g1, g2, ncol = 2, top = textGrob(title, gp = gpar(fontsize = 25)))
        grid.draw(g)
        d <- dev.off()
    }
}

for(cor.type in cor.types) {
    d_ply(res.cor,
          .variables = c("repo_name", "subchallenge", "mixture.type"),
          .fun = function(df) {
              mt <- as.character(df$mixture.type[1])
              meth <- as.character(df$repo_name[1])
              sc <- as.character(df$subchallenge[1])
              cor.var <- "null"
              if(cor.type == "Pearson") {
                  cor.var <- "cor.p"
              }
              if(cor.type == "Spearman") {
                  cor.var <- "cor.s"
              }

              g <- ggplot(df, aes_string(x = "spike.in.percent", y = cor.var))
              g <- g + geom_boxplot()
              g <- g + geom_beeswarm()
              g <- g + facet_wrap(~ cell.type, scales = "free")
              g <- g + xlab("Cancer Spike-in (%)") + ylab(paste0(firstup(cor.type), " Correlation"))
              sz <- 15
              g <- g + theme(text = element_text(size = sz), title = element_text(size = sz))
              g <- g + ggtitle(paste0(meth, ": ", mt, "\n(", firstup(sc), "-Grained Sub-Challenge)"))
              g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))              
              file <- paste0("validation-cancer-method-", meth, "-", sc, "-", mt, "-", cor.type, ".png")
              print(file)
              png(file)
              print(g)
              d <- dev.off()
          })
}

for(cor.type in cor.types) {
    d_ply(res.cor,
          .variables = c("repo_name", "subchallenge"),
          .fun = function(df) {
              mt <- as.character(df$mixture.type[1])
              meth <- as.character(df$repo_name[1])
              sc <- as.character(df$subchallenge[1])
              cor.var <- "null"
              if(cor.type == "Pearson") {
                  cor.var <- "cor.p"
              }
              if(cor.type == "Spearman") {
                  cor.var <- "cor.s"
              }
              sz <- 15

              df.rand <- subset(df, mixture.type == "Random")
              df.bio <- subset(df, mixture.type == "Biological")              

              g1 <- ggplot(df.rand, aes_string(x = "spike.in.percent", y = cor.var))
              g1 <- g1 + geom_boxplot()
              g1 <- g1 + geom_beeswarm()
              g1 <- g1 + facet_wrap(~ cell.type, scales = "free")
              g1 <- g1 + xlab("Cancer Spike-in (%)") + ylab(paste0(firstup(cor.type), " Correlation"))
              g1 <- g1 + theme(text = element_text(size = sz), title = element_text(size = sz))
              g1 <- g1 + ggtitle(paste0(meth, ": ", mt, "\n(", firstup(sc), "-Grained Sub-Challenge)"))
              g1 <- g1 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
              g1 <- g1 + ggtitle("Random Admixtures")

              g2 <- ggplot(df.bio, aes_string(x = "spike.in.percent", y = cor.var))
              g2 <- g2 + geom_boxplot()
              g2 <- g2 + geom_beeswarm()
              g2 <- g2 + facet_wrap(~ cell.type, scales = "free")
              g2 <- g2 + xlab("Cancer Spike-in (%)") + ylab(paste0(firstup(cor.type), " Correlation"))
              g2 <- g2 + theme(text = element_text(size = sz), title = element_text(size = sz))
              g2 <- g2 + ggtitle(paste0(meth, ": ", mt, "\n(", firstup(sc), "-Grained Sub-Challenge)"))
              g2 <- g2 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
              g2 <- g2 + ggtitle("Biological Admixtures")

              file <- paste0("validation-cancer-method-", meth, "-", sc, "-", cor.type, ".png")
              print(file)
              png(file, width = 2 * 480)
              title <- paste0(meth, " (", firstup(sc), "-Grained Sub-Challenge)")
              g <- grid.arrange(g1, g2, ncol = 2, top = textGrob(title, gp=gpar(fontsize=20)))
              ## print(g)
              d <- dev.off()
          })
}

d_ply(res.cor,
      .variables = c("cell.type", "repo_name", "subchallenge"),
      .fun = function(df) {
          ct <- as.character(df$cell.type[1])
          meth <- as.character(df$repo_name[1])
          sc <- as.character(df$subchallenge[1])
          cor.var <- NULL
          df.melt <- data.frame(mixture.type = paste0(df$mixture.type, " Admixture"),
                                correlation = df$cor.p,
                                correlation.type = "Pearson",
                                spike.in.percent = df$spike.in.percent)
          df.melt <- rbind(df.melt,
                           data.frame(mixture.type = paste0(df$mixture.type, " Admixture"),
                                      correlation = df$cor.s,
                                      correlation.type = "Spearman",
                                      spike.in.percent = df$spike.in.percent))
          g <- ggplot(df.melt, aes(x = spike.in.percent, y = correlation))
          ## g <- g + geom_point()
          g <- g + geom_boxplot()
          ## g <- g + geom_violin()
          g <- g + geom_beeswarm()
          g <- g + facet_grid(correlation.type ~ mixture.type)
          g <- g + xlab("Cancer Spike-in (%)") + ylab("Correlation")
          sz <- 25
          g <- g + theme(text = element_text(size = sz), title = element_text(size = sz))
          g <- g + ggtitle(paste0(meth, ": ", ct, " (", firstup(sc), "-Grained Sub-Challenge)"))
          file <- paste0("validation-cancer-cell-type-", ct, "-", meth, "-", sc, ".png")
          print(file)
          png(file, width = 2 * 480, height = 2 * 480)
          print(g)
          d <- dev.off()
      })

sub <- subset(res, cell.type == "B.cells" & subchallenge == "coarse" & repo_name == "MCP-counter" &
                   spike.in.prop == 0 & mixture.type == "Biological")

g <- ggplot(sub, aes(x = measured, y = prediction))
g <- g + geom_point()
g <- g + facet_wrap(~ dataset.name, scales = "free")
g <- g + geom_smooth(method = "lm")
g <- g + ggtitle(paste0("MCP-counter B-cells (0% Cancer; Biological)"))
g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
png("validation-cancer-MCP-counter-B-cell-no-cancer-biological.png")
print(g)
d <- dev.off()

sub <- subset(res, cell.type == "B.cells" & subchallenge == "coarse" & repo_name == "MCP-counter" &
                   spike.in.prop == 0 & mixture.type == "Random")

g <- ggplot(sub, aes(x = measured, y = prediction))
g <- g + geom_point()
g <- g + facet_wrap(~ dataset.name, scales = "free")
g <- g + geom_smooth(method = "lm")
g <- g + ggtitle(paste0("MCP-counter B-cells (0% Cancer; Random)"))
g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
png("validation-cancer-MCP-counter-B-cell-no-cancer-random.png")
print(g)
d <- dev.off()

cat("Exiting successfully\n")
stop("stop")
q(status=0)

for(meth in unique(res.cor$repo_name)) {
    for(cor.type in cor.types) {
        for(sc in unique(res.cor$subchallenge)) {
            g1 <- plts[[cor.type]][[paste0(sc, ".", meth, ".", "Random")]]
            g1 <- g1 + ggtitle("Random Admixtures")
            sz <- 15
            g1 <- g1 + theme(text = element_text(size = sz), title = element_text(size = sz))
            g2 <- plts[[cor.type]][[paste0(sc, ".", meth, ".", "Biological")]]
            g2 <- g2 + ggtitle("Biological Admixtures")
            g2 <- g2 + theme(text = element_text(size = sz), title = element_text(size = sz))
            title <- paste0(meth, " (", firstup(sc), "-Grained Sub-Challenge)")
            g <- grid.arrange(g1, g2, ncol = 2, top = textGrob(title, gp=gpar(fontsize=20)))
            file <- paste0("validation-cancer-", meth, "-", sc, "-", cor.type, ".png")
            png(file, width = 2 * 480, height = 2 * 480)
            grid.draw(g)
            d <- dev.off()
        }
    }
}



## res.coarse <- subset(res, dataset.name %in% c("G", "H"))

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

