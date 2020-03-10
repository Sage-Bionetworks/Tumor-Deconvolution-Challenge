suppressPackageStartupMessages(library(pacman))

suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(tidyr))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(ComplexHeatmap))
suppressPackageStartupMessages(p_load(reshape2))

synLogin()

file <- "all-leaderboard-results.tsv"
synId <- "syn21714436"
file <- synGet(synId, downloadFile = TRUE)$path
res <- read.table(file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)

## Limit to the final (of two) submissions for each group
## and the baseline methods (which were all submitted by Andrew L
## and hence only one of which is_latest)
res <- subset(res, ( is_latest == TRUE ) | ( grepl(repo_name, pattern = "baseline")))
res <- subset(res, !is.na(measured))

res$modelId <- paste0(as.character(res$repo_name), "-", as.character(res$submitterId))

sub.challenges <- list("coarse" = "coarse", "fine" = "fine")
rounds <- unique(res$round)
names(rounds) <- rounds

n.bootstraps <- 100
bootstraps <- 1:n.bootstraps
names(bootstraps) <- bootstraps

set.seed(1234)

bootstrapped.cors <-
    ddply(res,
          .variables = c("subchallenge"),
          .fun = function(df.sc) {
              ddply(df.sc,
                    .variables = c("round"),
                    .fun = function(df.rd) {
                        
                        tmp <-
                            ldply(bootstraps,
                              .fun = function(i) {
                                  ret <-
                                      ddply(df.rd, .variables = c("dataset.name"),
                                            .fun = function(df.ds) {
                                                sample.ids <- unique(df.ds$sample.id)
                                                ## NB: use the same bootstrap samples across datasets and celltypes (nested below)
                                                boot.samples <- sample(sample.ids, size = length(sample.ids), replace = TRUE)
                                                eps.vec <- runif(length(boot.samples), min=10^-10, max=10^-9)
                                                ret.ds <-
                                                    ddply(df.ds, .variables = c("modelId"),
                                                          .fun = function(df.method) {
                                                              ret.ct <-
                                                                  ddply(df.method, .variables = c("cell.type"),
                                                                        .fun = function(df.ct) {
                                                                            rownames(df.ct) <- df.ct$sample.id
                                                                            ## This condition is necessary because a few datasets lack
                                                                            ## _measurements_ (and hence predictions) for some cell types in
                                                                            ## some (but not all) samples
                                                                            ev <- eps.vec[boot.samples %in% rownames(df.ct)]
                                                                            boot.samples <- boot.samples[boot.samples %in% rownames(df.ct)]
                                                                            if(length(boot.samples) < 2) { return(NULL) }
                                                                            boot.df <- df.ct[boot.samples, c("prediction", "measured")]
                                                                            val <- 0
                                                                            ## In some cases all measured values are zero, if so just add
                                                                            ## a little noise
                                                                            if(!all(boot.df$measured == boot.df$measured[1])) { ev <- 0 }
                                                                            if(var(boot.df$prediction) != 0) {
                                                                                val <- cor(boot.df$prediction, boot.df$measured + ev)
                                                                                if(is.na(val)) {
                                                                                    print(boot.df)
                                                                                    stop("stop")
                                                                                }
                                                                            }
                                                                            data.frame(cor = val)
                                                                        })
                                                          })
                                            })
                              })
                        colnames(tmp)[1] <- "bootstrap"
                        tmp
                    })
          })

write.table(file = "bootstraps.tsv", bootstrapped.cors, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

## bootstrapped.cors <- read.table("bootstraps.tsv", sep = "\t", header = TRUE)

bootstrapped.scores <-
    ddply(bootstrapped.cors,
          .variables = c("subchallenge", "round", "bootstrap", "modelId"),
          .fun = function(df) {
              res.ds <-
                  ddply(df, .variables = c("dataset.name"),
                        .fun = function(df.ds) {
                            data.frame(ds.score = mean(df$cor))
                        })
              data.frame(score = mean(res.ds$ds.score))
          })

bootstrapped.scores.sv <- bootstrapped.scores

model.name.tbl <- merge(unique(bootstrapped.scores[, c("modelId"), drop = FALSE]),
                        unique(res[, c("modelId", "submitterId")]))
model.name.tbl <- unique(model.name.tbl[, c("submitterId"), drop = FALSE])
my.letters <- c(LETTERS, paste0("A", LETTERS))
model.name.tbl$id <- my.letters[1:nrow(model.name.tbl)]
model.name.tbl <- merge(model.name.tbl, unique(res[, c("submitterId", "modelId")]))
colnames(model.name.tbl) <- c("submitterId", "id", "modelId")
model.name.tbl <- merge(model.name.tbl, unique(res[, c("repo_name", "modelId")]))
flag <- grepl(model.name.tbl$repo_name, pattern = "baseline_method")
## model.name.tbl[flag, "id"] <- gsub(model.name.tbl[flag, "repo_name"], pattern = "method", replacement = "")

## Translate the baseline models
translate.baseline.models <- TRUE
if(!translate.baseline.models) {
    model.name.tbl[flag, "id"] <- gsub(model.name.tbl[flag, "repo_name"], pattern = "method", replacement = "")
} else {
    model.name.tbl[flag, "id"] <- model.name.tbl[flag, "repo_name"]
    trans <-
        list("baseline_method1" = "CIBERSORT",
             "baseline_method2" = "MCP-counter",
             "baseline_method3" = "quanTIseq",
             "baseline_method4" = "xCell",
             "baseline_method5" = "EPIC (Sage)",
             "baseline_method6" = "TIMER",
             "CelEst" = "EPIC")
    for(nm in names(trans)) {
        flag <- grepl(model.name.tbl$id, pattern = nm)
        model.name.tbl[flag, "id"] <- trans[[nm]]
    }
}



bootstrapped.scores <- merge(bootstrapped.scores, unique(model.name.tbl[, c("modelId", "id")]))

facet.labels <- list("1" = "Round 1", "2" = "Round 2", "3" = "Round 3")
my.facet.labeller <- function(value) {
    paste0("Round ", value)
}

plot.bootstraps <- function(df) {
    means <- ddply(df, .variables = c("round", "modelId"),
                   .fun = function(tmp) data.frame(mean = mean(tmp$score, na.rm=TRUE)))
    means <- merge(means, unique(model.name.tbl[, c("modelId", "id")]))
    means <- means[order(means$round, means$mean),]
    means <- means[!duplicated(means$id, fromLast = TRUE),]
    df$id <- factor(df$id, levels = unique(means$id))
    g1 <- ggplot(data = df)
    g1 <- g1 + geom_boxplot(aes(x = id, y = score))
    g1 <- g1 + facet_wrap(~ round, labeller = labeller(round = my.facet.labeller))
    int.data <- ddply(means, .variables = c("round"), .fun = function(df.rd) data.frame(mean = mean(df.rd$mean)))
    g1 <- g1 + geom_hline(data = int.data, aes(yintercept = mean), linetype = "dashed")
    g1 <- g1 + coord_flip()
    g1 <- g1 + ylab("Pearson Correlation\nAvg over Dataset & Cell Type")
    g1 <- g1 + theme(text = element_text(size=20), axis.text.y = element_text(size = 20),
                     axis.text.x = element_text(angle = 45, hjust = 1),
                     panel.spacing = unit(1, "lines"))
    g1 <- g1 + xlab("Method")
    g1
}

df <- subset(bootstrapped.scores, subchallenge == "coarse")
g1 <- plot.bootstraps(df)
g1 <- g1 + ggtitle("Coarse-Grained Sub-Challenge (Leaderboard Rounds)")
## pdf("leaderboard-coarse-bootstrap.pdf", width = 14)
png("leaderboard-coarse-bootstrap.png", width = 480 * 2)
print(g1)
d <- dev.off()

df <- subset(bootstrapped.scores, subchallenge == "fine")
g2 <- plot.bootstraps(df)
g2 <- g2 + ggtitle("Fine-Grained Sub-Challenge (Leaderboard Rounds)")
## pdf("leaderboard-fine-bootstrap.pdf", width = 14)
png("leaderboard-fine-bootstrap.png", width = 480 * 2)
print(g2)
d <- dev.off()

mean.bootstrapped.scores <-
    ddply(bootstrapped.scores,
          .variables = c("subchallenge", "round", "modelId", "id"),
          .fun = function(df) {
              data.frame(mean.score = mean(df$score))
          })

bootstrapped.cors <- merge(bootstrapped.cors, unique(model.name.tbl[, c("modelId", "id")]))
means.by.cell.type.method <-
    ddply(bootstrapped.cors,
          .variables = c("subchallenge", "round", "modelId", "cell.type", "id"),
          .fun = function(df) {
              ## first, average over bootstrap
              ret <- ddply(df, .variables = c("dataset.name"),
                           .fun = function(df) {
                               data.frame(cor = mean(df$cor, na.rm=TRUE))
                           })
              ## now, average over dataset
              data.frame(cor = mean(ret$cor))
          })

plot.cell.type.correlation.heatmap <- function(df, show.corr.text = FALSE) {

    df <- df[, c("id", "cell.type", "cor")]
    df$id <- as.character(df$id)
    df$cell.type <- as.character(df$cell.type)
    
    cell.type.means <- ddply(df, .variables = c("cell.type"),
                             .fun = function(tmp) data.frame(id = "mean", cor = mean(tmp$cor, na.rm=TRUE)))
    cell.type.means <- cell.type.means[order(cell.type.means$cor),]
    
    method.means <- ddply(df, .variables = c("id"),
                   .fun = function(tmp) data.frame(cell.type = "mean", cor = mean(tmp$cor, na.rm=TRUE)))
    method.means <- method.means[order(method.means$cor),]


    cell.type.levels <- c(cell.type.means$cell.type[cell.type.means$cell.type != "mean"], "mean")
    id.levels <- c("mean", method.means$id[method.means$id != "mean"])

    df <- rbind(df, cell.type.means[, c("id", "cell.type", "cor")])
    df <- rbind(df, method.means[, c("id", "cell.type", "cor")])    
    df$cell.type <- factor(df$cell.type, levels = cell.type.levels)
    df$id <- factor(df$id, levels = id.levels)

    df$cor.label <- formatC(df$cor, format="f", digits=2)
    g <- ggplot(data = df, aes(y = id, x = cell.type, fill = cor))
    g <- g + geom_tile()
    if(show.corr.text) {
        g <- g + geom_text(aes(label = cor.label))
    }
    g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   axis.text.y = element_text(size = 5), title = element_text(size = 8))
    g <- g + ylab("Method") + xlab("")
    ## g <- g + scale_fill_continuous("Pearson\ncorrelation", limits = c(-1,1))
    ## g <- g + scale_fill_gradient2("Pearson\ncorrelation", limits = c(-1,1),
    ##                               low = "red", high = "blue", mid = "white", na.value = "black")
    g <- g + scale_fill_gradient2("Correlation", limits = c(-1,1), low = "red", high = "blue", mid = "white", na.value = "black")
    ## g <- g + theme(text = element_text(size=20))
    g
}

plot.cell.type.correlation.complex.heatmap <- function(df, show.corr.text = FALSE) {
    means <- ddply(df, .variables = c("cell.type"),
                   .fun = function(tmp) mean(tmp$cor, na.rm=TRUE))
    means <- means[order(means$V1),]
    cell.type.levels <- means$cell.type
    df$cell.type <- factor(df$cell.type, levels = cell.type.levels)
    
    means <- ddply(df, .variables = c("id"),
                   .fun = function(tmp) mean(tmp$cor, na.rm=TRUE))
    means <- means[order(means$V1),]
    id.levels <- means$id
    df$id <- factor(df$id, levels = id.levels)
    

    mat <- as.matrix(acast(df[, c("cell.type", "id", "cor")], id ~ cell.type))
    mat <- mat[id.levels, cell.type.levels]
    Heatmap(mat)

}

glist <-
    dlply(subset(means.by.cell.type.method, subchallenge == "coarse"),
          .variables = c("round"),
          .fun = function(df) {
              g <- plot.cell.type.correlation.heatmap(df)
              ## g <- g + ggtitle(paste0("Coarse-Grained Sub-Challenge\nRound ", df$round[1]))
              g <- g + ggtitle(paste0("Round ", df$round[1]))
              g
          })

png("leaderboard-coarse-cell-type.png")
do.call("grid.arrange", c(glist, nrow = 2, top = "Coarse-Grained Sub-Challenge (Leaderboard Rounds)"))
d <- dev.off()

d_ply(subset(means.by.cell.type.method, subchallenge == "coarse"),
      .variables = c("round"),
      .fun = function(df) {
          g <- plot.cell.type.correlation.heatmap(df, show.corr.text = TRUE)
          rnd <- df$round[1]
          g <- g + ggtitle(paste0("Coarse-Grained Sub-Challenge\nRound ", rnd))
          g <- g + theme(text = element_text(size = 15), title = element_text(size = 15),
                         axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
          png(paste0("leaderboard-coarse-cell-type-round-", rnd, ".png"))
          print(g)
          d <- dev.off()
      })

glist <-
    dlply(subset(means.by.cell.type.method, subchallenge == "fine"),
          .variables = c("round"),
          .fun = function(df) {
              g <- plot.cell.type.correlation.heatmap(df)
              ## g <- g + ggtitle(paste0("Fine-Grained Sub-Challenge\nRound ", df$round[1]))
              g <- g + ggtitle(paste0("Round ", df$round[1]))
              g
          })

png("leaderboard-fine-cell-type.png")
do.call("grid.arrange", c(glist, nrow = 2, top = "Fine-Grained Sub-Challenge (Leaderboard Rounds)"))
d <- dev.off()

d_ply(subset(means.by.cell.type.method, subchallenge == "fine"),
      .variables = c("round"),
      .fun = function(df) {
          g <- plot.cell.type.correlation.heatmap(df, show.corr.text = TRUE)
          rnd <- df$round[1]
          g <- g + ggtitle(paste0("Fine-Grained Sub-Challenge\nRound ", rnd))
          g <- g + theme(text = element_text(size = 15), title = element_text(size = 15),
                         axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
          png(paste0("leaderboard-fine-cell-type-round-", rnd, ".png"), width = 2 * 480)
          print(g)
          d <- dev.off()
      })

## Plot the rank of the top n (e.g., n = 5) across the 3 rounds
plot.ranks <- function(df, n.top = 5) {
    df$round <- as.character(df$round)
    tmp <-
        ddply(df,
              .variables = c("subchallenge", "round"),
              .fun = function(sub.df) {
                  sub.df$rank <- rank(- sub.df$mean.score)
                  sub.df
              })
    top.performers <- as.character(subset(tmp, rank <= n.top)$id)
    g <- ggplot(subset(tmp, id %in% top.performers), aes(x = round, y = rank, colour = id, group = id))
    g <- g + geom_line()
    g <- g + xlab("Round") + ylab("Rank") + scale_color_discrete(name = "Method")
    text.df <- subset(tmp, (id %in% top.performers) & (round == 3))
    g <- g + geom_text(data = text.df, aes(x = round, y = rank, label = id), hjust = -0.1)
    g <- g + theme(legend.position = "none")
    g
}

png("leaderboard-coarse-rank-rounds.png")
g <- plot.ranks(subset(mean.bootstrapped.scores, subchallenge == "coarse"), n.top = 5)
g <- g + ggtitle("Coarse-Grained Sub-Challenge (Leaderboard Rounds)")
print(g)
d <- dev.off()

png("leaderboard-fine-rank-rounds.png")
g <- plot.ranks(subset(mean.bootstrapped.scores, subchallenge == "fine"), n.top = 5)
g <- g + ggtitle("Fine-Grained Sub-Challenge (Leaderboard Rounds)")
print(g)
d <- dev.off()
