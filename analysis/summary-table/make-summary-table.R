suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(ggplot2))
# suppressPackageStartupMessages(p_load(tidyverse))
suppressPackageStartupMessages(p_load(reshape2))
suppressPackageStartupMessages(p_load(cowplot))
# devtools::install_github("r-lib/xml2")
suppressPackageStartupMessages(p_load(xml2))
suppressPackageStartupMessages(p_load(ComplexHeatmap))
suppressPackageStartupMessages(p_load(xlsx))

set.seed(1234)

source("../utils.R")

synLogin()

method.name.col <- "method.name"
cell.type.col <- "cell.type"
nbins <- 4
na.bin <- nbins + 1

figs.dir <- "figs/"
dir.create(figs.dir, showWarnings = FALSE)

## Read in the healthy and cancer predictions (created by score-cancer-datasets.R)
all.res <- read.table("../cancer-validation/figs/cancer-validation-dataset-all-scores.tsv", sep="\t", header=TRUE)

methods.to.exclude <- c("TIMER", "timer")
all.res <- all.res[!(all.res[, method.name.col] %in% methods.to.exclude),]
all.res <- correct.ictd.and.cancer.deconv(all.res, col = method.name.col)

fill.in.nas <- function(mat, row.name, col.name, value.var) {
  ret <- melt(acast(mat, as.formula(paste0(row.name, "~", col.name)), value.var = value.var))
  colnames(ret) <- c(row.name, col.name, value.var)
  ret
}

# from https://stackoverflow.com/questions/75489458/getting-cutoff-values-for-each-ntile-group
# The best bin value is 1.
bin.value <- function(mat, col.name, higher.is.better = FALSE, nbins = 4) {
  # mat$bin_number <- ntile(factor * mat[, col.name], nbins)
  #na.flag <- is.na(mat[, col.name])
  #mat$bin_number <- NA
  #mat[!na.flag, "bin_number"] <- ntile(rank(factor * mat[!na.flag, col.name], ties.method="min"), nbins)
  qs <- quantile(mat[, col.name], probs = seq(0, 1, length.out = nbins + 1), na.rm=TRUE)
  mat$cut_bin <- cut(mat[, col.name], breaks = qs[!duplicated(qs)], include.lowest = TRUE)
  mat$bin_number <- cut(mat[, col.name], breaks = qs[!duplicated(qs)], include.lowest = TRUE, labels=FALSE)
  if(higher.is.better) {
    mat$bin_number <- max(mat$bin_number, na.rm=TRUE) - mat$bin_number + 1
  }
  return(mat)
}

all.res <-
  ddply(all.res, .variables = c("dataset.name"),
        .fun = function(mat) fill.in.nas(mat, row.name = cell.type.col, col.name = method.name.col, value.var = "cor.p"))

# combine all the challenge validation datasets into a 'Healthy' dataset by taking the mean
flag <- all.res$dataset.name %in% c("Pelka", "Wu")
all.res$dataset.name <- as.character(all.res$dataset.name)
all.res[!flag,"dataset.name"] <- "Healthy"

#flag <- all.res$dataset.name %in% c("Wu")
#all.res[flag,"dataset.name"] <- "Wu (BRCA)"
#flag <- all.res$dataset.name %in% c("Pelka")
#all.res[flag,"dataset.name"] <- "Pelka (CRC)"
all.res <- 
  ddply(all.res, .variables = c("dataset.name", "method.name", "cell.type"),
        .fun = function(df) {
          data.frame(cor.p = mean(df$cor.p))
        })

# Now take the average and mean over the 3 datasets (healthy, wu, and pelka)
# Note that some datasets do not have 
mean.all.res <- 
  ddply(all.res, .variables = c("method.name", "cell.type"),
        .fun = function(df) {
          data.frame(mean.cor.p = mean(df$cor.p))
        })

mean.all.res <- 
  ddply(mean.all.res,
        .variables = c("cell.type"),
        .fun = function(df) bin.value(df, "mean.cor.p", higher.is.better = TRUE, nbins = nbins))
mean.all.res[is.na(mean.all.res[, "mean.cor.p"]), "bin_number"] <- na.bin

# Cell types that only occur in one dataset won't have an std dev. Remove those
n.datasets.by.cell.type <-
  ddply(unique(all.res[, c("cell.type", "dataset.name")]), .variables = c("cell.type"), 
        .fun = function(df) data.frame(n.datasets = nrow(df)))
cell.types.to.exclude <- as.character(subset(n.datasets.by.cell.type, n.datasets==1)$cell.type)

sd.all.res <- 
  ddply(subset(all.res, !(cell.type %in% cell.types.to.exclude)), .variables = c("method.name", "cell.type"),
        .fun = function(df) {
          data.frame(sd.cor.p = sd(df$cor.p))
        })

sd.all.res <- 
  ddply(sd.all.res,
        .variables = c("cell.type"),
        .fun = function(df) {
            # print(df[1,"cell.type"])
            bin.value(df, "sd.cor.p", higher.is.better = FALSE, nbins = nbins)
        })

sd.all.res[is.na(sd.all.res[, "sd.cor.p"]), "bin_number"] <- na.bin

# Process limit of detection / spike in
# Note: used random for spike-ins in manuscript, separated by coarse- and fine-grained.
# Need to combine. 
# See analysis/in-silico-admixtures/plot-sensitivity-results.R
synIds <- list("sensitivity-analysis-results" = "syn22953003")
tbls <-
  llply(synIds,
        .fun = function(synId) {
          file <- synGet(synId, downloadFile = TRUE)$path
          readRDS(file)
        })

coarse.spikein <- subset(tbls[["sensitivity-analysis-results"]][["1"]]$res.matrices$coarse, mixture.type == "Random")
fine.spikein <- subset(tbls[["sensitivity-analysis-results"]][["1"]]$res.matrices$fine, mixture.type == "Random")
coarse.spikein$lod <- coarse.spikein$min.diff.prop * 100
coarse.spikein$subchallenge <- "coarse"
fine.spikein$lod <- fine.spikein$min.diff.prop * 100
fine.spikein$subchallenge <- "fine"
stopifnot(length(unique(coarse.spikein$dataset.name)) == 1)
stopifnot(length(unique(fine.spikein$dataset.name)) == 1)

coarse.spikein <- fill.in.nas(coarse.spikein, row.name = cell.type.col, col.name = method.name.col, value.var = "lod")
fine.spikein <- fill.in.nas(fine.spikein, row.name = cell.type.col, col.name = method.name.col, value.var = "lod")

all.spikein <- rbind(coarse.spikein, fine.spikein)
all.spikein <- subset(all.spikein, !(method.name %in% c("TIMER", "timer")))
# Merge coarse and fine-grained results by taking the mean
all.spikein <-
  ddply(all.spikein,
        .variables = c("cell.type", "method.name"),
        .fun = function(df) data.frame(lod = mean(df$lod)))

all.spikein <- fill.in.nas(all.spikein, row.name = cell.type.col, col.name = method.name.col, value.var = "lod")

all.spikein <- 
  ddply(all.spikein,
        .variables = c("cell.type"),
        .fun = function(df) {
          # print(df[1,"cell.type"])
          bin.value(df, "lod", higher.is.better = FALSE, nbins = nbins)
        })

all.spikein[is.na(all.spikein[, "lod"]), "bin_number"] <- na.bin

# Process  spillover
# specificity-analysis/figs/spillover-summary-round-1.tsv
spillover.res <- read.table("../specificity-analysis/figs/spillover-summary-round-1.tsv", sep="\t", header=TRUE)
spillover.res <- fill.in.nas(spillover.res, row.name = cell.type.col, col.name = method.name.col, value.var = "spillover")
spillover.res <- subset(spillover.res, !(method.name %in% c("TIMER", "timer")))
spillover.res <- 
  ddply(spillover.res,
        .variables = c("cell.type"),
        .fun = function(df) {
          # print(df[1,"cell.type"])
          bin.value(df, "spillover", higher.is.better = FALSE, nbins = nbins)
        })

spillover.res[is.na(spillover.res[, "spillover"]), "bin_number"] <- na.bin

# Process sample-level analyses
sample.level.coarse.res <- read.table("../sample-level-analysis/figs/sample-level-round-1-metric.sum.res-coarse.tsv", sep="\t", header=TRUE)
sample.level.fine.res <- read.table("../sample-level-analysis/figs/sample-level-round-1-metric.sum.res-fine.tsv", sep="\t", header=TRUE)

method.anno <- get.method.annotations()

subchallenge.col <- "subchallenge"
round.col <- "submission"
method.anno.round <- get.round.specific.annotations(method.anno, round = "1")
# Exclude score-based output
method.anno.round <- subset(method.anno.round, Output != "Score")
sample.level.coarse.res <- 
  subset(sample.level.coarse.res, method.name %in% 
           subset(method.anno.round, is.na(subchallenge) | (subchallenge == "coarse"))$method.name)
sample.level.fine.res <- 
  subset(sample.level.fine.res, method.name %in% 
           subset(method.anno.round, is.na(subchallenge) | (subchallenge == "coarse"))$method.name)

sample.level.coarse.bin <- bin.value(sample.level.coarse.res, "Pearson", higher.is.better = TRUE, nbins = nbins)
sample.level.fine.bin <- bin.value(sample.level.fine.res, "Pearson", higher.is.better = TRUE, nbins = nbins)


# Process run times
timing.tbl <- read.table("../timing/method-run-times.tsv", sep="\t", header=TRUE)
timing.tbl <- subset(timing.tbl, !(method.name %in% c("TIMER", "timer")))
timing.coarse.raw <- subset(timing.tbl, (submission == 1) & (subchallenge == "coarse"))
timing.fine.raw <- subset(timing.tbl, (submission == 1) & (subchallenge == "fine"))

timing.coarse.bin <- bin.value(timing.coarse.raw, "exec.time.seconds", higher.is.better = FALSE, nbins = nbins)
timing.fine.bin <- bin.value(timing.fine.raw, "exec.time.seconds", higher.is.better = FALSE, nbins = nbins)


# Assemble all metrics
raw.bins <- 
  list("mean" = mean.all.res,
       "sd" = sd.all.res,
       "lod" = all.spikein,
       "spillover" = spillover.res,
       "sample_coarse" = sample.level.coarse.bin,
       "sample_fine" = sample.level.fine.bin,
       "timing_coarse" = timing.coarse.bin,
       "timing_fine" = timing.fine.bin)

for(nm in names(raw.bins)) {
  write.table(file = paste0(figs.dir, "/", "binned-score-", nm, ".tsv"), raw.bins[[nm]], sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}

all.bin.cutoffs <- 
  ldply(raw.bins, .fun = function(df) {
    cols <- c("cell.type", "cut_bin", "bin_number")
    ret <- df[, cols[cols %in% colnames(df)]]
    if(!("cell.type" %in% colnames(df))) {
      ret <- cbind(cell.type = "NA", ret)
    }
    ret <- ret[, cols]
    unique(ret)
  })
colnames(all.bin.cutoffs)[1] <- "evaluation"
#all.bin.cutoffs <- all.bin.cutoffs[order(all.bin.cutoffs[, c("evaluation", "cell.type", "bin_number")]),]
o <- order(all.bin.cutoffs$evaluation, all.bin.cutoffs$cell.type, all.bin.cutoffs$bin_number)
all.bin.cutoffs <- all.bin.cutoffs[o, ]

write.table(file = paste0(figs.dir, "/", "all-bin-limits.tsv"), all.bin.cutoffs, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


                                                             
nms <- names(raw.bins)
names(nms) <- nms
bin.mats <-
  llply(nms,
        .fun = function(nm) {
          mat <- raw.bins[[nm]]
          if("cell.type" %in% colnames(mat)) {
            ret <- acast(mat, method.name ~ cell.type, value.var = "bin_number")
            colnames(ret) <- paste0(colnames(ret), "_", nm)
          } else {
            ret <- mat[, c("bin_number"), drop=FALSE]
            colnames(ret)[1] <- nm
            rownames(ret) <- mat$method.name
          }
          ret <- data.frame("method.name" = rownames(ret), ret)
          ret
        })

# Order by phenotypes
cell.type.levels <- c(
    "B.cells",
    "memory.B.cells",
    "naive.B.cells",
    "CD4.T.cells",
    "memory.CD4.T.cells",
    "naive.CD4.T.cells",
#    "regulatory.T.cells",
    "Tregs",
    "CD8.T.cells",
    "memory.CD8.T.cells",
    "naive.CD8.T.cells",
    "NK.cells",
    "neutrophils",
    "monocytic.lineage",
    "monocytes",
#    "myeloid.dendritic.cells",
    "myeloid.DCs",    
    "macrophages",
    "endothelial.cells",
    "fibroblasts"
)

cell.type.cols <- gsub(cell.type.levels, pattern=".cells", replacement="")

cell.type.cols <- as.vector(t(outer(cell.type.cols, c("_mean", "_sd", "_lod", "_spillover"), FUN="paste0")))

all.bins <- Reduce(function(...) merge(..., all=TRUE, by="method.name"), bin.mats)
rownames(all.bins) <- all.bins$method.name
all.bins[is.na(all.bins)] <- na.bin
cell.type.cols <- cell.type.cols[(cell.type.cols %in% colnames(all.bins))]
colnames(all.bins) <- c("method.name", cell.type.cols, "sample_coarse", "sample_fine", "timing_coarse", "timing_fine")

cell.type.cols <- grep(colnames(all.bins), pattern="lod|mean|sd|spillover", value=TRUE)
only.do.cell.type.cols <- FALSE
if(only.do.cell.type.cols) {
  all.bins <- all.bins[, c("method.name", cell.type.cols)]
}

# Create a table where NAs in any _cell type_ column are set to zero (i.e., will not)
# contribute to score. 
# First, do this for all methods
all.bins.ignore.nas <- all.bins
# Second, do this only for comparator methods
comparators <- get.comparators.cap()
comparators <- sort(comparators)
all.bins.ignore.nas.comparator <- all.bins
#all.bins.ignore.nas[is.na(all.bins.ignore.nas)] <- 0
# Only ignore NAs in the cell type-specific columns
#stop("stop")
for(col in cell.type.cols) {
  flag <- all.bins.ignore.nas[,col] == na.bin
  all.bins.ignore.nas[flag,col] <- 0
  flag <- (all.bins.ignore.nas.comparator[,col] == na.bin) & (rownames(all.bins.ignore.nas.comparator) %in% comparators)
  all.bins.ignore.nas.comparator[flag,col] <- 0
}
# Also ignore NAs in the cross-sample comparison if the method would otherwise be comparable
# across samples, but was not because not all cell types were reported.
# But only do this for comparator methods. The logic is that contributed/participant methods
# had the opportunity to retrain on all cell types.
# Note the logical: NAs occur because not all cell types are reported and/or the method is scored-based.
# Hence, if the method is not scored-based, we know the NA results from not all cell types being reported, and we
# can ignore the NA
method.annos <- get.method.annotations()
method.annos <- unique(subset(method.annos, is.na(submission) | (submission==1))[, c("method.name", "output.type")])
#method.annos <- subset(method.annos, method.name %in% comparators)
#rownames(method.annos) <- method.annos$method.name
#stopifnot(!any(duplicated(method.annos$method.name)))
score.based.methods <- unique(subset(method.annos, output.type=="score")$method.name)
for(col in c("sample_coarse", "sample_fine")) {
  flag <- (all.bins.ignore.nas[,col] == na.bin) & !(rownames(all.bins.ignore.nas) %in% score.based.methods) & (rownames(all.bins.ignore.nas) %in% comparators)
  all.bins.ignore.nas[flag,col] <- 0
  flag <- (all.bins.ignore.nas.comparator[,col] == na.bin) & !(rownames(all.bins.ignore.nas.comparator) %in% score.based.methods) & (rownames(all.bins.ignore.nas) %in% comparators)
  all.bins.ignore.nas.comparator[flag,col] <- 0
}
for(col in colnames(all.bins.ignore.nas)) { if(any(all.bins.ignore.nas[,col] == na.bin)) { print(col)}}
for(col in colnames(all.bins.ignore.nas.comparator)) { if(any(all.bins.ignore.nas.comparator[,col] == na.bin)) { print(col)}}
#all.bins.ignore.nas[all.bins.ignore.nas[,cell.type.cols] == na.bin, cell.type.cols] <- 0


all.bins <- all.bins[, !(colnames(all.bins) == "method.name")]
all.bins.ignore.nas <- all.bins.ignore.nas[, !(colnames(all.bins.ignore.nas) == "method.name")]
all.bins.ignore.nas.comparator <- all.bins.ignore.nas.comparator[, !(colnames(all.bins.ignore.nas.comparator) == "method.name")]



sorted.method.scores <- sort(rowSums(all.bins))
sorted.method.scores.ignore.nas <- sort(rowSums(all.bins.ignore.nas)/rowSums(all.bins.ignore.nas > 0))
sorted.method.scores.ignore.nas.comparator <- sort(rowSums(all.bins.ignore.nas.comparator)/rowSums(all.bins.ignore.nas.comparator > 0))

sorted.scores <- sorted.method.scores

sorted.methods.alpha <- sorted.scores[order(names(sorted.scores))]


# stop("stop")

if(FALSE) {
# The following code separates coarse and fine-grained cell types and orders alphabetically within each.
# Instead order by phenotypes
cell.type.cols <- unlist(lapply(colnames(all.bins), function(str) strsplit(str, split="_")[[1]][1]))
coarse.cell.types <- gsub(coarse.cell.types, pattern="\\.", replacement = " ")
coarse.cell.types <- gsub(coarse.cell.types, pattern=" cells", replacement = "")
coarse.cell.types <- coarse.cell.types[!(grepl(x=coarse.cell.types, pattern="sample|timing"))]
cell.type.cols <- gsub(cell.type.cols, pattern="\\.", replacement = " ")

coarse.cols <- cell.type.cols %in% coarse.cell.types
fine.cols <- !(cell.type.cols %in% coarse.cell.types) & !(grepl(cell.type.cols, pattern="sample|timing"))

coarse.cell.type.cols <- cell.type.cols[coarse.cols]
fine.cell.type.cols <- cell.type.cols[fine.cols]

coarse.bins <- all.bins[, coarse.cols]
coarse.bins <- coarse.bins[, order(coarse.cell.type.cols)]
fine.bins <- all.bins[, fine.cols]
fine.bins <- fine.bins[, order(fine.cell.type.cols)]

coarse.analysis.cols <- unlist(lapply(colnames(coarse.bins), function(str) strsplit(str, split="_")[[1]][2]))
fine.analysis.cols <- unlist(lapply(colnames(fine.bins), function(str) strsplit(str, split="_")[[1]][2]))

coarse.bins.ignore.nas <- all.bins.ignore.nas[, coarse.cols]
coarse.bins.ignore.nas <- coarse.bins.ignore.nas[, order(coarse.cell.type.cols)]
fine.bins.ignore.nas <- all.bins.ignore.nas[, fine.cols]
fine.bins.ignore.nas <- fine.bins.ignore.nas[, order(fine.cell.type.cols)]

coarse.bins.ignore.nas.comparator <- all.bins.ignore.nas.comparator[, coarse.cols]
coarse.bins.ignore.nas.comparator <- coarse.bins.ignore.nas.comparator[, order(coarse.cell.type.cols)]
fine.bins.ignore.nas.comparator <- all.bins.ignore.nas.comparator[, fine.cols]
fine.bins.ignore.nas.comparator <- fine.bins.ignore.nas.comparator[, order(fine.cell.type.cols)]

if(only.do.cell.type.cols) {
  all.bins <- cbind(coarse.bins, fine.bins)
  all.bins.ignore.nas <- cbind(coarse.bins.ignore.nas, fine.bins.ignore.nas)
  all.bins.ignore.nas.comparator <- cbind(coarse.bins.ignore.nas.comparator, fine.bins.ignore.nas.comparator)
} else {
  all.bins <- cbind(coarse.bins, fine.bins, all.bins[, c("sample_coarse", "sample_fine", "timing_coarse", "timing_fine")])
  all.bins.ignore.nas <- cbind(coarse.bins.ignore.nas, fine.bins.ignore.nas, all.bins.ignore.nas[, c("sample_coarse", "sample_fine", "timing_coarse", "timing_fine")])
  all.bins.ignore.nas.comparator <- cbind(coarse.bins.ignore.nas.comparator, fine.bins.ignore.nas.comparator, all.bins.ignore.nas.comparator[, c("sample_coarse", "sample_fine", "timing_coarse", "timing_fine")])
}
} # end if(FALSE)

if(only.do.cell.type.cols) {
  exclude.cols <- c("sample_coarse", "sample_fine", "timing_coarse", "timing_fine")
  all.bins <- all.bins[, !(colnames(all.bins) %in% exclude.cols)]
  all.bins.ignore.nas <- all.bins.ignore.nas[, !(colnames(all.bins.ignore.nas) %in% exclude.cols)]
  all.bins.ignore.nas.comparator <- all.bins.ignore.nas.comparator[, !(colnames(all.bins.ignore.nas.comparator) %in% exclude.cols)]    
}

cell.type.cols <- unlist(lapply(colnames(all.bins), function(str) strsplit(str, split="_")[[1]][1]))
cell.type.cols <- gsub(cell.type.cols, pattern="\\.", replacement = " ")
cell.type.cols[grepl(x=cell.type.cols, pattern="sample")] <- "Sample-level"
cell.type.cols[grepl(x=cell.type.cols, pattern="timing")] <- "Runtime"
analysis.cols <- unlist(lapply(colnames(all.bins), function(str) strsplit(str, split="_")[[1]][2]))

analysis.cols[analysis.cols == "mean"] <- "Mean"
analysis.cols[analysis.cols == "sd"] <- "SD"
analysis.cols[analysis.cols == "lod"] <- "LoD"
analysis.cols[analysis.cols == "spillover"] <- "Spillover"
analysis.cols[analysis.cols == "coarse"] <- "Coarse"
analysis.cols[analysis.cols == "fine"] <- "Fine"

comparators <- comparators[comparators %in% names(sorted.scores)]

num.nas.per.method <- rowSums(all.bins==5)
# Note that some fine-grained methods do not report scores that are comparable
# across cell types. This will lead to NAs.
# Also, some crashed when we ran across the cancer data sets.
#fine.grained.methods <- names(num.nas.per.method[num.nas.per.method<15])
#fine.grained.methods <- fine.grained.methods[!(fine.grained.methods %in% comparators)]

#coarse.grained.methods <- rownames(all.bins)
#coarse.grained.methods <- coarse.grained.methods[!(coarse.grained.methods %in% c(comparators, fine.grained.methods))]

num.nas.per.method.fine <- rowSums(fine.bins==5)
coarse.grained.methods <- names(num.nas.per.method.fine[num.nas.per.method.fine == ncol(fine.bins)])
fine.grained.methods <- rownames(all.bins)
fine.grained.methods <- fine.grained.methods[!(fine.grained.methods %in% c(comparators, coarse.grained.methods))]

fine.grained.methods <- sort(fine.grained.methods)
coarse.grained.methods <- sort(coarse.grained.methods)
sorted.methods.by.class.alpha <- sorted.scores[c(comparators, fine.grained.methods, coarse.grained.methods)]



cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colors = c(cbbPalette[2:5], "#000000")
library(RColorBrewer)
cat.colors <- cbbPalette[2:(1+length(unique(analysis.cols)))]
cat.names <- c("Mean", "SD", "LoD", "Spillover", "Coarse", "Fine")
if(only.do.cell.type.cols) {
  cat.names <- c("Mean", "SD", "LoD", "Spillover")
}
names(cat.colors) <- cat.names

name.order <- names(sorted.scores)

make.heatmap <- function(data, name.order) {
  row.labels <- rownames(data[name.order,])
  flag <- row.labels %in% get.comparators.cap()
  # Make the comparator methods bold
  row.labels[flag] <- paste0("**", row.labels[flag], "**")
  
  colors <- c(rev(brewer.pal(n=4, name='RdBu')), "#000000")
  names(colors) <- as.character(1:5)
  
  if(any(data == 0)) {
    colors <- c(colors, "#FFFFFF")
    data[data == 0] <- "NA"
    names(colors) <- c(as.character(1:5),"NA")
    print(colors)
  }
  
  htmp <-
  Heatmap(data[name.order,], cluster_rows = FALSE, cluster_columns = FALSE, col = colors,
          name = "Bin",
          column_split = factor(cell.type.cols, levels = unique(cell.type.cols)), column_title_rot = 45,
          top_annotation = HeatmapAnnotation(Category=factor(analysis.cols, levels = cat.names), 
                                             annotation_legend_param = list(labels = cat.names),
                                             col = list(Category = cat.colors),
                                             annotation_name_gp= gpar(fontsize = 8)),
          show_column_names = FALSE, 
          column_title_gp = gpar(fontsize = 8),
          column_names_gp = grid::gpar(fontsize = 8),
          row_names_gp = grid::gpar(fontsize = 8),
          row_labels = gt_render(row.labels))
}

all.sorted.scores <- list("score-sorted-ignore-nas" = sorted.method.scores.ignore.nas,
            "score-sorted-ignore-nas-comparator" = sorted.method.scores.ignore.nas.comparator,
            "score-sorted" = sorted.scores,
            "alpha" = sorted.methods.alpha,
            "alpha-by-class" = sorted.methods.by.class.alpha)

all.data <- list("score-sorted-ignore-nas" = all.bins.ignore.nas,
                          "score-sorted-ignore-nas-comparator" = all.bins.ignore.nas.comparator,
                          "score-sorted" = all.bins,
                          "alpha" = all.bins,
                          "alpha-by-class" = all.bins)


for(nm in names(all.sorted.scores)) {
  htmp <- make.heatmap(all.data[[nm]], names(all.sorted.scores[[nm]]))

  png(paste0(figs.dir, "/", "binned-score-heatmap-", nm, ".png"))
  draw(htmp, merge_legend=TRUE)
  d <- dev.off()

  pdf(paste0(figs.dir, "/", "binned-score-heatmap-", nm, ".pdf"))
  draw(htmp, merge_legend=TRUE)
  d <- dev.off()
  
  cols <- colnames(all.bins)
  tot.df <- data.frame("total.score" = as.numeric(all.sorted.scores[[nm]]))
  rownames(tot.df) <- names(all.sorted.scores[[nm]])
  all.bins.with.tot.scores <- merge(all.bins, tot.df, by = "row.names")
  colnames(all.bins.with.tot.scores)[1] <- "method.name"
  all.bins.with.tot.scores <- all.bins.with.tot.scores[order(all.bins.with.tot.scores$total.score, decreasing=FALSE),]
  
  write.table(file = paste0(figs.dir, "/", "binned-scores-", nm, ".tsv"), all.bins.with.tot.scores, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
}


cat(paste("Exiting successfully\n"))
stop("stop")
q(status=0)

