suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
if(!suppressPackageStartupMessages(require(rlang))) {
  install.packages("rlang")
}
if(!suppressPackageStartupMessages(require(purrr))) {
  install.packages("purrr")
}
suppressPackageStartupMessages(p_load(purrr))
suppressPackageStartupMessages(p_load(ggplot2))
# suppressPackageStartupMessages(p_load(tidyverse))
suppressPackageStartupMessages(p_load(ggupset)) # for axis_combmatrix
#if(!suppressPackageStartupMessages(require(ComplexUpset))) {
#  devtools::install_github("krassowski/complex-upset")
#}
#suppressPackageStartupMessages(p_load(ComplexUpset))
suppressPackageStartupMessages(p_load(reshape2))
suppressPackageStartupMessages(p_load(cowplot))

set.seed(1234)

synLogin()

cat(paste0("Run on JAX HPC with conda env r3\n"))

method.name.col <- "method.name"

figs.dir <- "figs/"
dir.create(figs.dir, showWarnings = FALSE)

## Read in the rerun predicitons (i.e., where the coarse- and fine-grained datasets are the same)
synId <- "syn22320329"
obj <- synGet(synId, downloadFile=TRUE)
res.all <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)
res.all <- subset(res.all, submission == "1")

# this defines methods.to.exclude
source("../validation-analysis/methods-to-exclude.R")

# Let's also exclude Aboensis IV, as it crashed on the fine-grained simulated cancer datasets
# Error: Error in predict.randomForest(rf, testdata) : missing values in newdata
# Maybe it was expecting all fine-grained Challenge cell types, but some were missing from each dataset.
# Also exclude this from the cancer datasets
methods.to.exclude <- c(methods.to.exclude, "Aboensis IV")

flag <- res.all[,method.name.col] %in% methods.to.exclude
cat(paste0("Excluding methods: ", paste(unique(res.all[flag, method.name.col]), collapse = ", "), "\n"))
res.all <- res.all[!flag,]

comparator.methods <- unique(subset(res.all, comparator==TRUE)[, "method.name"])

# The above results already have ground truth appended ('measured' col).
# However, there were two versions of ground truth (and, of course, admixtures):
# 1. The original admixtures were different between fine- and coarse-grained
# 2. The revised and final admixtures were the same between fine- and coarse-grained.
# Let's make the results have those second set of final/revised ground truth values

# Original coarse-grained ground truth
# coarse.gt.synId <- "syn21820375"
# New coarse-grained ground truth based on fine-grained
coarse.gt.synId <- "syn22267267"
fine.gt.synId <- "syn21820376"
coarse.gt <- read.table(synGet(coarse.gt.synId)$path, header=TRUE, sep=",")
coarse.gt$subchallenge <- "coarse"
fine.gt <- read.table(synGet(fine.gt.synId)$path, header=TRUE, sep=",")
fine.gt$subchallenge <- "fine"
gt <- rbind(coarse.gt, fine.gt)
colnames(gt)[colnames(gt) == "measured"] <- "revised.measured"

m <- merge(res.all, gt, all.x=TRUE)
m <- m[!is.na(m$revised.measured),]
stopifnot(all(m$measured == m$revised.measured))

method.id.to.name.tbl <- unique(res.all[, c("objectId", "method.name")])

res.all <- res.all[, c("dataset.name", "subchallenge", "method.name", "cell.type", "sample.id", "prediction", "objectId", "measured")]

coarse.orig.preds <- subset(res.all, subchallenge == "coarse")
fine.orig.preds <- subset(res.all, subchallenge == "fine")

## Read in tumor predictions of submitted models
synId <- "syn51238861"
obj <- synGet(synId, downloadFile=TRUE)
tumor.res.all <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)
# I think the round field here is broken.
# Instead, match these results by objectId to those above
# This will effectively exclude Aboensis IV
tumor.res.all <- subset(tumor.res.all, objectId %in% res.all$objectId)



# Read in tumor predictions of comparator models
synId <- "syn51273581"
obj <- synGet(synId, downloadFile=TRUE)
tumor.res.comp.all <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)



## Keep only the first submission
#tumor.res.all <- subset(tumor.res.all, round == "1")

# Effectively rename team to method.name to be consistent with the above
tumor.res.all <- tumor.res.all[, c("dataset.name", "subchallenge", "team", "cell.type", "sample.id", "prediction", "objectId", "measured")]
colnames(tumor.res.all) <- c("dataset.name", "subchallenge", "method.name", "cell.type", "sample.id", "prediction", "objectId", "measured")

tumor.res.comp.all <- tumor.res.comp.all[, c("dataset.name", "subchallenge", "team", "cell.type", "sample.id", "prediction", "objectId", "measured")]
colnames(tumor.res.comp.all) <- c("dataset.name", "subchallenge", "method.name", "cell.type", "sample.id", "prediction", "objectId", "measured")

# Name the methods based on objectId
method.names <- 
  list("cibersort" = "CIBERSORT", "epic" = "EPIC", "mcpcounter" = "MCP-counter", "quantiseq" = "quanTIseq",
       "timer" = "TIMER", "xcell" = "xCell")
for(mn in names(method.names)) {
  flag <- grepl(tumor.res.comp.all$objectId, pattern=mn)
  tumor.res.comp.all[flag, "method.name"] <- method.names[[mn]]
}


# Note that the 'team' column in these results seems to be the name of the submitter.
# Let's translate those, as we did with the original results, to be the method.name
# That is just shortening names, e.g., "Northwestern Polytechnical University" -> "NPU"
source("../utils.R")
tumor.res.all <- simplify.submitter.names(tumor.res.all, "method.name")

tumor.res.all <- rbind(tumor.res.all, tumor.res.comp.all)

# Restrict the tumor datasets to those that only mixed challenge cell types (+ tumor) together,
# as opposed to mixing unknown cell types. This is most consistent with what we did in the
# in vitro and in silico challenge.
tumor.res.all <- subset(tumor.res.all, dataset.name %in% c("pelka-only.challenge.types", "wu-challenge-cells"))
flag <- tumor.res.all$dataset.name == "pelka-only.challenge.types"
tumor.res.all[flag, "dataset.name"] <- "Pelka"
flag <- tumor.res.all$dataset.name == "wu-challenge-cells"
tumor.res.all[flag, "dataset.name"] <- "Wu"

flag <- grepl(tumor.res.all$subchallenge, pattern="coarse")
tumor.res.all[flag, "subchallenge"] <- "coarse"
flag <- grepl(tumor.res.all$subchallenge, pattern="fine")
tumor.res.all[flag, "subchallenge"] <- "fine"

tumor.res.all$subchallenge <- gsub(x=tumor.res.all$subchallenge, pattern="*.coarse", replacement="coarse", perl=FALSE)
tumor.res.all$subchallenge <- gsub(x=tumor.res.all$subchallenge, pattern="*.fine", replacement="fine", perl=FALSE)

tumor.res.all <- tumor.res.all[!is.na(tumor.res.all$measured),]
res.all <- res.all[!is.na(res.all$measured),]

tumor.gt <- unique(tumor.res.all[, c("dataset.name", "subchallenge", "cell.type", "sample.id", "measured")]) 

# Recall that in the cancer admixtures, some populations were not included in the single-cell data.
# These will show up as having all measured values == 0.
# This will give us NAs when we try to compute correlations (and are dummy variables anyways)
# Detect and exclude them
tmp <-
  ddply(tumor.gt,
        .variables = c("dataset.name", "subchallenge", "cell.type"),
        .fun = function(df) {
          data.frame(all.zero = all(df$measured == 0))
        })
non.zero.populations <- subset(tmp, all.zero==FALSE)
tumor.res.all <- merge(tumor.res.all, non.zero.populations)
tumor.res.all <- tumor.res.all[, !(colnames(tumor.res.all) %in% "all.zero")]

## Read in all the CIBERSORTx results, merge gt, and append to above tumor.res.all results

# Read in the tumor predictions for cibersortX
# Pelka/CRC: https://www.synapse.org/#!Synapse:syn51277135
# Pelka/CRC coarse
synId <- "syn51277165"
obj <- synGet(synId, downloadFile=TRUE)
csx.pelka.crc.coarse <- read.table(obj$path, sep="\t", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)

# Pelka/CRC fine
synId <- "syn51277164"
obj <- synGet(synId, downloadFile=TRUE)
csx.pelka.crc.fine <- read.table(obj$path, sep="\t", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)

# Wu/BRCA:  https://www.synapse.org/#!Synapse:syn51277079
# Wu/BRCA coarse
synId <- "syn51277148"
obj <- synGet(synId, downloadFile=TRUE)
csx.wu.brca.coarse <- read.table(obj$path, sep="\t", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)

# Wu/BRCA fine
synId <- "syn51277133"
obj <- synGet(synId, downloadFile=TRUE)
csx.wu.brca.fine <- read.table(obj$path, sep="\t", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)

csx.all <- rbind(csx.pelka.crc.coarse, csx.pelka.crc.fine, csx.wu.brca.coarse, csx.wu.brca.fine)
csx.all$objectId <- "CIBERSORTx"

csx.all <- subset(csx.all, dataset.name %in% c("pelka-only.challenge.types", "wu-challenge-cells"))
flag <- csx.all$dataset.name == "pelka-only.challenge.types"
csx.all[flag, "dataset.name"] <- "Pelka"
flag <- csx.all$dataset.name == "wu-challenge-cells"
csx.all[flag, "dataset.name"] <- "Wu"

csx.all <- merge(csx.all, non.zero.populations)
nr <- nrow(csx.all)
csx.all <- merge(csx.all, unique(tumor.res.all[, c("dataset.name", "sample.id", "cell.type", "subchallenge", "measured")]))
stopifnot(nr == nrow(csx.all))

csx.all <- csx.all[, colnames(tumor.res.all)]
tumor.res.all <- rbind(tumor.res.all, csx.all)

# Memory B cells are not in any of the original data, so exclude from tumor data
tumor.res.all <- subset(tumor.res.all, !(cell.type == "memory.B.cells"))

stopifnot(all(sort(unique(tumor.res.all$method.name)) == sort(unique(tumor.res.all$method.name))))

tumor.scores <-
  ddply(tumor.res.all, .variables = c("method.name", "subchallenge"),
        .fun = function(method.preds) {
          score.preds(method.preds)$all
        })

orig.scores <-
  ddply(res.all, .variables = c("method.name", "subchallenge"),
        .fun = function(method.preds) {
          score.preds(method.preds)$all
        })

common.methods <- intersect(unique(tumor.scores$method.name), unique(orig.scores$method.name))

tumor.scores <- subset(tumor.scores, method.name %in% common.methods)
orig.scores <- subset(orig.scores, method.name %in% common.methods)

all.scores <- rbind(tumor.scores, orig.scores)



in.silico.validation.metadata.synId <- "syn22013519"
obj <- synGet(in.silico.validation.metadata.synId, downloadFile=TRUE)
in.silico.validation.metadata <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)
in.silico.validation.metadata <- in.silico.validation.metadata[, c("tumor.type", "mixture.type", "dataset.name")]
colnames(in.silico.validation.metadata) <- c("tumor.type", "distribution.type", "dataset.name")
in.silico.validation.metadata$mixture.type <- "In Silico"

suppressPackageStartupMessages(p_load("openxlsx"))
in.vitro.validation.metadata <- get.validation.metadata()
in.vitro.validation.metadata <- unique(in.vitro.validation.metadata[, c(2,4,5)])
colnames(in.vitro.validation.metadata) <- c("tumor.type", "distribution.type", "dataset.name")
in.vitro.validation.metadata$mixture.type <- "In Vitro"

validation.metadata <- rbind(in.silico.validation.metadata, in.vitro.validation.metadata)

validation.metadata$immune.origin <- "Healthy"
validation.metadata <- rbind(validation.metadata, c("BRCA", "Biological", "Wu", "In Silico", "Tumor"))
validation.metadata <- rbind(validation.metadata, c("CRC", "Biological", "Pelka", "In Silico", "Tumor"))

validation.metadata <- subset(validation.metadata, (dataset.name %in% res.all$dataset.name) | (dataset.name %in% tumor.res.all$dataset.name) )

validation.metadata <- validation.metadata %>%
  mutate(Label = pmap(list(tumor.type, distribution.type, mixture.type, immune.origin), 
                      function(tumor.type, distribution.type, mixture.type, immune.origin){
                        c(tumor.type, distribution.type, mixture.type, immune.origin)
                      }),
         Label.str = unlist(pmap(list(tumor.type, distribution.type, mixture.type, immune.origin), 
                          function(tumor.type, distribution.type, mixture.type, immune.origin){
                            paste(c(tumor.type, distribution.type, mixture.type, immune.origin), collapse="-")
                          })))


# Average scores across sub-challenges for cell types common to both
all.scores <-
  ddply(all.scores, .variables = c("method.name", "cell.type", "dataset.name"),
        .fun = function(df) {
          data.frame(cor.p = mean(df$cor.p),
                     cor.s = mean(df$cor.s),
                     rmse = mean(df$rmse))
        })

# Drop any NA
na.scores <-
  ddply(all.scores,
        .variables = c("method.name", "cell.type"),
        .fun = function(df) {
          data.frame(na.scores = any(is.na(df$cor.p)))
        })

all.scores <- merge(all.scores, na.scores)
all.scores <- subset(all.scores, na.scores == FALSE)

all.scores <- merge(all.scores, validation.metadata)

compute.stats_ <- function(df) {
  per.cell.mean <-
    ddply(df, 
          .variables = c("cell.type"),
          .fun = function(df.ct) {
            data.frame(mean = mean(na.omit(subset(df.ct, grepl(dataset.name, pattern="AA"))$cor.p)))
          })
  df <- merge(df, per.cell.mean)
  print(per.cell.mean)
  df$cor.p <- df$cor.p - df$mean
  
  ddply(df,
        .variables = c("cell.type"),
        .fun = function(df.ct) {
          print(df.ct[1,"cell.type"])
          all.factors <- c("mixture.type", "distribution.type", "dataset.name", "tumor.type")
          flags <- unlist(lapply(all.factors, function(fct) length(unique(df.ct[,fct])) > 1))
          factors <- all.factors[flags]
          mdl.str <- paste0("cor.p ~ ", paste(factors, collapse=" + "))
          mdl.str <- paste0("cor.p ~ 0 + dataset.name + method.name")
          # mdl.str <- paste0("cor.p ~ 0 + dataset.name")
          cat(paste0(df.ct[1,"cell.type"], ": ", mdl.str, "\n"))
          mdl <- as.formula(mdl.str)
          lm.fit <- lm(mdl, data=df.ct)
          sm <- summary(lm.fit)
          cf <- coef(sm)
          #flag <- grepl(rownames(cf), pattern="type") | grepl(rownames(cf), pattern="dataset") | grepl(rownames(cf), pattern="method")
          #ret.df <- cf[flag,]
          ret.df <- cf
          ret.df <- cbind(variable = rownames(ret.df), ret.df)
          pval <- pf(sm$fstatistic[1],sm$fstatistic[2],sm$fstatistic[3],lower.tail=FALSE)
          ret.df <- rbind(ret.df, c("F-statistic", as.numeric(sm$fstatistic[1]), NA, NA, pval))
          ret.df
        })
}

# Compare stats by controlling for method
compute.stats <- function(df) {
  # Motivation: compute significant differences from per-method, per-cell type mean across datasets
  # to see what datasets are above or below that mean.
  # 1. Control for method by subtracting off the method mean (for each cell type) across datasets
  ## No! Method mean is just zero (after subtracting off method/cell type mean).
  ## Also, this doesn't make much intuitive sense.
  ## 2. Subtract off overall (method-controlled) mean (across methods for each cell type)
  # 3. Use a linear model _without_ intercept (since the mean goes through zero)

  # Compute mean across datasets for each method and cell type
  per.method.cell.type.mean <-
    ddply(df, 
          .variables = c("method.name", "cell.type"),
          .fun = function(df.m) {
            data.frame(method.cell.type.mean = mean(df.m$cor.p))
          })

  print(head(per.method.cell.type.mean))
  #per.method.cell.type.mean <- merge(per.method.cell.type.mean, per.cell.type.mean)
  
  # Subtract the mean across datasets for each method and cell type from the
  # per-method, per-cell type, per-dataset score.
  df <- merge(df, per.method.cell.type.mean)
  df$cor.p.method.adj <- df$cor.p - df$method.cell.type.mean
  
  subtract.overall.mean <- FALSE
  cor.p.var <- "cor.p.method.adj"
  if(subtract.overall.mean) {
    per.cell.type.mean <-
      ddply(df,
            .variables = c("cell.type"),
            .fun = function(df.ct) {
              data.frame(cell.type.mean = mean(df.ct$cor.p.method.adj))
            })
    
    # Duh! These cell.type.means are zero, of course.
    print(head(per.cell.type.mean))
    
    df <- merge(df, per.cell.type.mean)
    df$cor.p.cell.type.method.adj <- df$cor.p.method.adj - df$cell.type.mean
    cor.p.var <- "cor.p.cell.type.method.adj"
  }
  
  ret <-
    ddply(df,
          .variables = c("cell.type"),
          .parallel = TRUE,
          .fun = function(df.ct) {
            print(df.ct[1,"cell.type"])
            all.factors <- c("mixture.type", "distribution.type", "dataset.name", "tumor.type")
            flags <- unlist(lapply(all.factors, function(fct) length(unique(df.ct[,fct])) > 1))
            factors <- all.factors[flags]
            mdl.str <- paste0(cor.p.var, " ~ ", paste(factors, collapse=" + "))
            # mdl.str <- paste0(cor.p.var, " ~ 0 + dataset.name + method.name")
            mdl.str <- paste0(cor.p.var, " ~ 0 + dataset.name")
            cat(paste0(df.ct[1,"cell.type"], ": ", mdl.str, "\n"))
            mdl <- as.formula(mdl.str)
            lm.fit <- lm(mdl, data=df.ct)
            sm <- summary(lm.fit)
            cf <- coef(sm)
	    # Add the 2.5-97.5% confidence interval
	    cf <- cbind(cb, confint(lm.fit))
            
            # One-sided p-values are discussed here:
            # https://stats.stackexchange.com/questions/325354/if-and-how-to-use-one-tailed-testing-in-multiple-regression
            
            # For H1: beta < 0
            cf <- cbind(cf, "dof" = lm.fit$df, "one.sided.pval" = unlist(lapply(as.numeric(cf[, 3]), function(x) pt(x, lm.fit$df, lower = TRUE))))

            #flag <- grepl(rownames(cf), pattern="type") | grepl(rownames(cf), pattern="dataset") | grepl(rownames(cf), pattern="method")
            #ret.df <- cf[flag,]
            ret.df <- as.data.frame(cf)
            ret.df <- cbind(variable = rownames(ret.df), ret.df)
            pval <- pf(sm$fstatistic[1],sm$fstatistic[2],sm$fstatistic[3],lower.tail=FALSE)
            ret.df <- rbind(ret.df, c("F-statistic", as.numeric(sm$fstatistic[1]), NA, NA, NA, pval))
            ret.df
          })
  return(list("stats"=ret, "adj.df" = df))
}


ret <- compute.stats(all.scores)
all.stats <- ret[["stats"]]
colnames(all.stats)[colnames(all.stats) == "Pr(>|t|)"] <- "two.sided.pval"
adj.df <- ret[["adj.df"]]
adj.df <- adj.df[, !(colnames(adj.df) %in% c("Label"))]

# Compute median across methods for all datasets
median.stats <- 
  ddply(adj.df,
        .variables = c("cell.type", "dataset.name", "tumor.type", "distribution.type", "mixture.type", "immune.origin", "Label.str"),
        .fun = function(df) {
                 data.frame(cor.p = median(df$cor.p, na.rm=TRUE),
                            cor.s = median(df$cor.s, na.rm=TRUE),
                            rmse = median(df$rmse, na.rm=TRUE))
        })


write.table(file = paste0(figs.dir, "/cancer-validation-dataset-pvals.tsv"), all.stats, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(file = paste0(figs.dir, "/cancer-validation-median-correlations.tsv"), median.stats, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(file = paste0(figs.dir, "/cancer-validation-correlations.tsv"), adj.df, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


# Pick out stats for Pelka, Wu, and overall.
# We should only look at Pelka and Wu if overall is significant.
overall.stats <- subset(all.stats, grepl(variable, pattern="F-statistic"))
overall.stats <- overall.stats[, c("cell.type", "Estimate", "two.sided.pval")]
colnames(overall.stats) <- c("cell.type", "F.stat.estimate", "F.stat.p")

all.dataset.stats <- subset(all.stats, grepl(variable, pattern="dataset.name"))
p.val.col <- "two.sided.pval"
p.val.col <- "one.sided.pval"
all.dataset.stats$two.sided.pval.holm <- p.adjust(all.dataset.stats[,"two.sided.pval"], method="holm")
all.dataset.stats$one.sided.pval.holm <- p.adjust(all.dataset.stats[,"one.sided.pval"], method="holm")
wu.pelka.stats <- subset(all.dataset.stats, grepl(variable, pattern="Pelka") | grepl(variable, pattern="Wu"))
non.wu.pelka.stats <- subset(all.dataset.stats, !grepl(variable, pattern="Pelka") & !grepl(variable, pattern="Wu"))
sig.dataset.stats <- subset(all.dataset.stats, two.sided.pval.holm < 0.05)
print(sig.dataset.stats[order(sig.dataset.stats$variable),c("cell.type", "variable", "Estimate", "two.sided.pval.holm")])
write.table(file = paste0(figs.dir, "/cancer-validation-dataset-comparison-pvals.tsv"), all.dataset.stats, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


wu.pelka.stats <- merge(wu.pelka.stats, overall.stats)
wu.pelka.stats[, p.val.col] <- as.numeric(wu.pelka.stats[, p.val.col])
wu.pelka.stats$F.stat.p <- as.numeric(wu.pelka.stats$F.stat.p)
wu.pelka.stats$Estimate <- as.numeric(wu.pelka.stats$Estimate)
stopifnot(all(wu.pelka.stats$F.stat.p < 0.05))
cat(paste0("All Anovas are significant for Wu/Pelka datasets max pvalue: ", max(wu.pelka.stats$F.stat.p), "\n"))

non.wu.pelka.stats <- merge(non.wu.pelka.stats, overall.stats)
non.wu.pelka.stats[, p.val.col] <- as.numeric(non.wu.pelka.stats[, p.val.col])
non.wu.pelka.stats$F.stat.p <- as.numeric(non.wu.pelka.stats$F.stat.p)
non.wu.pelka.stats$Estimate <- as.numeric(non.wu.pelka.stats$Estimate)
# Several are not significant
#stopifnot(all(non.wu.pelka.stats$F.stat.p < 0.05))
#cat(paste0("All Anovas are significant for non Wu/Pelka max pvalue: ", max(non.wu.pelka.stats$F.stat.p), "\n"))


all.stats <- merge(all.stats, overall.stats)


# axis_combmatrix respects the order of levels, whereas scale_x_upset sorts the columns based on freq, etc.
# validation.metadata <- validation.metadata[order(validation.metadata$immune.origin, validation.metadata$tumor.type, validation.metadata$distribution.type),]
validation.metadata <- validation.metadata[order(validation.metadata$immune.origin, validation.metadata$mixture.type, validation.metadata$tumor.type, validation.metadata$distribution.type),]

dataset.levels <- validation.metadata$dataset.name
label.str.levels <- validation.metadata$Label.str


all.scores$dataset.name <- factor(all.scores$dataset.name, levels = dataset.levels)
all.scores <- all.scores[order(all.scores$dataset.name),]
all.scores$Label.str <- factor(all.scores$Label.str, levels = label.str.levels)



sets <- c("BRCA", "CRC", "Biological", "Random", "In Vitro", "In Silico", "Healthy", "Tumor")
sets <- c("Healthy", "Tumor", "In Vitro", "In Silico", "BRCA", "CRC", "Biological", "Random")

all.scores <- rename.cell.types(all.scores, from.col = "cell.type", to.col = "cell.type")

# Order cell types according to phenotype, rather than alphabetically
# No! Order by a quantitative metric
#cell.type.levels <- c(
#"B", "naive B", 
#"CD8 T", "naive CD8 T", "memory CD8 T", 
#"CD4 T", "naive CD4 T", "memory CD4 T", "Tregs", 
#"NK", "neutrophils", 
#"monocytic lineage", "monocytes", "myeloid DCs", "macrophages", 
#"endothelial", "fibroblasts")

means.across.all.8.datasets <-
  ddply(all.scores,
        .variables = c("method.name", "cell.type"),
        .fun = function(df) {
          data.frame(cor.p = mean(df$cor.p))
        })
means.across.all.8.datasets$dataset.name <- "Mean"

cell.type.maxs.across.all.8.datasets <- ddply(means.across.all.8.datasets, .variables=c("cell.type"), .fun = function(df) data.frame(cor.p = max(df$cor.p)))
cell.type.maxs.across.all.8.datasets <- cell.type.maxs.across.all.8.datasets[order(cell.type.maxs.across.all.8.datasets$cor.p,decreasing=TRUE),]

#method.means.across.all.8.datasets <- ddply(means.across.all.8.datasets, .variables=c("method.name"), .fun = function(df) data.frame(cor.p = mean(df$cor.p)))
#method.means.across.all.8.datasets <- method.means.across.all.8.datasets[order(method.means.across.all.8.datasets$cor.p,decreasing=FALSE),]


#cell.type.levels <- calculate.cell.type.levels(all.scores, id.var = "method.name", cell.type.var = "cell.type", cor.var = "cor.p")
#print(cell.type.levels)
cell.type.levels <- cell.type.maxs.across.all.8.datasets$cell.type
print(cell.type.levels)
#method.levels <- method.means.across.all.8.datasets$method.name


all.scores$cell.type <- factor(all.scores$cell.type, levels = rev(cell.type.levels))
#all.scores$method.name <- factor(all.scores$method.name, levels = method.levels)
g <- ggplot(all.scores, aes(x=Label.str, y=cor.p)) + geom_boxplot() 
# g <- ggplot(all.scores, aes(x=Label, y=cor.p)) + geom_boxplot() 
# g <- g + geom_jitter(width=0.1) 
# g <- g + scale_x_upset(position = "top", sets = sets)
g <- g + axis_combmatrix(sep = "-", levels = sets) + scale_x_discrete(position = "top")
g <- g + theme_combmatrix(combmatrix.panel.line.size=0)
g <- g + facet_wrap(~ cell.type)
g <- g + xlab("Dataset") + ylab("Pearson Correlation")

png(paste0(figs.dir, "fig-cancer-validation", ".png"), width = 1 * 480, height = 1 * 480)                    
print(g)
d <- dev.off()

pdf(paste0(figs.dir, "fig-cancer-validation", ".pdf"), width = 1 * 7, height = 1 * 7)                    
print(g)
d <- dev.off()

all.scores.simplified <- all.scores
all.scores.simplified <- rename.cell.types(all.scores.simplified, from.col = "cell.type", to.col = "cell.type")
all.scores.simplified$dataset.name <- as.character(all.scores.simplified$dataset.name)
flag <- all.scores.simplified$dataset.name %in% c("Pelka", "Wu")
all.scores.simplified[!flag,"dataset.name"] <- "Healthy"
flag <- all.scores.simplified$dataset.name %in% c("Wu")
all.scores.simplified[flag,"dataset.name"] <- "Wu (BRCA)"
flag <- all.scores.simplified$dataset.name %in% c("Pelka")
all.scores.simplified[flag,"dataset.name"] <- "Pelka (CRC)"
all.scores.simplified <- 
  ddply(all.scores.simplified, .variables = c("dataset.name", "method.name", "cell.type"),
        .fun = function(df) {
          data.frame(cor.p = mean(df$cor.p))
        })

# g <- plot.strip.plots(all.scores.simplified, id.var="method.name", cell.type.var="cell.type", var="cor.p")

means <-
  ddply(all.scores.simplified,
        .variables = c("method.name", "cell.type"),
        .fun = function(df) {
          data.frame(cor.p = mean(df$cor.p))
        })
means$dataset.name <- "Mean"

flag <- duplicated(all.scores.simplified[, c("method.name", "cell.type")], fromLast = TRUE) |
  duplicated(all.scores.simplified[, c("method.name", "cell.type")], fromLast = FALSE) 
# For those cell types / methods for which we have results across multiple datassets, calculate the mean across those datasets.
# We will only show the mean when there is more than one dataset. This does not impact the calculation of any stats.
means.of.repeated <-
  ddply(all.scores.simplified[flag,],
        .variables = c("method.name", "cell.type"),
        .fun = function(df) {
          data.frame(cor.p = mean(df$cor.p))
        })
means.of.repeated$dataset.name <- "Mean"

cell.type.means <- ddply(means, .variables=c("cell.type"), .fun = function(df) data.frame(cor.p = mean(df$cor.p)))
cell.type.means <- cell.type.means[order(cell.type.means$cor.p,decreasing=TRUE),]

# No typo here -- order cell types by their max (over methods) of their means (over datasets)
cell.type.maxs <- ddply(means, .variables=c("cell.type"), .fun = function(df) data.frame(cor.p = max(df$cor.p)))
cell.type.maxs <- cell.type.maxs[order(cell.type.maxs$cor.p,decreasing=TRUE),]

method.means <- ddply(means, .variables=c("method.name"), .fun = function(df) data.frame(cor.p = mean(df$cor.p)))
method.means <- method.means[order(method.means$cor.p,decreasing=FALSE),]

best.performing <-
  ddply(means,
        .variables = c("cell.type"),
        .fun = function(df) {
          df[which.max(df$cor.p),]
        })

subset.methods <- 
  unique(c(get.comparators.cap(), best.performing$method.name, c("Aginome-XMU", "Biogem", "mitten_TDC19", "DA_505")))
subset.methods <- subset.methods[!(subset.methods == "TIMER")]

all.scores.simplified <- rbind(all.scores.simplified[, c("dataset.name", "method.name", "cell.type", "cor.p")],
                               means.of.repeated[, c("dataset.name", "method.name", "cell.type", "cor.p")])

all.scores.simplified$dataset.name <- factor(all.scores.simplified$dataset.name, levels = c("Healthy", "Pelka (CRC)", "Wu (BRCA)", "Mean"))
# all.scores.simplified <- subset(all.scores.simplified, dataset.name == "Mean")

# Order cell types by their max (over methods) of their means (over datasets)
# all.scores.simplified$cell.type <- factor(all.scores.simplified$cell.type, levels = cell.type.means$cell.type)
all.scores.simplified$cell.type <- factor(all.scores.simplified$cell.type, levels = rev(cell.type.maxs$cell.type))
# Order methods by their mean (over cell type) of their means (over datasets)
all.scores.simplified$method.name <- factor(all.scores.simplified$method.name, levels = method.means$method.name)
# Append means

df <- subset(all.scores.simplified, method.name %in% subset.methods)

lvls <- levels(df[,method.name.col])
lvls <- lvls[lvls %in% df[,method.name.col]]
y.bold.labels <- ifelse(lvls %in% comparator.methods, yes = "bold", no = "plain")

datasets <- as.character(unique(df$dataset.name))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
g <- ggplot(data = df,
            aes(x = cor.p, y = method.name, colour = dataset.name)) + facet_wrap(~ cell.type, nrow=3)
g <- g + geom_point() + ylab("") + xlab("Pearson Correlation") + labs(colour = "Dataset")
g <- g + theme(axis.text.y = element_text(face = y.bold.labels))
g <- g + scale_color_manual(breaks = c(datasets[datasets != "Mean"], "Mean"),
                            values=c(cbbPalette[2:length(datasets)], "#000000"))

strip.text.sz <- 8 
g <- g + theme(axis.text.y = element_text(size = 6), text = element_text(size=15), strip.text = element_text(size = strip.text.sz))
g <- g + theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 15), strip.text = element_text(size = 15))
g <- g +  theme(legend.position="top")

png(paste0(figs.dir, "fig-cancer-validation-per-cell-type.png"), width = 2 * 480, height = 1 * 480)                    
print(g)
d <- dev.off()

pdf(paste0(figs.dir, "fig-cancer-validation-per-cell-type.pdf"), width = 2 * 7, height = 1 * 7)                    
print(g)
d <- dev.off()

all.scores <- all.scores[, !(colnames(all.scores) %in% "Label")]

write.table(file = paste0(figs.dir, "/cancer-validation-dataset-all-scores.tsv"), all.scores, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


all.scores.wu <- subset(all.scores, dataset.name=="Wu")
all.scores.pelka <- subset(all.scores, dataset.name=="Pelka")

all.scores.cancer.merged <- rbind(all.scores.wu, all.scores.pelka)
all.scores.cancer.merged <-
  ddply(all.scores.cancer.merged,
        .variables = c("method.name", "cell.type"),
        .fun = function(df) {
          data.frame(cor.p = mean(df$cor.p, na.rm=TRUE),
                     cor.s = mean(df$cor.s, na.rm=TRUE),
                     rmse = mean(df$rmse, na.rm=TRUE))
        })

write.table(file = paste0(figs.dir, "/cancer-validation-dataset-wu-scores.tsv"), all.scores.wu, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(file = paste0(figs.dir, "/cancer-validation-dataset-pelka-scores.tsv"), all.scores.pelka, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(file = paste0(figs.dir, "/cancer-validation-dataset-wu-pelka-merged-scores.tsv"), all.scores.cancer.merged, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

titles <- list("wu" = "Wu (BRCA)", "pelka" = "Pelka (CRC)", "all" = "Merged Wu (BRCA) and Pelka (CRC)")


results <- list("wu" = all.scores.wu, 
                "pelka" = all.scores.pelka, 
                "all" = all.scores.cancer.merged)

nms <- names(titles)
names(nms) <- nms

gs <- llply(nms,
            .fun = function(nm) {
              cell.type.levels <- calculate.cell.type.levels(results[[nm]], id.var = "method.name", cell.type.var = "cell.type", cor.var = "cor.p")
              method.levels <- calculate.method.levels(results[[nm]], id.var = "method.name", cell.type.var = "cell.type", cor.var = "cor.p")
              g <- plot.cell.type.correlation.heatmap(results[[nm]], show.corr.text = TRUE,
                                                      id.var = "method.name", cell.type.var = "cell.type", cor.var = "cor.p",
                                                      ids.to.bold = comparator.methods,
                                                      second.col.summary.fun = "mean",
                                                      cell.type.levels = rev(cell.type.levels),
                                                      method.levels = method.levels)
              g <- g + theme(plot.title = element_text(hjust = 0.5))
              g <- g + ggtitle(titles[[nm]])
              png(paste0(figs.dir, "fig-cancer-validation-heatmap-", nm, ".png"), width = 2 * 480, height = 1 * 480)                    
              print(g)
              d <- dev.off()
              
              pdf(paste0(figs.dir, "fig-cancer-validation-heatmap-", nm, ".pdf"), width = 2 * 7, height = 1 * 7)                    
              print(g)
              d <- dev.off()
              
              g
            })

g <- plot_grid(gs[[1]], gs[[2]], nrow=2, labels="AUTO")
png(paste0(figs.dir, "fig-cancer-validation-heatmap-wu-and-pelka.png"), width = 2 * 480, height = 2 * 480)                    
print(g)
d <- dev.off()

pdf(paste0(figs.dir, "fig-cancer-validation-heatmap-wu-and-pelka.pdf"), width = 2 * 7, height = 2 * 7)                    
print(g)
d <- dev.off()

cat("Exiting successfully")
q(status=0)

validation.upset <- validation.metadata[, c("dataset"), drop=F]
cols <- list("tumor.type" = "BRCA", "distribution.type" = "Biological", "mixture.type" = "In Vitro", "immune.origin" = "Healthy")
for(col in names(cols)) {
  validation.upset[cols[[col]]] <- FALSE
  flag <- validation.metadata[, col] == cols[[col]]
  validation.upset[flag, cols[[col]]] <- TRUE
}
rownames(validation.upset) <- validation.upset$dataset
validation.upset <- validation.upset[, -1]

validation.tidy <- validation.upset %>%
  as_tibble(rownames = "dataset") %>%
  gather(Label, Member, -dataset) %>%
  filter(Member) %>%
  select(- Member)

validation.tidy <- validation.tidy %>%
  group_by(dataset) %>%
  summarize(Label = list(Label))



g <- ggplot(data = subset(all.scores, cell.type == "CD4.T.cells"))
g <- g + geom_boxplot(aes(x = dataset.name, y = cor.p))
g <- g + geom_point(aes(x = dataset.name, y = cor.p))

g <- ggplot(data = all.scores)
g <- g + geom_boxplot(aes(x = dataset.name, y = cor.p))
g <- g + geom_point(aes(x = dataset.name, y = cor.p))
g <- g + facet_wrap(~ cell.type)


"CIBERSORT"
"EPIC"
"MCP-counter"
"quanTIseq"
"TIMER"
"xCell" 

"CIBERSORTx" 



merge.predictions.and.ground.truth <- function(preds, ground.truth) {
  preds <- merge(preds, ground.truth, all=TRUE)
  preds <- preds[!is.na(preds$measured),]
}


submitterIds <- as.character(unique(na.omit(res.all$submitterId)))
names(submitterIds) <- submitterIds
names <- unlist(lapply(submitterIds, translate.submitterId))
df <- data.frame(team = as.character(names), submitterId = as.numeric(names(names)))
df <- merge(df, unique(res.all[, c("submitterId", "method.name")]), all.x=TRUE)



coarse.gt <- read.table(synGet(coarse.gt.synId)$path, header=TRUE, sep=",")
fine.gt <- read.table(synGet(fine.gt.synId)$path, header=TRUE, sep=",")

datasets <- c("AA", "AB", "AE", "AF", "DS1", "DS2", "DS3", "DS4")
aginome.coarse.preds <- merge.predictions.across.datasets("aginome", datasets, "coarse", coarse.gt)
aginome.fine.preds <- merge.predictions.across.datasets("aginome", datasets, "fine", fine.gt)

aginome.coarse.scores <- score.preds(aginome.coarse.preds)
aginome.fine.scores <- score.preds(aginome.fine.preds)


aginome.coarse.orig.preds <- subset(coarse.orig.preds, method.name == "Aginome-XMU")
aginome.coarse.orig.preds <- merge.predictions.and.ground.truth(aginome.coarse.orig.preds, coarse.gt)
aginome.coarse.orig.scores <- score.preds(aginome.coarse.orig.preds)

aginome.fine.orig.preds <- subset(fine.orig.preds, method.name == "Aginome-XMU")
aginome.fine.orig.preds <- merge.predictions.and.ground.truth(aginome.fine.orig.preds, fine.gt)
aginome.fine.orig.scores <- score.preds(aginome.fine.orig.preds

stop("stop")


# Get the input file for the challenge.
# Note that this is the input using the same admixtures for both fine- and coarse-grained,
# i.e., the 'input-translated-from-fine.csv'
# Originally, we had different admixtures for fine- and for coarse-grained.
input.synId <- "syn22267272"

# Read in the input file
input.obj <- synGet(input.synId, downloadFile=TRUE)
input.tbl <- read.table(input.obj$path, header=TRUE, sep=",")

# Get the folder of the input file, which will hold the data files listed within it.
data.folder.synId <- input.obj$properties$parentId
children <- synGetChildren(data.folder.synId)
l <- as.list(children)
df <- do.call(rbind.data.frame, l)

# Download the data files of interest (e.g., symbol or ensembl id-based, TPM or count-based)
# "hugo.expr.file"                        
# "hugo.expr.est.counts.file"             
# "ensg.expr.file"                        
# "ensg.expr.est.counts.file"       
files <- input.tbl[, c("dataset.name", "hugo.expr.file")]
colnames(files) <- c("dataset.name", "file")

files <- merge(files, df[, c("name", "id")], by.x = c("file"), by.y = c("name"))
for(col in colnames(files)) { files[, col] <- as.character(files[, col]) }

for(i in 1:nrow(files)) {

  # Download the data file
  mat <- read.table(synGet(files[i, "id"], downloadFile=TRUE)$path, sep=",", header=TRUE)

  # Make the Gene column the rowname
  mat$Gene <- as.character(mat$Gene)
  rownames(mat) <- mat$Gene
  mat <- mat[, !(colnames(mat) == "Gene")]
  
  # Save the file
  ofile <- paste0(files[i, "dataset.name"], "_symbol_tpm.csv")
  write.table(file=ofile, mat, row.names=TRUE, col.names=TRUE, quote=FALSE)
}

stop("stop")

datasets <- names(files)
names(datasets) <- datasets

results <- list()
results[["coarse"]] <- list()
results[["coarse"]][["xmu"]] <- "xmu-coarse-challenge-admixtures-prediction.csv"
results[["fine"]] <- list()
results[["fine"]][["xmu"]] <- "xmu-fine-challenge-admixtures-prediction.csv"

