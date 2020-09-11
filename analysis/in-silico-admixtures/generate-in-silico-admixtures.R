usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, repos = "http://cran.us.r-project.org", dep = TRUE)
  require(p, character.only = TRUE)
}

usePackage("pacman")

home_dir <- "../"
cur_dir <- getwd()
cache_dir <- paste0(home_dir, "/admixtures/new-admixtures/")

setwd(cache_dir)
source("utils.R")
source("load-data.R")
setwd(cur_dir)
source("../define-biological-models.R")

## Spike-in dilution experiment. 

cat(paste0("This script (may) require a lot of memory\n"))
cat(paste0("I've run in successfully on the bassh2 node\n"))
cat(paste0("i.e., bwhite-beat-aml-bioc-again\n"))
cat(paste0("Private IP: 10.23.19.117\n"))
cat(paste0("Instance type r3.4xlarge\n"))
cat(paste0("ami-08333c98cd3b5ada7\n"))

suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(reshape2))
## install.packages("dirichlet", repos="http://R-Forge.R-project.org")
suppressPackageStartupMessages(p_load(dirichlet))
suppressPackageStartupMessages(p_load(hitandrun))

suppressPackageStartupMessages(p_load("foreach"))
suppressPackageStartupMessages(p_load("parallel"))

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(p_load("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

synLogin()

## Load the TPM validation data
##synId <- "syn21574299"
##obj <- synGet(synId, downloadFile = TRUE)
##cpm.expr <- read.table(obj$path, sep = "\t", header = TRUE)

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

coarse.grained.map <- list(
    "Breast" = "Breast",
    "CRC" = "CRC",
    "Dendritic_cells" = "monocytic.lineage",
    "Endothelial_cells" = "endothelial.cells",
    "Fibroblasts" = "fibroblasts",
    "Macrophages" = "monocytic.lineage",
    "Memory_CD4_T_cells" = "CD4.T.cells",
    "Memory_CD8_T_cells" = "CD8.T.cells",
    "Monocytes" = "monocytic.lineage",
    "NK_cells" = "NK.cells",
    "Naive_B_cells" = "B.cells",
    "Naive_CD4_T_cells" = "CD4.T.cells",
    "Naive_CD8_T_cells" = "CD8.T.cells",
    "Neutrophils" = "neutrophils",
    "Tregs" = "CD4.T.cells")
    

fine.grained.map.df <- data.frame(population = names(fine.grained.map), challenge.population = unname(unlist(fine.grained.map)))
fine.grained.pop.df <- merge(pop.df, fine.grained.map.df)

coarse.grained.map.df <- data.frame(population = names(coarse.grained.map), challenge.population = unname(unlist(coarse.grained.map)))
coarse.grained.pop.df <- merge(pop.df, coarse.grained.map.df)

## Separate the samples into two batches
fine.grained.pop.df <- fine.grained.pop.df[order(fine.grained.pop.df$sample),]
fine.grained.pop1.df <- ddply(fine.grained.pop.df, .variables = c("population"), .fun = function(df) df[1,,drop=F])
fine.grained.pop2.df <- ddply(fine.grained.pop.df, .variables = c("population"), .fun = function(df) df[min(2,nrow(df)),,drop=F])

coarse.grained.pop.df <- coarse.grained.pop.df[order(coarse.grained.pop.df$sample),]
coarse.grained.pop1.df <- ddply(coarse.grained.pop.df, .variables = c("population"), .fun = function(df) df[1,,drop=F])
coarse.grained.pop2.df <- ddply(coarse.grained.pop.df, .variables = c("population"), .fun = function(df) df[min(2,nrow(df)),,drop=F])

immune.fine.grained.populations <- unique(as.character(fine.grained.pop.df$challenge.population))
immune.fine.grained.populations <- immune.fine.grained.populations[immune.fine.grained.populations != "CRC"]
immune.fine.grained.populations <- immune.fine.grained.populations[immune.fine.grained.populations != "Breast"]
num.immune.fine.grained.populations <- length(immune.fine.grained.populations)

immune.coarse.grained.populations <- unique(as.character(coarse.grained.pop.df$challenge.population))
immune.coarse.grained.populations <- immune.coarse.grained.populations[immune.coarse.grained.populations != "CRC"]
immune.coarse.grained.populations <- immune.coarse.grained.populations[immune.coarse.grained.populations != "Breast"]
num.immune.coarse.grained.populations <- length(immune.coarse.grained.populations)

create.rand.admixture.wo.cancer <- function(populations, spike.ins) {

    num.pops <- length(populations)
    admixs <-
        ldply(spike.ins,
              .fun = function(spike.in) {
                  dummy <- 1:(num.pops-1)
                  df <- generate.random.uniform.admixtures(dummy, n = 5, min.prop = 0, max.prop = 1 - spike.in, step = 0.001)
                  t(df)
              })
    admixs[,1] <- as.numeric(admixs[,1])

    ret <- ldply(1:num.pops,
                 .fun = function(i) {
                     tmp <- admixs
                     colnames(tmp) <- populations
                     tmp.name <- colnames(tmp)[1]
                     spike.in.pop <- colnames(tmp)[i]
                     colnames(tmp)[1] <- spike.in.pop
                     colnames(tmp)[i] <- tmp.name
                     tmp[, populations]
                 })
    spike.in.df <-
        ldply(1:num.pops,
                 .fun = function(i) {
                     tmp <- admixs
                     colnames(tmp) <- populations
                     data.frame(spike.in.pop = colnames(tmp)[i], spike.in.prop = tmp[,1])
                 })
    return(list("ret" = ret, "spike.in.df" = spike.in.df))
}

create.rand.admixture.w.cancer <- function(non.cancer.populations, cancer, spike.ins) {

    num.pops <- length(non.cancer.populations)
    admixs <-
        ldply(spike.ins,
              .fun = function(spike.in) {
                  dummy <- 1:num.pops
                  df <- generate.random.uniform.admixtures(dummy, n = 15, min.prop = 0, max.prop = 1 - spike.in, step = 0.001)
                  t(df)
              })
    admixs[,1] <- as.numeric(admixs[,1])

    colnames(admixs) <- c(cancer, non.cancer.populations)

    spike.in.df <-
        data.frame(spike.in.pop = cancer, spike.in.prop = admixs[,1])
    
    return(list("admixs" = admixs, "spike.in.df" = spike.in.df))
}

generate.random.admixtures.wo.cancer <- function(immune.populations, population.maps, spike.ins) {

    rand.admixtures <-
        ldply(population.maps,
              .fun = function(pop.df) {
                  rt <- create.rand.admixture.wo.cancer(immune.populations, spike.ins)
                  ret <- rt[["ret"]]
                  spike.in.df <- rt[["spike.in.df"]]
                  spike.in.df$admixture <- 1:nrow(spike.in.df)
                  long <- melt(as.matrix(ret))
                  colnames(long) <- c("admixture", "challenge.population", "prop")
                  long <- merge(long, spike.in.df)
                  long$challenge.population <- as.character(long$challenge.population)

                  long$spike.in.pop <- as.character(long$spike.in.pop)
                  flag <- long$spike.in.pop == long$challenge.population
                  long$spike.in.prop <- NA
                  long[flag, "spike.in.prop"] <- long[flag, "prop"]

                  pop.df$challenge.population <- as.character(pop.df$challenge.population)
                  col <- "sample"
                  long <- merge(long, pop.df[, c("challenge.population", col)], all.x = TRUE)
                  long <-
                      ddply(long, .variables = c("challenge.population", "admixture"),
                            .fun = function(df) {
                                if(nrow(df) == 1) { return(df) }
                                ## Distribute the proportion across the samples mapped to this
                                ## challenge population
                                ## (e.g., in the coarse-grained challenge memory CD4, naive CD4,
                                ## and Tregs map to CD4)
                                props <- as.numeric(rdirichlet(n = 1, alpha = rep(1/nrow(df), nrow(df))))
                                df$prop <- df$prop * props
                                df
                            })
                  ## long <- long[, c("admixture", col, "prop")]
                  long <- long[, c("admixture", col, "prop", "spike.in.pop", "spike.in.prop")]
                  long[,col] <- as.character(long[,col])
                  long[,"prop"] <- as.numeric(long[,"prop"])
                  long[,"spike.in.prop"] <- as.numeric(long[,"spike.in.prop"])                  
                  long[,"spike.in.pop"] <- as.character(long[,"spike.in.pop"])                  
                  long
              })
    
    rand.admixtures$admixture <- paste0(rand.admixtures$.id, rand.admixtures$admixture)
    rand.admixtures <- rand.admixtures[, !(colnames(rand.admixtures) == ".id")]

    ## Confirm all admixtures sum to one
    sums <-
        ddply(rand.admixtures, .variables = c("admixture"),
              .fun = function(df) sum(df$prop))
    eps <- 10^-5
    flag <- abs(sums$V1 - 1) < eps 
    if(!(all(flag))) {
        print(head(sums[!flag,]))
        stop("Some admixtures don't sum to one\n")
    }

    rand.admixtures
}

generate.random.admixtures.w.cancer <- function(immune.populations, population.maps, spike.ins, cancers = c("Breast", "CRC")) {

    names(cancers) <- c("C", "D")
    rand.admixtures <-
        ldply(cancers,
              .fun = function(cancer) {
                  tmp <- ldply(population.maps,
                        .fun = function(pop.df) {
                            rt <- create.rand.admixture.w.cancer(immune.populations, cancer, spike.ins)
                            ret <- rt[["admixs"]]
                            spike.in.df <- rt[["spike.in.df"]]
                            spike.in.df$admixture <- 1:nrow(spike.in.df)
                            long <- melt(as.matrix(ret))
                            colnames(long) <- c("admixture", "challenge.population", "prop")
                            long <- merge(long, spike.in.df)
                            long$challenge.population <- as.character(long$challenge.population)
                            long$spike.in.pop <- as.character(long$spike.in.pop)

                            flag <- long$spike.in.pop == long$challenge.population
                            long$spike.in.prop <- NA
                            long[flag, "spike.in.prop"] <- long[flag, "prop"]
                            
                            pop.df$challenge.population <- as.character(pop.df$challenge.population)
                            col <- "sample"
                            long <- merge(long, pop.df[, c("challenge.population", col)], all.x = TRUE)
                            long <-
                                ddply(long, .variables = c("challenge.population", "admixture"),
                                      .fun = function(df) {
                                          if(nrow(df) == 1) { return(df) }
                                          ## Distribute the proportion across the samples mapped to this
                                          ## challenge population
                                          ## (e.g., in the coarse-grained challenge memory CD4, naive CD4,
                                          ## and Tregs map to CD4)
                                          props <- as.numeric(rdirichlet(n = 1, alpha = rep(1/nrow(df), nrow(df))))
                                          df$prop <- df$prop * props
                                          df
                                      })
                            ## long <- long[, c("admixture", col, "prop")]
                            long <- long[, c("admixture", col, "prop", "spike.in.pop", "spike.in.prop")]
                            long[,col] <- as.character(long[,col])
                            long[,"prop"] <- as.numeric(long[,"prop"])
                            long[,"spike.in.prop"] <- as.numeric(long[,"spike.in.prop"])                  
                            long[,"spike.in.pop"] <- as.character(long[,"spike.in.pop"])                  
                            long
                        })
                  colnames(tmp)[1] <- "pop.id"
                  tmp
              })
    colnames(rand.admixtures)[1] <- "cancer.id"
    
    rand.admixtures$admixture <- paste0(rand.admixtures$pop.id, rand.admixtures$cancer.id, rand.admixtures$admixture)
    rand.admixtures <- rand.admixtures[, !(colnames(rand.admixtures) %in% c("cancer.id", "pop.id"))]

    ## Confirm all admixtures sum to one
    sums <-
        ddply(rand.admixtures, .variables = c("admixture"),
              .fun = function(df) sum(df$prop))
    eps <- 10^-5
    flag <- abs(sums$V1 - 1) < eps 
    if(!(all(flag))) {
        print(head(sums[!flag,]))
        stop("Some admixtures don't sum to one\n")
    }

    rand.admixtures
}

fine.grained.pop.maps <- list("A" = fine.grained.pop1.df, "B" = fine.grained.pop2.df)
coarse.grained.pop.maps <- list("A" = coarse.grained.pop1.df, "B" = coarse.grained.pop2.df)

spike.ins <- c(1/(2^seq(from=2,to=12,by=1)), 0)
spike.ins <- unique(c(0, seq(from=0,to=0.1,by=0.01), seq(from=0.1,to=1,by=0.1), seq(from=1,to=20,by=1), seq(from=20,to=40,by=2))) / 100
spike.ins <- sort(spike.ins, decreasing=TRUE)
print(spike.ins)
names(spike.ins) <- spike.ins

cancer.spike.ins <- c(3/4, 1/2, 1/4, 1/8, 1/16, 1/32, 0)
names(cancer.spike.ins) <- cancer.spike.ins

set.seed(100)

coarse.grained.rand.admixtures.wo.cancer <-
    generate.random.admixtures.wo.cancer(immune.coarse.grained.populations, coarse.grained.pop.maps, spike.ins)

fine.grained.rand.admixtures.wo.cancer <-
    generate.random.admixtures.wo.cancer(immune.fine.grained.populations, fine.grained.pop.maps, spike.ins)

num.w.cancer.datasets <- 5

fine.grained.rand.admixtures.w.cancers <-
    llply(1:num.w.cancer.datasets,
          function(i) 
              generate.random.admixtures.w.cancer(immune.fine.grained.populations, fine.grained.pop.maps, cancer.spike.ins,
                                                  cancers = c("CRC", "Breast")))

coarse.grained.rand.admixtures.w.cancers <-
    llply(1:num.w.cancer.datasets,
          function(i) 
              generate.random.admixtures.w.cancer(immune.coarse.grained.populations, coarse.grained.pop.maps, cancer.spike.ins,
                                                  cancer = c("CRC", "Breast")))

flat.model <- define.biological.model(unname(unlist(populations)), cache.dir = cache_dir)

generate.biological.admixtures_ <- function(model, spike.in.population, other.populations, spike.ins, n = 20) {
    model <- model[, other.populations]
    ndim <- ncol(model)
    n.tot <- 1000
    ##N <- 1000*ndim^3
    ## hitandrun documentation says that, since har uses MCMC which is time-correlated,
    ## need to run for O(n^3) steps.
    prefactor <- 100
    if(prefactor < ndim) { stop("Expected prefator >= ndim to ensure O(n^3) steps\n") }
    N <- prefactor*(ndim)^3
    N.thin <- N
    n.rand <- 5
    N <- 10*ndim^3
    N.thin <- 10*ndim^3
    ## despite what documentation says, output has length N / N.thin (not N * N.thin)
    N <- 1000*(ndim)^3
    N.thin <- (ndim)^3
    names(spike.ins) <- spike.ins
    admixtures <-
        ldply(spike.ins, .parallel = TRUE,
              .fun = function(spike.in) {

                  eq.constr <- list(
                      constr = matrix(rep(1,ndim),nrow=1,ncol=(ndim)),
                      dir = '=',
                      rhs = 1-spike.in)

                  ineq.constr <- list(
                      constr = rbind(-diag(ndim),diag(ndim)),
                      dir = rep('<=',(ndim)*2),
                      rhs = c(-as.numeric(model["min",]),as.numeric(model["max",])))

                  basis <- solution.basis(eq.constr)
                  transform <- createTransform(basis)
                  constr <- transformConstraints(transform,ineq.constr)

                  x0 <- createSeedPoint(constr,homogeneous=TRUE,randomize=TRUE)
                  df <- har(x0,constr,N,thin=N.thin,transform=transform, homogeneous=TRUE)$samples
                  eps <- 10^-4
                  flag <- abs( rowSums(df) - (1 - spike.in) ) < eps
                  df <- df[flag, ,drop=F]
                  df <- df[!duplicated(df),,drop=F]
                  colnames(df) <- colnames(model)
                  df
              })
    colnames(admixtures)[1] <- spike.in.population
    admixtures[, spike.in.population] <- as.numeric(admixtures[, spike.in.population])
    ret <-
        dlply(admixtures, .variables = spike.in.population,
              .fun = function(df) {
                  find.extremal.points(t(df), n)
              })
    t(do.call(cbind, ret))
}

generate.biological.admixtures <- function(model, spike.in.population, other.populations, population.maps, spike.ins, n = 10) {

    admixtures <-
        ldply(population.maps,
              .fun = function(pop.df) {
                  ret <- generate.biological.admixtures_(model, spike.in.population, other.populations, spike.ins, n = n)
                  long <- melt(as.matrix(ret))
                  colnames(long) <- c("admixture", "challenge.population", "prop")
                  long$challenge.population <- as.character(long$challenge.population)
                  long$spike.in.pop <- as.character(spike.in.population)
                  flag <- long$spike.in.pop == long$challenge.population
                  long[flag, "spike.in.prop"] <- long[flag, "prop"]
                  pop.df$challenge.population <- as.character(pop.df$challenge.population)
                  col <- "sample"
                  long <- merge(long, pop.df[, c("challenge.population", col)], all.x = TRUE)
                  long <-
                      ddply(long, .variables = c("challenge.population", "admixture"),
                            .fun = function(df) {
                                if(nrow(df) == 1) { return(df) }
                                ## Distribute the proportion across the samples mapped to this
                                ## challenge population
                                ## (e.g., in the coarse-grained challenge memory CD4, naive CD4,
                                ## and Tregs map to CD4)
                                props <- as.numeric(rdirichlet(n = 1, alpha = rep(1/nrow(df), nrow(df))))
                                old <- df$prop[1]
                                if(!all(df$prop == df$prop[1])) {
                                    stop(paste0("Not all props are the same for pop ", df$challenge.population[1], "\n"))
                                }
                                eps <- 10^-4
                                df$prop <- df$prop * props
                                if(abs(sum(df$prop) - old) > eps) {
                                    print(props)
                                    stop(paste0(sum(df$prop), " != ", old, "\n"))
                                }
                                df
                            })
                  ## long <- long[, c("admixture", col, "prop")]
                  long <- long[, c("admixture", col, "prop", "spike.in.pop", "spike.in.prop")]
                  long[,col] <- as.character(long[,col])
                  long[,"prop"] <- as.numeric(long[,"prop"])
                  long[,"spike.in.prop"] <- as.numeric(long[,"spike.in.prop"])                  
                  long[,"spike.in.pop"] <- as.character(long[,"spike.in.pop"])                  
                  long
              })
    colnames(admixtures)[1] <- "pop.id"
    
    admixtures$admixture <- paste0(admixtures$pop.id, admixtures$admixture)
    admixtures <- admixtures[, !(colnames(admixtures) %in% c("pop.id"))]

    ## Confirm all admixtures sum to one
    sums <-
        ddply(admixtures, .variables = c("admixture"),
              .fun = function(df) sum(df$prop))
    eps <- 10^-5
    flag <- abs(sums$V1 - 1) < eps 
    if(!(all(flag))) {
        print(head(sums[!flag,]))
        stop("Some admixtures don't sum to one\n")
    }

    admixtures
}

flat.model.no.cancer <- flat.model[, !(colnames(flat.model) %in% "tumor.fraction")]
tmp <- subset(fine.grained.map.df, population %in% colnames(flat.model.no.cancer))
flat.model.no.cancer <- flat.model.no.cancer[, as.character(tmp$population)]
colnames(flat.model.no.cancer) <- as.character(tmp$challenge.population)

flat.model.no.cancer.coarse <- flat.model.no.cancer
flat.model.no.cancer.coarse <- cbind(flat.model.no.cancer.coarse, "monocytic.lineage" =
    flat.model.no.cancer.coarse[, "myeloid.dendritic.cells"] +
    flat.model.no.cancer.coarse[, "macrophages"] +
    flat.model.no.cancer.coarse[, "monocytes"])
flat.model.no.cancer.coarse <- cbind(flat.model.no.cancer.coarse, "CD4.T.cells" =
    flat.model.no.cancer.coarse[, "memory.CD4.T.cells"] +
    flat.model.no.cancer.coarse[, "naive.CD4.T.cells"] +
    flat.model.no.cancer.coarse[, "regulatory.T.cells"])
flat.model.no.cancer.coarse <- cbind(flat.model.no.cancer.coarse, "CD8.T.cells" =
    flat.model.no.cancer.coarse[, "memory.CD8.T.cells"] +
    flat.model.no.cancer.coarse[, "naive.CD8.T.cells"])
flat.model.no.cancer.coarse <- cbind(flat.model.no.cancer.coarse, "B.cells" =
                                                                      flat.model.no.cancer.coarse[, "naive.B.cells"])
## all of the population minimums were set to 0.01, so when we add those we get something higher.
## reset to 0.1
flat.model.no.cancer.coarse["min", "monocytic.lineage"] <- 0.01
flat.model.no.cancer.coarse["min", "CD4.T.cells"] <- 0.01
flat.model.no.cancer.coarse["min", "CD8.T.cells"] <- 0.01

tmp <- subset(coarse.grained.map.df, challenge.population %in% colnames(flat.model.no.cancer.coarse))
flat.model.no.cancer.coarse <- flat.model.no.cancer.coarse[, unique(as.character(tmp$challenge.population))]

## Generate fine-grained biological admixtures w/o cancer
## Generate 5 replicates for each of the two sets of samples
spike.in.populations <- colnames(flat.model.no.cancer)
names(spike.in.populations) <- LETTERS[1:length(spike.in.populations)]

my.letters <- c(LETTERS, paste0("A", LETTERS))

cat("Generating fine-grained biological admixtures w/o cancer\n")
fine.grained.bio.admixtures.wo.cancer <-
    ldply(spike.in.populations,
          .parallel = TRUE,
          .fun = function(spike.in.population) {
              other.populations <-
                  colnames(flat.model.no.cancer)[!(colnames(flat.model.no.cancer) %in% c(spike.in.population, "tumor.fraction"))]
              ret <- generate.biological.admixtures(flat.model.no.cancer, spike.in.population, other.populations,
                                                    fine.grained.pop.maps, spike.ins, n = 5)
          })

if(!(".id" %in% colnames(fine.grained.bio.admixtures.wo.cancer))) {
    stop(paste0(".id not in fine.grained.bio.admixtures.wo.cancer\n"))
}
fine.grained.bio.admixtures.wo.cancer$admixture <- paste0(fine.grained.bio.admixtures.wo.cancer$.id, fine.grained.bio.admixtures.wo.cancer$admixture)
fine.grained.bio.admixtures.wo.cancer <- fine.grained.bio.admixtures.wo.cancer[, !(colnames(fine.grained.bio.admixtures.wo.cancer) %in% ".id")]

spike.in.populations <- colnames(flat.model.no.cancer.coarse)
names(spike.in.populations) <- LETTERS[1:length(spike.in.populations)]

## Generate coarse-grained biological admixtures w/o cancer
cat("Generating coarse-grained biological admixtures w/o cancer\n")
coarse.grained.bio.admixtures.wo.cancer <-
    ldply(spike.in.populations,
          .parallel = TRUE,
          .fun = function(spike.in.population) {
              other.populations <-
                  colnames(flat.model.no.cancer.coarse)[!(colnames(flat.model.no.cancer.coarse) %in% c(spike.in.population, "tumor.fraction"))]
              ret <- generate.biological.admixtures(flat.model.no.cancer.coarse, spike.in.population, other.populations,
                                                    coarse.grained.pop.maps, spike.ins, n = 5)
          })
if(!(".id" %in% colnames(coarse.grained.bio.admixtures.wo.cancer))) {
    stop(paste0(".id not in coarse.grained.bio.admixtures.wo.cancer\n"))
}
coarse.grained.bio.admixtures.wo.cancer$admixture <- paste0(coarse.grained.bio.admixtures.wo.cancer$.id, coarse.grained.bio.admixtures.wo.cancer$admixture)
coarse.grained.bio.admixtures.wo.cancer <- coarse.grained.bio.admixtures.wo.cancer[, !(colnames(coarse.grained.bio.admixtures.wo.cancer) %in% ".id")]


## Generate fine-grained biological admixtures w/ cancer
## For each cancer, generate 10 replicates for each of two sets of samples
spike.in.populations <- list("C" = "CRC", "D" = "Breast")

cat("Generating fine-grained biological admixtures w/ cancer\n")
fine.grained.bio.admixtures.w.cancers <-
    llply(1:num.w.cancer.datasets,
          function(i) {
              fine.grained.bio.admixtures.w.cancer <-
                  ldply(spike.in.populations,
                        .parallel = TRUE,
                        .fun = function(spike.in.population) {
                            other.populations <-
                                colnames(flat.model.no.cancer)[!(colnames(flat.model.no.cancer) %in% c(spike.in.population, "tumor.fraction"))]
                            ret <- generate.biological.admixtures(flat.model.no.cancer, spike.in.population, other.populations,
                                                                  fine.grained.pop.maps, cancer.spike.ins, n = 10)
                        })
              
              if(!(".id" %in% colnames(fine.grained.bio.admixtures.w.cancer))) {
                  stop(paste0(".id not in fine.grained.bio.admixtures.w.cancer\n"))
              }
              fine.grained.bio.admixtures.w.cancer$admixture <- paste0(fine.grained.bio.admixtures.w.cancer$.id, fine.grained.bio.admixtures.w.cancer$admixture)
              fine.grained.bio.admixtures.w.cancer <- fine.grained.bio.admixtures.w.cancer[, !(colnames(fine.grained.bio.admixtures.w.cancer) %in% ".id")]
              fine.grained.bio.admixtures.w.cancer
          })

## Generate coarse-grained biological admixtures w/ cancer
cat("Generating coarse-grained biological admixtures w/ cancer\n")
coarse.grained.bio.admixtures.w.cancers <-
    llply(1:num.w.cancer.datasets,
          function(i) {
              coarse.grained.bio.admixtures.w.cancer <-
                  ldply(spike.in.populations,
                        .parallel = TRUE,
                        .fun = function(spike.in.population) {
                            other.populations <-
                                colnames(flat.model.no.cancer.coarse)[!(colnames(flat.model.no.cancer.coarse) %in% c(spike.in.population, "tumor.fraction"))]
                            ret <- generate.biological.admixtures(flat.model.no.cancer.coarse, spike.in.population, other.populations,
                                                                  coarse.grained.pop.maps, cancer.spike.ins, n = 10)
                        })
              if(!(".id" %in% colnames(coarse.grained.bio.admixtures.w.cancer))) {
                  stop(paste0(".id not in coarse.grained.bio.admixtures.w.cancer\n"))
              }
              coarse.grained.bio.admixtures.w.cancer$admixture <- paste0(coarse.grained.bio.admixtures.w.cancer$.id, coarse.grained.bio.admixtures.w.cancer$admixture)
              coarse.grained.bio.admixtures.w.cancer <- coarse.grained.bio.admixtures.w.cancer[, !(colnames(coarse.grained.bio.admixtures.w.cancer) %in% ".id")]
              coarse.grained.bio.admixtures.w.cancer
          })

nxt.letter <- 1
datasets <- list()

datasets[[my.letters[nxt.letter]]] <-
    list("data" = fine.grained.rand.admixtures.wo.cancer, "mixture.type" = "Random", "subchallenge" = "fine", "tumor.type" = NA)
nxt.letter <- nxt.letter + 1

datasets[[my.letters[nxt.letter]]] <-
    list("data" = fine.grained.bio.admixtures.wo.cancer, "mixture.type" = "Biological", "subchallenge" = "fine", "tumor.type" = NA)
nxt.letter <- nxt.letter + 1    

datasets[[my.letters[nxt.letter]]] <-
    list("data" = coarse.grained.rand.admixtures.wo.cancer, "mixture.type" = "Random", "subchallenge" = "coarse", "tumor.type" = NA)
nxt.letter <- nxt.letter + 1

datasets[[my.letters[nxt.letter]]] <-
    list("data" = coarse.grained.bio.admixtures.wo.cancer, "mixture.type" = "Biological", "subchallenge" = "coarse", "tumor.type" = NA)
nxt.letter <- nxt.letter + 1    

if(FALSE) {
l <- fine.grained.rand.admixtures.w.cancers
for(ds in l) {
    flag <- ds$spike.in.pop == "CRC"
    crc.ds <- ds[flag,]
    flag <- ds$spike.in.pop == "Breast"    
    breast.ds <- ds[flag,]
    datasets[[my.letters[nxt.letter]]] <- list("data" = crc.ds, "mixture.type" = "Random", "subchallenge" = "fine", "tumor.type" = "CRC")
    nxt.letter <- nxt.letter + 1
    datasets[[my.letters[nxt.letter]]] <- list("data" = breast.ds, "mixture.type" = "Random", "subchallenge" = "fine", "tumor.type" = "BRCA")
    nxt.letter <- nxt.letter + 1
}

l <- fine.grained.bio.admixtures.w.cancers
for(ds in l) {
    flag <- ds$spike.in.pop == "CRC"
    crc.ds <- ds[flag,]
    flag <- ds$spike.in.pop == "Breast"    
    breast.ds <- ds[flag,]
    datasets[[my.letters[nxt.letter]]] <- list("data" = crc.ds, "mixture.type" = "Biological", "subchallenge" = "fine", "tumor.type" = "CRC")
    nxt.letter <- nxt.letter + 1
    datasets[[my.letters[nxt.letter]]] <- list("data" = breast.ds, "mixture.type" = "Biological", "subchallenge" = "fine", "tumor.type" = "BRCA")
    nxt.letter <- nxt.letter + 1
}

l <- coarse.grained.rand.admixtures.w.cancers
for(ds in l) {
    flag <- ds$spike.in.pop == "CRC"
    crc.ds <- ds[flag,]
    flag <- ds$spike.in.pop == "Breast"    
    breast.ds <- ds[flag,]
    datasets[[my.letters[nxt.letter]]] <- list("data" = crc.ds, "mixture.type" = "Random", "subchallenge" = "coarse", "tumor.type" = "CRC")
    nxt.letter <- nxt.letter + 1
    datasets[[my.letters[nxt.letter]]] <- list("data" = breast.ds, "mixture.type" = "Random", "subchallenge" = "coarse", "tumor.type" = "BRCA")
    nxt.letter <- nxt.letter + 1
}

l <- coarse.grained.bio.admixtures.w.cancers
for(ds in l) {
    flag <- ds$spike.in.pop == "CRC"
    crc.ds <- ds[flag,]
    flag <- ds$spike.in.pop == "Breast"    
    breast.ds <- ds[flag,]
    datasets[[my.letters[nxt.letter]]] <- list("data" = crc.ds, "mixture.type" = "Biological", "subchallenge" = "coarse", "tumor.type" = "CRC")
    nxt.letter <- nxt.letter + 1
    datasets[[my.letters[nxt.letter]]] <- list("data" = breast.ds, "mixture.type" = "Biological", "subchallenge" = "coarse", "tumor.type" = "BRCA")
    nxt.letter <- nxt.letter + 1
}
}

metadata <- ldply(datasets, .fun = function(entry) as.data.frame(entry[names(entry) != "data"], stringsAsFactors = FALSE))
colnames(metadata)[1] <- "dataset.name"
coarse.datasets <- subset(metadata, subchallenge == "coarse")$dataset.name
fine.datasets <- subset(metadata, subchallenge == "fine")$dataset.name

metadata.file <- "in-silico-metadata.csv"
write.table(file = metadata.file, metadata, row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)

## Create the gold standards
## gold standard files should be csvs with columns: dataset.name,sample.id,cell.type,measured

## NB: these need to be sum'ed across cell populations that map to the same challenge subtype
## NB: ideally, include both in the ground truth

coarse.gs <- ldply(datasets, .fun = function(df) df$data)
colnames(coarse.gs) <- c("dataset.name", "sample.id", "sample", "measured", "spike.in.pop", "spike.in.prop")
coarse.gs <- subset(coarse.gs, dataset.name %in% coarse.datasets)
old.sz <- nrow(coarse.gs)
coarse.gs <- merge(coarse.gs, coarse.grained.pop.df[, c("sample", "challenge.population")], by = "sample")
new.sz <- nrow(coarse.gs)
if(old.sz != new.sz) { stop(paste0("Size changed from ", old.sz, " to ", new.sz, "\n")) }
coarse.gs <- coarse.gs[, c("dataset.name", "sample.id", "challenge.population", "sample", "measured", "spike.in.pop", "spike.in.prop")]
colnames(coarse.gs) <- c("dataset.name", "sample.id", "sample", "vendor.sample", "measured", "spike.in.pop", "spike.in.prop")

fine.gs <- ldply(datasets, .fun = function(df) df$data)
colnames(fine.gs) <- c("dataset.name", "sample.id", "sample", "measured", "spike.in.pop", "spike.in.prop")
fine.gs <- subset(fine.gs, dataset.name %in% fine.datasets)
old.sz <- nrow(fine.gs)
fine.gs <- merge(fine.gs, fine.grained.pop.df[, c("sample", "challenge.population")], by = "sample")
new.sz <- nrow(fine.gs)
if(old.sz != new.sz) { stop(paste0("Size changed from ", old.sz, " to ", new.sz, "\n")) }
fine.gs <- fine.gs[, c("dataset.name", "sample.id", "challenge.population", "sample", "measured", "spike.in.pop", "spike.in.prop")]
colnames(fine.gs) <- c("dataset.name", "sample.id", "sample", "vendor.sample", "measured", "spike.in.pop", "spike.in.prop")

non.spike.in.cols <- c("dataset.name", "sample.id", "sample", "measured")
write.table(file = "in-silico-val-fine.csv", fine.gs[, non.spike.in.cols], row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)
write.table(file = "in-silico-val-coarse.csv", coarse.gs[, non.spike.in.cols], row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)

write.table(file = "in-silico-val-fine-spike-in-annotations.csv", fine.gs, row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)
write.table(file = "in-silico-val-coarse-spike-in-annotations.csv", coarse.gs, row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)

## this is symbol_tpm.csv
synId <- "syn21576632"
obj <- synGet(synId, downloadFile = TRUE)
cpm.expr <- read.table(obj$path, sep = ",", header = TRUE)

rownames(cpm.expr) <- as.character(cpm.expr$Gene)
cpm.expr <- cpm.expr[, !(colnames(cpm.expr) %in% c("Gene"))]

## this is ensg_tpm.csv
synId <- "syn21576631"
obj <- synGet(synId, downloadFile = TRUE)
ensg.cpm.expr <- read.table(obj$path, sep = ",", header = TRUE)

rownames(ensg.cpm.expr) <- as.character(ensg.cpm.expr$Gene)
ensg.cpm.expr <- ensg.cpm.expr[, !(colnames(ensg.cpm.expr) %in% c("Gene"))]

## Load the counts
## this is symbol_counts.csv
synId <- "syn21576630"
obj <- synGet(synId, downloadFile = TRUE)
cnts.expr <- read.table(obj$path, sep = ",", header = TRUE)

rownames(cnts.expr) <- as.character(cnts.expr$Gene)
cnts.expr <- cnts.expr[, !(colnames(cnts.expr) %in% c("Gene"))]

## this is ensg_counts.csv
synId <- "syn21576629"
obj <- synGet(synId, downloadFile = TRUE)
ensg.cnts.expr <- read.table(obj$path, sep = ",", header = TRUE)

rownames(ensg.cnts.expr) <- as.character(ensg.cnts.expr$Gene)
ensg.cnts.expr <- ensg.cnts.expr[, !(colnames(ensg.cnts.expr) %in% c("Gene"))]

coarse.gs <- ldply(datasets, .fun = function(df) df$data)
colnames(coarse.gs) <- c("dataset.name", "sample.id", "sample", "measured", "spike.in.pop", "spike.in.prop")
coarse.gs <- subset(coarse.gs, dataset.name %in% coarse.datasets)

fine.gs <- ldply(datasets, .fun = function(df) df$data)
colnames(fine.gs) <- c("dataset.name", "sample.id", "sample", "measured", "spike.in.pop", "spike.in.prop")
fine.gs <- subset(fine.gs, dataset.name %in% fine.datasets)

all.gs <- rbind(coarse.gs, fine.gs)

## Create the expression files
## First column is Gene

if(FALSE) {
admixtures <-
    dlply(all.gs,
          .parallel = TRUE,
          .variables = c("dataset.name"),
          .fun = function(df) {
              mxs <-
                  dlply(df,
                        .variables = c("sample.id"),
                        .fun = function(sdf) {
                            cols <- as.character(sdf$sample)
                            fracs <- as.numeric(sdf$measured)
                            if(!all(cols %in% colnames(cpm.expr))) {
                                print(cols)
                                stop("Missing some columns\n")
                            }
                            mat <- cpm.expr[, cols]
                            ret <- as.matrix(mat) %*% fracs
                            colnames(ret) <- sdf$sample.id[1]
                            ret
                        })
              mxs <- do.call(cbind, mxs)
              mxs <- cbind(Gene = rownames(cpm.expr), mxs)
              mxs
          })
print(warnings())
}

admixtures <-
    dlply(all.gs,
          .parallel = TRUE,
          .variables = c("dataset.name"),
          .fun = function(df) {
              mat <- as.matrix(cpm.expr)
              mxs <- create.in.silico.admixtures(mat, df,
                                                 sample.col = "sample.id", cell.type.col = "sample",
                                                 fraction.col = "measured")
              mxs <- cbind(Gene = rownames(mat), mxs)
              mxs
          })

ensg.admixtures <-
    dlply(all.gs,
          .parallel = TRUE,
          .variables = c("dataset.name"),
          .fun = function(df) {
              mat <- as.matrix(ensg.cpm.expr)
              mxs <- create.in.silico.admixtures(mat, df,
                                                 sample.col = "sample.id", cell.type.col = "sample",
                                                 fraction.col = "measured")
              mxs <- cbind(Gene = rownames(mat), mxs)
              mxs
          })

## Create the same admixtures, but for count matrices
purified.samples <- sort(as.character(unique(all.gs$sample)))

purified.cnts.expr <- cnts.expr[, purified.samples]
rownames(purified.cnts.expr) <- cnts.expr$Gene
tot.purified.cnts <- colSums(purified.cnts.expr)
median.tot.purified.cnts <- median(tot.purified.cnts)

norm.purified.cnts.expr <- sweep(purified.cnts.expr, 2, colSums(purified.cnts.expr),`/`) * median.tot.purified.cnts
rownames(norm.purified.cnts.expr) <- rownames(cnts.expr)

ensg.purified.cnts.expr <- ensg.cnts.expr[, purified.samples]
rownames(ensg.purified.cnts.expr) <- ensg.cnts.expr$Gene
ensg.tot.purified.cnts <- colSums(ensg.purified.cnts.expr)
ensg.median.tot.purified.cnts <- median(ensg.tot.purified.cnts)

ensg.norm.purified.cnts.expr <- sweep(ensg.purified.cnts.expr, 2, colSums(ensg.purified.cnts.expr),`/`) * ensg.median.tot.purified.cnts
rownames(ensg.norm.purified.cnts.expr) <- rownames(ensg.cnts.expr)

cnts.admixtures <-
    dlply(all.gs,
          .parallel = TRUE,
          .variables = c("dataset.name"),
          .fun = function(df) {
              mat <- as.matrix(norm.purified.cnts.expr)
              mxs <- create.in.silico.admixtures(mat, df,
                                                 sample.col = "sample.id", cell.type.col = "sample",
                                                 fraction.col = "measured")
              mxs <- cbind(Gene = rownames(mat), mxs)
              mxs
          })

ensg.cnts.admixtures <-
    dlply(all.gs,
          .parallel = TRUE,
          .variables = c("dataset.name"),
          .fun = function(df) {
              mat <- as.matrix(ensg.norm.purified.cnts.expr)
              mxs <- create.in.silico.admixtures(mat, df,
                                                 sample.col = "sample.id", cell.type.col = "sample",
                                                 fraction.col = "measured")
              mxs <- cbind(Gene = rownames(mat), mxs)
              mxs
          })


print(warnings())

## Create the input.csv file

cat("Writing input.csv")

input.tbl <-
    ldply(1:nrow(metadata),
          .fun = function(i) {
              ds <- as.character(metadata[i, "dataset.name"])
              ct <- as.character(metadata[i, "tumor.type"])
              params <- list(
                  "dataset.name" = ds,
                  "cancer.type" = ct,
                  "platform" = "Illumina",
                  "scale" = "Linear",
                  "normalization" = "TPM",
                  "native.probe.type" = "ENSG",
                  "symbol.compression.function" = "colMeans",
                  "ensg.compression.function" = "colMeans",
                  "symbol.to.native.mapping.file" = "native_to_hugo.tsv",
                  "ensg.to.native.mapping.file" = "native_to_ensg.tsv",
                  "hugo.expr.file" = paste0(ds, "_symbol_val_tpm.csv"),
                  "hugo.expr.est.counts.file" = paste0(ds, "_symbol_val_cnts.csv"),
                  "ensg.expr.file" = paste0(ds, "_ensg_val_tpm.csv"),
                  "ensg.expr.est.counts.file" = paste0(ds, "_ensg_val_cnts.csv"),
                  "fastq.samples" = NA,                  
                  "fastq1.files" = NA,
                  "fastq2.files" = NA,
                  "native.expr.file" = paste0(ds, "_ensg_val_tpm.csv"),
                  "native.expr.est.counts.file" = paste0(ds, "_ensg_val_cnts.csv"),
                  "symbol.compression.est.counts.function" = "colMeans",
                  "ensg.compression.est.counts.function" = "colMeans"
              )
              as.data.frame(params)
          })

parent.id <- "syn22361008"

file <- "input.csv"
write.table(file = file, input.tbl, col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
cat(paste0("Storing ", file, " to synapse\n"))
f <- File(file, parentId = parent.id, synapseStore = TRUE)
synStore(f)

file <- "fine-input.csv"
write.table(file = file, subset(input.tbl, dataset.name %in% fine.datasets), col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
cat(paste0("Storing ", file, " to synapse\n"))
f <- File(file, parentId = parent.id, synapseStore = TRUE)
synStore(f)

file <- "coarse-input.csv"
write.table(file = file, subset(input.tbl, dataset.name %in% coarse.datasets), col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
cat(paste0("Storing ", file, " to synapse\n"))
f <- File(file, parentId = parent.id, synapseStore = TRUE)
synStore(f)

nms <- names(admixtures)
names(nms) <- nms
for(nm in nms) {
    cat(paste0("Writing ", nm, " : ", nrow(admixtures[[nm]]), " x ", ncol(admixtures[[nm]]), "\n"))
    file <- paste0(nm, "_symbol_val_tpm.csv")
    write.table(file = file, admixtures[[nm]], sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
    cat(paste0("Storing ", nm, " to synapse\n"))
    f <- File(file, parentId = parent.id, synapseStore = TRUE)
    synStore(f)
}

nms <- names(ensg.admixtures)
names(nms) <- nms
for(nm in nms) {
    cat(paste0("Writing ", nm, " : ", nrow(ensg.admixtures[[nm]]), " x ", ncol(ensg.admixtures[[nm]]), "\n"))
    file <- paste0(nm, "_ensg_val_tpm.csv")
    write.table(file = file, ensg.admixtures[[nm]], sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
    cat(paste0("Storing ", nm, " to synapse\n"))
    f <- File(file, parentId = parent.id, synapseStore = TRUE)
    synStore(f)
}

nms <- names(cnts.admixtures)
names(nms) <- nms
for(nm in nms) {
    cat(paste0("Writing ", nm, " : ", nrow(cnts.admixtures[[nm]]), " x ", ncol(cnts.admixtures[[nm]]), "\n"))
    file <- paste0(nm, "_symbol_val_cnts.csv")
    write.table(file = file, cnts.admixtures[[nm]], sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
    cat(paste0("Storing ", nm, " to synapse\n"))
    f <- File(file, parentId = parent.id, synapseStore = TRUE)
    synStore(f)
}

nms <- names(ensg.cnts.admixtures)
names(nms) <- nms
for(nm in nms) {
    cat(paste0("Writing ", nm, " : ", nrow(ensg.cnts.admixtures[[nm]]), " x ", ncol(ensg.cnts.admixtures[[nm]]), "\n"))
    file <- paste0(nm, "_ensg_val_cnts.csv")
    write.table(file = file, ensg.cnts.admixtures[[nm]], sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
    cat(paste0("Storing ", nm, " to synapse\n"))
    f <- File(file, parentId = parent.id, synapseStore = TRUE)
    synStore(f)
}

cat("Done\n")
save.image(".Rdata")
cat("Done saving\n")

file <- metadata.file
f <- File(file, parentId = parent.id, synpaseStore = TRUE)
synStore(f)

## Write out the admixtures and store to synapse
l <- list("in-silico-val-fine.csv" = fine.gs, "in-silico-val-coarse.csv" = coarse.gs,
          "in-silico-val-fine-spike-in-annotations.csv" = fine.gs, "in-silico-val-coarse-spike-in-annotations.csv" = coarse.gs)
for(nm in names(l)) {
    file <- nm
    f <- File(file, parentId = parent.id, synapseStore = TRUE)
    synStore(f)
}

cat("Exiting\n")
q(status = 0)
