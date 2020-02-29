library(pacman)

p_load(plyr)
p_load(dplyr)
p_load(ggplot2)

files <- list("coarse-round1" = "round1-coarse-results.tsv",
              "fine-round1" = "round1-fine-results.tsv")


files <- list("coarse-round2" = "round2-coarse-results.tsv",
              "fine-round2" = "round2-fine-results.tsv")

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

## I confirmed the moment-matching equations here:
## https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
## where a = alpha = shape1
## and   b = beta = shape2
## a = mu^2 * [ ( (1 - mu)/(sigma^2) ) - (1 / mu) ]
## b = a* (1 - mu) / mu

fit.beta <- function(x) {
    mu <- mean(x)
    var <- var(x)
    a = (1 - mu) / (var^2)
    a = a - (1 / mu)
    a = a * mu^2
    b = a * (1 - mu) / mu
    
    ##    list(shape1 = a, shape2 = b)
    method <- "mle"
    fw <- fitdist(x, "beta", start=list(shape1 = a, shape2 = b), method = method)
    list(shape1 = as.numeric(fw[[1]]["shape1"]), shape2 = as.numeric(fw[[1]]["shape2"]))

}

fit.norm <- function(x) {
    mu <- mean(x)
    std.dev <- sd(x)
    list(mean = mu, sd = std.dev)
}

dbeta.wrapper <- function(xs, params) {
    ys <- dbeta(xs, params$shape1, params$shape2)
    ys
}

dnorm.wrapper <- function(xs, params) {
    ys <- dnorm(xs, params$mean, params$sd)
    ys
}


library(fitdistrplus) ## for fitdist

fit.weibull <- function(x) {
    fw <- fitdist(x, "weibull")
    list(shape = as.numeric(fw[[1]]["shape"]), scale = as.numeric(fw[[1]]["scale"]))
}

dweibull.wrapper <- function(x, params) {
    ys <- dweibull(x, shape = params$shape, scale = params$shape, log = FALSE)
    ys
}

library(extraDistr) # for dkumar

library(VGAM) # for kumar

## kumaraswamy
fit.kumaraswamy <- function(x) {
    kdata <- data.frame(y = x)
    fit <- vglm(y ~ 1, kumar(lshape1 = "identitylink", lshape2 = "identitylink"), data = kdata, trace = FALSE)
    ## shape1 is in the column 1; shape2 is in column 2
    shape1 <- coef(fit, matrix=TRUE)[1,1]
    shape2 <- coef(fit, matrix=TRUE)[1,2]
    ##    list(a = shape1, b = shape2)
    list(shape1 = shape1, shape2 = shape2)
}

dkumaraswamy.wrapper <- function(x, params) {
    ## ys <- dkumar(x, a = params$a, b = params$b, log = FALSE)
    ys <- VGAM::dkumar(x, shape1 = params$shape1, shape2 = params$shape2, log = FALSE)
    ys
}



nms <- names(files)
for(nm in nms) {
    orig.tsv <- read.table(files[[nm]], sep="\t", header=TRUE, as.is=TRUE)
    for(scoring.metric in c("pearson", "spearman")) {
        tsv <- subset(orig.tsv, is_latest == "true" | grepl(repo_name, pattern = "baseline"))
        tsv <- subset(tsv, metric == scoring.metric)

        fit.methods <- list("norm" = "fit.norm")
        prob.densities <- list("norm" = "dnorm.wrapper")

        fit.methods <- list("weibull" = "fit.weibull")
        prob.densities <- list("weibull" = "dweibull.wrapper")
        
        fit.methods <- list("kumaraswamy" = "fit.kumaraswamy")
        prob.densities <- list("kumaraswamy" = "dkumaraswamy.wrapper")

        fit.methods <- list("none" = NULL)
        prob.densities <- list("none" = NULL) 


        fit.methods <- list("kumaraswamy" = "fit.kumaraswamy",
                            "norm" = "fit.norm",
                            "weibull" = "fit.weibull",
                            "beta" = "fit.beta",
                            "none" = NULL)

        prob.densities <- list("kumaraswamy" = "dkumaraswamy.wrapper",
                               "norm" = "dnorm.wrapper",
                               "weibull" = "dweibull.wrapper",
                               "beta" = "dbeta.wrapper",
                               "none" = NULL)

        ## Stay away from the boundary by add eps
        eps <- 10^-5
        shifts <- list("beta" = 1+eps, "norm" = 0, "weibull" = 1+eps, "kumaraswamy" = 1+eps, "none" = 0)
        scales <- list("beta" = 2+2*eps, "norm" = 1, "weibull" = 1, "kumaraswamy" = 2+2*eps, "none" = 1)
        
        for(dist in names(fit.methods)) {
            fits <- NULL
            if(dist != "none") {
                fits <-
                    dlply(tsv,
                          .variables = c("dataset"),
                          .fun = function(df.ds) {
                              ret.ds <-
                                  dlply(df.ds,
                                        .variables = c("celltype"),
                                        .fun = function(df) {
                                            df
                                            ## params <- fit.beta(df$metric_value)
                                            vals.shifted <- (df$metric_value + shifts[[dist]]) / scales[[dist]]
                                            params <- do.call(fit.methods[[dist]], list(vals.shifted))
                                            params
                                            
                                        })
                          })
            }
            plts <-
                dlply(tsv,
                      .variables = c("dataset"),
                      .fun = function(df) {
                          ## df$metric_value <- ( df$metric_value + shifts[[dist]] ) / scales[[dist]]
                          ## fits.dist <- fits[[df$dataset[1]]]
                          fits.dist <- NULL
                          if(dist != "none") {
                              fits.dist <- fits[[df$dataset[1]]]
                          }
                          ## nms <- names(fits.dist)
                          nms <- as.character(unique(df$celltype))
                          names(nms) <- nms
                          hist.maxes <-
                              llply(nms,
                                    .fun = function(ct) {
                                        g <- ggplot(data = subset(df, celltype == ct))
                                        g <- g + xlim(c(-1,1))
                                        ## g <- g + xlim((c(-1,1)+shifts[[dist]]) / scales[[dist]])
                                        g <- g + geom_histogram(aes(metric_value))
                                        tbl <- ggplot_build(g)$data[[1]]
                                        max(tbl$y)
                                    })
                          all.fits <- NULL
                          if(!(dist == "none")) {
                              all.fits <-
                                  ldply(nms,
                                        .fun = function(ct) {
                                            ## xs <- ( seq(from = -1, to = 1, by = 0.02) + shifts[[dist]] ) / scales[[dist]]
                                            xs <- seq(from = -1, to = 1, by = 0.02)
                                            ## ys <- dbeta(xs, fits.dist[[ct]]$shape1, fits.dist[[ct]]$shape2)
                                            xs.shifted <- (xs + shifts[[dist]])/scales[[dist]]
                                            ys <- do.call(prob.densities[[dist]], list(xs.shifted, fits.dist[[ct]]))
                                            ## rescale
                                            ys <- ys / max(ys)
                                            ys <- ys * hist.maxes[[ct]]
                                            data.frame(x = xs, y = ys)
                                        })
                              colnames(all.fits)[1] <- "celltype"
                          }
                          g <- ggplot(df)
                          g <- g + geom_histogram(aes(metric_value))
                          ## g <- g + xlim(( c(-1,1)+shifts[[dist]] ) / scales[[dist]])
                          g <- g + xlim(c(-1,1))
                          g <- g + facet_wrap(~ celltype, nrow = 2, scale = "free_x")
                          g <- g + xlab(paste0(firstup(scoring.metric), " Correlation (Predicted vs Ground Truth)"))
                          title <- paste0(nm, " ", dist, " ", scoring.metric, " ", df$dataset[1], " (n = ", length(unique(df$repo_name)), ")") 
                          g <- g + ggtitle(title)
                          if(dist != "none") {
                              g <- g + geom_line(data = all.fits, aes(x = x, y = y))
                          }
                          file <- paste0(nm, "-", scoring.metric, "-", df$dataset[1], "-", dist, "-fits.pdf")
                          pdf(file)
                          print(g)
                          d <- dev.off()

                          file <- paste0(nm, "-", scoring.metric, "-", df$dataset[1], "-", dist, "-fits.png")
                          png(file)
                          print(g)
                          d <- dev.off()
                          
                      })
            }
    }
}

if(FALSE) {
r <- as.numeric(na.omit(pos.curves$l1.norm.norm))
r <- r / 100
d <- density(r)
fit <- fitdistr(r,"gamma", list(shape=1.5, rate=2))
dy <- dgamma(d$x, shape=1.5, rate=2)
dy <- dgamma(d$x, shape=fit$estimate[1], rate=fit$estimate[2]) 
fit <- fitdist(r, "gamma", method="qme", probs = c(0.25,0.75))
dy <- dgamma(d$x, shape=as.numeric(fit$estimate["shape"]), rate=as.numeric(fit$estimate["rate"]))
plot(d)
lines(d$x, dy * max(d$y) / max(dy))

fit <- fitdistr(r,"beta",list(shape1=0.5,shape2=0.5))
dy <- dbeta(d$x, shape1=fit$estimate[1], shape2=fit$estimate[2]) 
plot(d)
lines(d$x, dy * max(d$y) / max(dy))
}
