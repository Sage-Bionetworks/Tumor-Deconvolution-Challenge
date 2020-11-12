
if(FALSE) {
    both <- subset(res.all, ( method.name=="IZI" | method.name == "CIBERSORTx") & submission == "1" & subchallenge == "coarse")
    ## csx <- subset(res.all, method.name=="CIBERSORTx" & submission == "1" & subchallenge == "coarse")

    calc.noise <- function(tbl) {
        ddply(tbl,
              .variables = c(method.name.col, round.col, subchallenge.col, dataset.name.col, sample.id.col),
              .fun = function(df) {
                  data.frame(other = 1 - sum(df[, measured.col], na.rm=TRUE))
                             
              })
    }
    both.noise <- calc.noise(both)
    both.with.noise <- merge(both, both.noise)

    calc.cors <- function(tbl) {
        cors <-
            ddply(tbl,
                  .variables = c(method.name.col, round.col, subchallenge.col, dataset.name.col, cell.type.col),
                  .fun = function(df) {
                      if(nrow(df) <= 2) { return(NULL) }
                      data.frame(cor = cor(df[, prediction.col], df[, measured.col]))
                  })
    }

    cutoffs <- list("< 0.2" = 0.2, "< 0.4" = 0.4, "< 0.6" = 0.6, "< 0.8" = 0.8)

    noise.vs.cors <-
        ldply(cutoffs,
              .fun = function(cutoff) {
                  sb <- subset(both.with.noise, other < cutoff)
                  cors <- calc.cors(sb)
                  cors[!is.na(cors$cor),]
              })
    colnames(noise.vs.cors)[1] <- "noise"
    mean.cors <-
        ldply(cutoffs,
              .fun = function(cutoff) {
                  sb <- subset(both.with.noise, other < cutoff)
                  cors <- calc.cors(sb)                  
                  ret <- ddply(cors, .variables = c(method.name.col),
                               .fun = function(df) data.frame(mean.cor = mean(df$cor, na.rm=TRUE)))
                  colnames(ret)[1] <- "method.name"
                  ret
              })
    colnames(mean.cors)[1] <- "cutoff"
    g <- ggplot(data = noise.vs.cors, aes(x = cell.type, y = cor, colour = method.name))
    g <- g + geom_boxplot()
    ## g <- g + geom_beeswarm()
    g <- g + facet_wrap(~ noise)
    g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    g <- g + geom_abline(slope = 0, intercept = 0.9, linetype = "dashed")

    g2 <- ggplot(data = mean.cors, aes(y = mean.cor, x = cutoff))
    g2 <- g2 + geom_col()
    g2 <- g2 + facet_wrap(~ method.name)
    m <- merge(subset(noise.vs.cors, noise==0.8), subset(noise.vs.cors, noise==0.6), by = c("method.name", "submission", "subchallenge", "dataset.name", "cell.type"))
    plot(m$cor.x, m$cor.y)
    m[m$cor.x != m$cor.y,]
flag <- izi.cors$submission == 1
## plot(izi[flag, measured.col], izi[flag, prediction.col])
plot(izi.cors[flag, "other"], izi.cors[flag, "cor"])
}

inferred.output.type <-
    ddply(res.all,
          .variables = c(method.name.col, round.col, subchallenge.col),
          .fun = function(res.meth) {
              sums <-
                  ddply(res.meth,
                        .variables = c(round.col, dataset.name.col, sample.id.col),
                        .fun = function(res.meth.sample) {
                            sm <- sum(res.meth.sample[, prediction.col])
                            if( (sm > 1.2) & (res.meth[1,method.name.col] == target.meth)) {
                                print(res.meth.sample)
                                stop("stop")
                            }
                            data.frame(sums = sm, num = length(!is.na(res.meth.sample[,prediction.col])))
                        })
              na.rm <- FALSE
              min.sum <- min(sums$sum, na.rm=na.rm)
              max.sum <- max(sums$sum, na.rm=na.rm)
              min.num <- min(sums$num)
              eps <- 0.01
              min.pred <- min(res.meth[, prediction.col], na.rm=na.rm)
              max.pred <- max(res.meth[, prediction.col], na.rm=na.rm)
              output.type <- NA
              flag <- min.pred < -eps
              flag <- flag || (max.pred > (1 + eps))
              if(flag) {
                  output.type <- "score"
              } else if(abs(max.sum - 1) < eps) {
                  output.type <- "constrained.fraction"
              } else {
                  output.type <- "fraction"
              }
              data.frame(min.pred = min.pred,
                         max.pred = max.pred,
                         max.sum = max.sum,
                         min.sum = min.sum,
                         output.type = output.type,
                         non.negative = min.pred >= 0,
                         sums.to.one = ((min.sum >= (1-eps)) && (max.sum <= (1+eps))),
                         sum.lt.one = min.sum <= (1 - eps),
                         any.cell.gt.one = max.pred >= 1 + eps)
          })
inferred.output.type <- subset(inferred.output.type, !(submission == "latest"))
