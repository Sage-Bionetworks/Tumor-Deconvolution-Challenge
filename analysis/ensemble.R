library(ConsRank)
crs <- subset(results[["latest"]][["res.round"]], subchallenge == "coarse")
fine <- subset(results[["latest"]][["res.round"]], subchallenge == "fine")

res <- results[["1"]][["res.round"]]

vars <- c("subchallenge", "method.name", "submission", "dataset.name", 
          "sample.id", "cell.type")
vars <- colnames(res)
vars <- vars[!(vars %in% c("measured", "prediction", "sample.id"))]
tres <- 
  ddply(res, .variables = vars,
        .fun = function(df) {
                 data.frame(measured = df$measured, prediction = df$prediction,
                            pred.rank = rank(df$prediction, ties.method="first"),
                            sample.id = df$sample.id)
               })

exclude.vars <- c("measured", "prediction", "pred.rank", "method.name",
                  "objectId", "submitterId", "repo_name", "comparator")
vars2 <- colnames(res)[!(colnames(res) %in% exclude.vars)]

mean.res <-
  ddply(tres, .variables = vars2,
        .fun = function(df) {
                 m <- unique(df$measured)
                 if(length(m) != 1) { print(m); stop("") }
                 data.frame(measured = m[1], prediction = mean(df$pred.rank))
               })
vars <- colnames(mean.res)
vars <- vars[!(vars %in% c("measured", "prediction", "sample.id"))]
tmp <-
  ddply(mean.res, .variables = vars,
                  .fun = function(df) {
                           data.frame(spearman = cor(df$prediction, df$measured, method="spearman"))
                         })

# Average over cell type first, then over data set
vars <- colnames(tmp)
vars <- vars[!(vars %in% c("cell.type", "spearman"))]
scores <-
  ddply(tmp, .variables = vars,
             .fun = function(df) {
                      data.frame(spearman = mean(df$spearman))
                    })

vars <- colnames(scores)
vars <- vars[!(vars %in% c("cell.type", "spearman", "dataset.name", "tumor.type", "mixture.type", "distribution.type"))]
scores <-
  ddply(scores, .variables = vars,
             .fun = function(df) {
                      data.frame(spearman = mean(df$spearman))
                    })

y <- subset(tres, cell.type=="B" & dataset.name == "DS4")
x <- acast(y[, c("pred.rank", "method.name", "sample.id")], formula = "method.name ~ sample.id", value.var = "pred.rank")

res <- consrank(x, proc=TRUE)
res <- consrank(x, proc=TRUE, algorithm="quick")
