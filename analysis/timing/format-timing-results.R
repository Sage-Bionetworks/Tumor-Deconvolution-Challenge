## The predictions are those from the 'rerun-validation' / 'post-competitive' phase, i.e.,
## where the coarse- and fine-grained challenge use the same data.
## This is as opposed to the original competitive phase (against which the
## methods were ranked) where the coarse- and fine-grained challenges
## differed.

predictions.synId <- "syn22314641"

## Read in the "post-competitive" predictions 
obj <- synGet(predictions.synId, downloadFile=TRUE)
res.all <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)
