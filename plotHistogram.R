## plotHistogram.R is a script for Measurement Histogram
source("bmarray.R")
library(ggplot2)
## Function for arranging ggplots. use png(); arrange(p1, p2, ncol=1); dev.off() to save.
require(grid)

vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange_ggplot2 <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
  dots <- list(...)
  n <- length(dots)
  if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
  if(is.null(nrow)) { nrow = ceiling(n/ncol)}
  if(is.null(ncol)) { ncol = ceiling(n/nrow)}
  ## NOTE see n2mfrow in grDevices for possible alternative
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
  ii.p <- 1
  for(ii.row in seq(1, nrow)){
    ii.table.row <- ii.row	
    if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
    for(ii.col in seq(1, ncol)){
      ii.table <- ii.p
      if(ii.p > n) break
      print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
      ii.p <- ii.p + 1
    }
  }
}

plotHistogram <- function(object, type = c("M", "beta"), width = 0.01, cond = FALSE, sampleRatio = 3000) {
  type = match.arg(type)
  m <- getBeta(object, type = type)
  sample <- dim(m)[2]
  len    <- sampleRatio
  m <- m[sample(dim(m)[1])[1:len], , drop = FALSE]
  if(cond){
    df <- data.frame(M = as.vector(m), cond = factor(c(rep("tumor", len * sample / 2), rep("normal", len * sample / 2))))
    p  <- ggplot(df, aes(x=M, fill=cond)) + geom_histogram(binwidth=width, stat="density", alpha=.5) + geom_density(alpha=.6)
    return(p)
  }
  else{
    df <- data.frame(M = as.vector(m))
    p <- ggplot(df, aes(x = M)) + geom_histogram(binwidth=width, stat="density", aes(y = ..density..), colour="grey") + geom_density(alpha=.5, color = "blue")
    return(p)
  }
}

plotDiff <- function(object, type = c("M", "beta"), width = 0.01, sampleRatio = 3000) {
  type = match.arg(type)
  m <- getBeta(object, type = type)
  sample <- dim(m)[2]
  len    <- sampleRatio
  m <- m[sample(dim(m)[1])[1:len], , drop = FALSE]
  ssize <- sample / 2
  diff <- abs(m[, 1: ssize, drop = FALSE] - m[, (ssize+1):sample, drop = FALSE])
  df <- data.frame(D = as.vector(diff))
  p1 <- ggplot(df, aes(x = D)) + geom_histogram(binwidth=width, stat="density", aes(y = ..density..), colour="grey") + geom_density(alpha=.5, color = "blue")
  ### log lm
  hdat <- hist(diff, breaks=ceiling(1 / width), plot=FALSE)
  prob <- log(hdat$counts / sum(hdat$counts))  # probability of each region
  logbeta <- hdat$mids                         # mids of each region stands for its value
  idx = (prob != -Inf)
  df2 <- data.frame(logB = logbeta[idx], prob = prob[idx])
  p2  <- ggplot(df2, aes(x=logB, y=prob)) + geom_point(shape=1) + geom_smooth(method=lm) + ggtitle("measure/Log ~ probability")
  p   <- list(p1 = p1,p2 = p2, corr = lm(prob ~ logB, data = df2))
  return(p)
}

plotNull <- function(object, width = 0.01){
  nulldat <- abs(unlist(object$null$value)) * unlist(object$null$length)
  df <- data.frame(Null = nulldat)
  p <- ggplot(df, aes(x = Null)) + geom_histogram(binwidth=width, stat="density", aes(y = ..density..), colour="grey") + geom_density(alpha=.5, color = "red")
  return(p)
}
