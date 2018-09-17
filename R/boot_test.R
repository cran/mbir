#' Bootstrap Confidence Intervals via Resampling
#'
#'Provides nonparametric confidence intervals via percentile-based resampling.
#'
#'@param x,y numeric vectors of data values
#'@param conf.int (optional) confidence level of the interval. Defaults to \code{0.90}
#'@param resample (optional) number of resamples. Defaults to 10,000
#'@param med (optional) number indicating true difference in medians to test against. Defaults to zero.
#'@details Refer to vignette for further information.
#'@examples require(graphics)
#'
#'@examples a <- rnorm(25, 80, 35)
#'@examples b <- rnorm(25, 100, 50)
#'
#'@examples boot_test(a, b, 0.95, 10000)
#'@export

boot_test<-function(x,y,conf.int,resample, med){

  if(is.character(x) == TRUE || is.factor(x) == TRUE || is.character(y) == TRUE || is.factor(y) == TRUE){
    error<-"Sorry, data must be numeric or integer values."
    stop(error)
  }

  if(missing(conf.int)){
    conf.int<-.9
  }

  if(missing(resample)){
    resample<-10000
  }

  if(missing(med)){
    med<-0
  }
  medx <- medy <- numeric(resample)

  for (i in 1:resample) {
    medx[i] <- stats::median(sample(x, replace=TRUE),na.rm = T)
    medy[i] <- stats::median(sample(y, replace=TRUE),na.rm = T)
  }

  Diff <- medx - medy

  interval<-(1-conf.int)/2

  level<-paste(as.character(100*conf.int),"%",sep = "")

  cat("Median of ",deparse(substitute(x))," = ",round(stats::median(x,na.rm = T),digits = 2),"; ","Median of ",deparse(substitute(y))," = ",round(stats::median(y,na.rm = T),digits = 2),"\n\n",sep = "")

  cat(level," Bootstrap Confidence Interval\n",format(resample, big.mark=",", scientific=FALSE)," Resamples\n\n",sep = "")

  ci<-round(stats::quantile(Diff, c(interval, 1-interval),na.rm=T),digits = 3)
  print(ci)
  if(unname(ci[1]) < med && unname(ci[2]) > med){
    cat("\nDifference in Medians: Lacking Evidence (CI contains ",med,").",sep = "")
    Inference <- "Lacking Evidence"}
  else {
    cat("\nDifference in Medians: Evidence Present (CI does not contain ",med,").",sep = "")
    Inference <- "Evidence Present"}

  graphics::hist(Diff, col='gray', border='white', las=1, xlab = "Difference in Medians", main = " ")
  graphics::abline(v=med, lty=2)
  graphics::abline(v=c(unname(ci[1])))
  graphics::abline(v=unname(ci[2]))

  invisible(list(
    med = med,
    med.diff = stats::median(Diff),
    b.LL = unname(ci[1]),
    b.UL = unname(ci[2]),
    Inference = Inference)
  )
}
