#' Smallest Worthwhile Change: Individual
#'
#'Provides longitudinal magnitude-based inferences for an individual's change from previous time point and magnitude of deviation from trend line.
#'
#'@param x numeric vectors of data values
#'@param swc smallest worthwhile change
#'@param type (character) indicator specifying which type of analysis: "previous" or "trend"
#'@param ts (required if \code{type = "trend"}) target slope
#'@param te (optional) typical error. Defaults to typical error of the estimate
#'@param main (optional) plot title. Defaults to blank
#'@param xlab (optional) x-axis label. Defaults to "Measurement"
#'@param ylab (optional) y-axis label. Defaults to name of \code{x}
#'@details Refer to vignette for further information.
#'@references Hopkins WG. (2017). A spreadsheet for monitoring an individual's changes and trend. \emph{Sportscience} 21, 5-9. sportsci.org/2017/wghtrend.htm
#'@examples df<-c(97.5,99.9,100.2,101,101.2,99.8)
#'
#'@examples swc_ind(x = df, swc = 0.5, te = 1, ts = 0.25, type = "trend")
#'@export

swc_ind<-function(x, swc, type=c("previous","trend"), ts, te, main, xlab, ylab){


  if(sum(x < 0) > 0){
    error<-"Sorry, positive values only."
    stop(error)
  }

  if(is.character(x) == TRUE || is.factor(x) == TRUE || is.character(swc) == TRUE || is.factor(swc) == TRUE){
    error<-"Sorry, data must be numeric or integer values."
    stop(error)
  }

  if(length(x) < 4){
    error<-"Sorry, not enough data."
    stop(error)
  }

  if(missing(type)){
    error<-"Missing Argument: Please specify type of analysis: previous or trend."
    stop(error)
  }

  if(type == "previous" || type == "trend"){
    type<-type
  } else {
    error<-"Incorrect Argument. Types: previous or trend"
    stop(error)
  }

  if(missing(swc)){
    error<-"Missing Argument: Please specify Smallest Worthwhile Change."
    stop(error)
  }

  if(type=="trend" && missing(ts)){
    error<-"Missing Argument: Please specify Target Slope."
    stop(error)
  }

  if(missing(ts)){
    ts<-0
  } else { ts<-ts}

  if(ts < 0){
    error<-"Sorry, target slope must be positive value."
    stop(error)
  }

  if(missing(main)){
    main<-" "
  }

  if(missing(xlab)){
    xlab<-"Measurement"
  }

  if(missing(ylab)){
    ylab<-deparse(substitute(x))
  }

  times<-c(1:length(x))
  mylen<-times[-utils::tail(times,n=1)]
  times2<-c(2:length(x))

  df<-length(x)-1
  offset<-data.frame(x,blank=1)
  offset$diff <- stats::ave(offset$x, offset$blank, FUN=function(x) c(0, diff(x)))
  offset$diffn<-offset$diff*-1
  point<-offset$diff[-1]
  point2<-offset$diffn[-1]
  prev<-offset$x[1:length(offset$x)-1]

  X <- cbind(1, times)
  colnames(X)[1] <- "Intercept"
  y <- x
  A <- base::solve(t(X) %*% X)
  beta <- (A %*% t(X) %*% y)
  slope <- beta[[2]]
  resid <- (y - X %*% beta)
  fit <- (y - resid)
  residual_var <- base::sum(resid^2) / (length(x) - 1 - 1)
  residual.scale <- sqrt(residual_var)
  rss<-base::sum(resid^2)
  mss <- base::sum((fit - mean(fit))^2)
  r.squared <- mss/(mss + rss)
  adj.r.squared <- 1 - (1 - r.squared) * ((length(x) - 1)/(length(x) - 1 - 1))
  df.res <- (length(x) - 1 - 1)
  beta_covar <- (residual_var * A)
  beta_SE <- base::sqrt(base::diag(beta_covar))
  T.stat <- slope/beta_SE
  p.value <- 2 * stats::pt(-base::abs(T.stat), length(x)-2)
  tha.p.value <- p.value[[2]]
  f.val <- (mss/(1))/residual_var

  trivband.plus<-fit+swc*fit/100+residual.scale*sqrt((1/length(times))+(1/(length(times)-1))*(times-mean(times))^2/stats::sd(times)^2)
  trivband.minus<-fit-swc*fit/100-residual.scale*sqrt((1/length(times))+(1/(length(times)-1))*(times-mean(times))^2/stats::sd(times)^2)


  if(missing(te)){
    te<-residual.scale
  } else {
    te<-te
    }

  yrange<-ifelse(swc >= te,swc,te)

  if(type == "previous"){

    pos<-round(100*(ifelse(point > 0,stats::pt(abs((point-swc)/sqrt((te^2)+(te^2))),df),stats::pt(-abs((point-swc)/sqrt((te^2)+(te^2))),df))))
    neg<-round(100*(ifelse(point2 > 0,stats::pt(abs((point2-swc)/sqrt((te^2)+(te^2))),df),stats::pt(-abs((point2-swc)/sqrt((te^2)+(te^2))),df))))

    pos<-ifelse(point > 0 & swc > point,100-pos,pos)
    neg<-ifelse(point < 0 & swc > abs(point),100-neg,neg)

    triv<-100-pos-neg

  mbi<-ifelse(pos > 10 & neg > 10,"---",
              ifelse(pos <= 10 & neg <= 10,"Trivial",
              ifelse(neg > 10 & triv >= 10,"Possible Decrease",
                     ifelse(pos > 10 & triv >= 10,"Possibe Increase",
                            ifelse(pos < 90 & triv < 10 & neg < 10,"Likely Increase",
                                   ifelse(pos < 10 & triv < 10 & neg < 90,"Likely Decrease",
                                          ifelse(pos <= 10 & neg >= 90," Very Likely Decrease",
                                                 ifelse(pos >= 90 & neg <= 10," Very Likely Increase"," "))))))))

  noway<-data.frame(mylen,"-",times2," ",point," ",neg,triv,pos," ",mbi)
  names(noway)<-c(" "," "," "," ","Diff"," ","N","T","P"," ","MBI")
  cat("   MBI From Previous Time Point:\n\n")
  print(noway,row.names=F)
  graphics::plot(times,x,ylim=range(c(x-yrange, x+yrange)),ylab = ylab,xlab=xlab,main=main,pch=19,type="b")
  graphics::lines(times,trivband.plus,lty=2,col="#666666")
  graphics::lines(times,trivband.minus,lty=2,col="#666666")
  graphics::arrows(times, x-te, times, x+te, length=0.05, angle=90, code=3)
  }


  if(type == "trend"){

  pred.err<-sqrt(te^2+residual.scale^2*((1/length(times)+1/(length(times)-1))*(times-mean(times))^2/sd(times)^2))

  degfree<-(residual.scale^2*(((1/length(times)+1/(length(times)-1))*(times-mean(times))^2/stats::sd(times)^2)+(te*fit/100)^2)^2)/((residual.scale^2*((1/length(times)+1/(length(times)-1))*(times-mean(times))^2/stats::sd(times)^2))^2/(mean(times)-2)+(te*fit/100)^4/(df.res-1))
  diff<-(x-fit)
  diffo<-(fit-x)

  if(length(x) <= 8){
  pos<-round(100*(ifelse(diff < 0,stats::pt(-abs((diff-(swc*fit/100))/sqrt(te/2^2+te/2^2)),degfree),1-stats::pt(-abs((diff-(swc*fit/100))/sqrt(te^2+te/2^2)),degfree))),0)
  neg<-round(100*(ifelse(diff < 0,stats::pt(abs((diffo-(swc*fit/100))/sqrt(te/2^2+te/2^2)),degfree),1-stats::pt(abs((diffo-(swc*fit/100))/sqrt(te^2+te/2^2)),degfree))),0)
  } else {
  pos<-round(100*(ifelse(diff < 0,stats::pt(-abs((diff-(swc*fit/100))/sqrt(residual.scale^2+pred.err^2)),degfree),1-stats::pt(-abs((diff-(swc*fit/100))/sqrt(residual.scale^2+pred.err^2)),degfree))),0)
  neg<-round(100*(ifelse(diff < 0,stats::pt(abs((diffo-(swc*fit/100))/sqrt(residual.scale^2+pred.err^2)),degfree),1-stats::pt(abs((diffo-(swc*fit/100))/sqrt(residual.scale^2+pred.err^2)),degfree))),0)
  }


  pos<-ifelse(diff > 0 & swc > diff,100-pos,pos)

  neg<-ifelse(diff < 0 & swc > abs(diff),100-neg,neg)


  triv<-round(100-pos-neg,0)
  neg<-ifelse(triv < 0,abs(triv),neg)
  triv<-round(100-pos-neg,0)

  diff<-round(x-fit,1)
  mbi<-ifelse(pos > 10 & neg > 10,"---",
              ifelse(neg <= 10 & pos <= 10,"Trivial",
              ifelse(neg > 10 & triv >= 10,"Possibly Lower",
                     ifelse(pos > 10 & triv >= 10,"Possibly Higher",
                            ifelse(pos < 90 & triv < 10 & neg < 10,"Likely Higher",
                                   ifelse(pos < 10 & triv < 10 & neg < 90,"Likely Lower",
                                          ifelse(pos <= 10 & neg >= 90," Very Likely Lower",
                                                 ifelse(pos >= 90 & neg <= 10," Very Likely Higher"," "))))))))

  cat("   Trend Parameters:")
  table1<-matrix(c(round(unname(slope),digits = 3),
                   round(unname(adj.r.squared),digits=3),
                   round(unname(f.val),digits=3),
                   round(unname(tha.p.value),digits = 3)),ncol =  1,byrow = T)
  rownames(table1)<-c("Slope","R-squared","F stat","P value")
  colnames(table1)<-c(" ")
  print(table1)
  cat("\n")
  diff2<-(unname(slope)-ts)
  diff3<-(unname(-slope)-ts)

  if(ts < abs(slope)){
  pos2<-round(100*ifelse(slope > 0,stats::pt(abs((diff2)/tha.p.value),df.res),1-stats::pt(abs((diff2)/tha.p.value),df.res)),0)
  neg2<-round(100*ifelse(slope > 0,stats::pt(-abs((diff3)/tha.p.value),df.res),1-stats::pt(-abs((diff3)/tha.p.value),df.res)),0)
  } else {
    pos2<-round(100*ifelse(slope > 0,1-stats::pt(abs((diff2)/tha.p.value),df.res),stats::pt(abs((diff2)/tha.p.value),df.res)),0)
    neg2<-round(100*ifelse(slope > 0,stats::pt(-abs((diff3)/tha.p.value),df.res),1-stats::pt(abs((diff3)/tha.p.value),df.res)),0)
  }



  triv2<-round(100-pos2-neg2,digits = 1)
  neg2<-ifelse(triv2 < 0,abs(triv2),neg2)
  triv2<-round(100-pos2-neg2,0)


  decrease<-ifelse(neg2 < 10,"Very Unlikely",
               ifelse(neg2 < 25,"Unlikely",
                      ifelse(neg2 < 75,"Possibly",
                             ifelse(neg2 < 90,"Likely",
                                    ifelse(neg2 >= 90,"Very Likely")))))
  none<-ifelse(triv2 < 10,"Very Unlikely",
                   ifelse(triv2 < 25,"Unlikely",
                          ifelse(triv2 < 75,"Possibly",
                                 ifelse(triv2 < 90,"Likely",
                                        ifelse(triv2 >= 90,"Very Likely")))))

  increase<-ifelse(pos2 < 10,"Very Unlikely",
                   ifelse(pos2 < 25,"Unlikely",
                          ifelse(pos2 < 75,"Possibly",
                                 ifelse(pos2 < 90,"Likely",
                                        ifelse(pos2 >= 90,"Very Likely")))))
  table2<-matrix(c("Decrease","Trivial","Increase",neg2,triv2,pos2),nrow = 2,byrow = T)
  rownames(table2)<-c(" ","MBI (%)")
  colnames(table2)<-c(decrease,none,increase)
  print(table2)
  graphics::plot(times,x,ylim=range(c(x-yrange, x+yrange)),ylab = ylab,xlab=xlab,main=main,pch=19,type="b")
  graphics::abline(beta[1],slope,col="red")
  graphics::lines(times,trivband.plus,lty=2)
  graphics::lines(times,trivband.minus,lty=2)
  graphics::arrows(times, x-te, times, x+te, length=0.05, angle=90, code=3)
  noway<-data.frame(times," ",diff," ",neg,triv,pos," ",mbi)
  names(noway)<-c("Point"," ","Diff"," ","N","T","P"," ","MBI")
  cat("\n   MBI From Trend Line:\n\n")
  print(noway,row.names=F)
  rval <- list(slope=table1[1], rsq=table1[2],fstat=table1[3],p.value=table1[4])
  }
}




