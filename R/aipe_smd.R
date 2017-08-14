#' Accuracy in Parameter Estimation: Standardized Mean Difference
#'
#'Estimates sample size for paired or independent, two-sample study desings via Accuracy in Parameter Estimation. Calculates \emph{n} so a given study is likely to obtain margin of error no larger than chosen target margin of error.
#'
#'@param moe target margin of error in standard deviation units
#'@param paired (character) logical indicator specifying if \code{x} and \code{y} are paired \code{(TRUE)} or independent \code{(FALSE)}
#'@param conf.int (optional) confidence level of the interval. Defaults to \code{0.90}
#'@param assur.lvl (optional) desired level of \emph{assurance} (percent experiments whose MOE is less than target MOE). Defaults to \code{0.99}
#'@param r (required if \code{paired = TRUE}) population correlation between the two measures
#'@details Refer to vignette for further information.
#'@references Maxwell SE, Kelley K & Rausch JR. (2008). Sample size planning for statistical power and accuracy in parameter estimation. \emph{Annual Review of Psychology}, 59, 537-563.
#'@references Kelley K & Rausch JR. (2006). Sample size planning for the standardized mean difference: Accuracy in parameter estimation via narrow confidence intervals. \emph{Psychological Methods}, 11, 363–385.
#'@examples aipe_smd(moe = 0.55, paired = TRUE, conf.int = .9, assur.lvl = .99, r = 0.75)
#'@export

aipe_smd<-function(moe, paired=c(TRUE,FALSE), conf.int, assur.lvl, r){

  if(is.character(moe) == TRUE || is.factor(moe) == TRUE){
    error<-"Sorry, data must be numeric or integer values."
    stop(error)
  }

  if(missing(paired)){
    error<-"Missing Argument: Please specify if data is paired = TRUE/FALSE."
    stop(error)
  }

  if(paired==TRUE && missing(r)){
    error<-"Please specify r for paired designs."
    stop(error)
  }

  if(missing(r)){
    r<-0
  } else { r<-r}

  if(abs(r) > 1){
    error<-"Please double check. Correlation cannot surpass 1."
    stop(error)
  }

  if(missing(conf.int)){
    conf.int<-.9
  }

  if(missing(assur.lvl)){
    assur.lvl<-.99
  }

  if(assur.lvl < 0.8 || assur.lvl > 0.99){
    error<-"Assurance level must be between 0.8 - 0.99"
    stop(error)
  }

  zcrit<-abs(stats::qnorm((1-conf.int)/2))
  level<-paste(as.character(100*conf.int),"%",sep = "")
  level2<-paste(as.character(100*assur.lvl),"%",sep = "")

  if(paired==FALSE){
    n<-max(c(2*round((zcrit/moe)^2,0)),3)
    avg<-sqrt(2/n)*abs(stats::qt((1-conf.int)/2, n-1))
    df<-(2*ceiling(max(c(2*(zcrit/moe)^2),3)*(stats::qchisq(assur.lvl, max(c(2*(zcrit/moe)^2),3))/ceiling(max(c(2*(zcrit/moe)^2),3)))))
    ass<-(ceiling((2*(abs(stats::qt(((1-conf.int)/2),df))/moe)^2)*(stats::qchisq(assur.lvl, df)/df))+1)
    ass.avg<-sqrt(2/n)*abs(stats::qt((1-conf.int)/2, n-1))


    while(avg > moe){
      n<-n+1
      avg<-sqrt(2/n)*abs(stats::qt((1-conf.int)/2, n-1))
      if(avg < moe)
        break}

    while(ass.avg > moe){
      ass<-ass+1
      ass.avg<-sqrt(2/n)*abs(stats::qt((1-conf.int)/2, n-1))
      if(ass.avg < moe)
        break}
  }

  else{
    rmod<-2*(1-r)
    n<-max(c(round(rmod*(zcrit/moe)^2,0)),3)
    avg<-(1/sqrt(n))*abs(stats::qt((1-conf.int)/2, n-1))*sqrt(rmod)
    df<-round((n*stats::qchisq(assur.lvl,n)/n),0)
    ass<-ceiling((rmod*(abs(stats::qt(((1-conf.int)/2),df))/moe)^2)*(stats::qchisq(assur.lvl, df)/df))
    ass.avg<-(1/sqrt(ass))*abs(stats::qt((1-conf.int)/2, ass-1))*sqrt(rmod)

    while(avg > moe){
      n<-n+1
      avg<-(1/sqrt(n))*abs(stats::qt((1-conf.int)/2, n-1))*sqrt(rmod)
      if(avg < moe)
        break}

    while(ass.avg > moe){
      ass<-ass+1
      ass.avg<-(1/sqrt(ass))*abs(stats::qt((1-conf.int)/2, ass-1))*sqrt(rmod)
      if(ass.avg < moe)
        break}

  }

  if(paired==FALSE){
    table1<-matrix(c(moe,conf.int),ncol=1,byrow = T)
    rownames(table1)<-c("MOE","CI")
    colnames(table1)<-c(" ")
    cat("Independent Design")
    print(table1)
    cat("\nEstimated Sample Size (per group)\n","n = ",n,sep = "")
    cat("\n\nWith ",level2," Assurance\n","n = ",ass,sep = "")
  }

  else{
    table1<-matrix(c(moe,conf.int,r),ncol=1,byrow = T)
    rownames(table1)<-c("MOE","CI","r")
    colnames(table1)<-c(" ")
    cat("Paired Design")
    print(table1)
    cat("\nEstimated Total Sample Size\n","n = ",n,sep = "")
    cat("\n\nWith ",level2," Assurance\n","n = ",ass,sep = "")
  }
}
