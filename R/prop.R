#'Test of Two Proportions
#'
#'Provides magnitude-based inferences upon given proportions and sample sizes. Based upon WG Hopkins Microsoft Excel spreadsheet.
#'
#'@param p1 proportion of group 1
#'@param n1 sample size of group 1
#'@param p2 proportion of group 2
#'@param n2 sample size of group 2
#'@param conf.int (optional) confidence level of the interval. Defaults to \code{0.90}
#'@details Refer to vignette for further information.
#'@references Hopkins WG. (2007). A spreadsheet for deriving a confidence interval, mechanistic inference and clinical inference from a \emph{p} value. \emph{Sportscience} 11, 16-20. sportsci.org/2007/wghinf.htm
#'@examples prop(p1 = 0.7, n1 = 25, p2 = 0.5, n2 = 20)
#'@export



prop<-function(p1, n1, p2, n2, conf.int){

  if(is.character(p1) == TRUE || is.factor(p1) == TRUE || is.character(n1) == TRUE || is.factor(n1) == TRUE){
    error<-"Sorry, data must be numeric or integer values."
    stop(error)
  }

  if(is.character(p2) == TRUE || is.factor(p2) == TRUE || is.character(n2) == TRUE || is.factor(n2) == TRUE){
    error<-"Sorry, data must be numeric or integer values."
    stop(error)
  }

  if(length(p1) > 1 || length(n1) > 1 || length(p2) > 1 || length(n2) > 1){
    error<-"Please enter only one effect size."
    stop(error)
  }

  if(missing(conf.int)){
    conf.int<-.9
  }

  diff<-p2-p1
  zcrit<-abs(stats::qnorm((1-conf.int)/2))
  n<-n1+n2
  x1<-p1*n1
  x2<-p2*n2
  a<-2*x1+zcrit^2
  b<-zcrit*sqrt((zcrit^2)+4*x1*(1-p1))
  c<-2*(n1+zcrit^2)
  p1.LL<-(a-b)/c
  p1.UL<-(a+b)/c
  a2<-2*x2+zcrit^2
  b2<-zcrit*sqrt((zcrit^2)+4*x2*(1-p2))
  c2<-2*(n2+zcrit^2)
  p2.LL<-(a2-b2)/c2
  p2.UL<-(a2+b2)/c2
  diff.LL<-diff-sqrt((p2-p2.LL)^2+(p1.UL-p1)^2)
  diff.UL<-diff+sqrt((p1-p1.LL)^2+(p2.UL-p2)^2)
  nx1<-n1-x1
  nx2<-n2-x2
  tab<-as.table(rbind(c(x1,x2),c(nx1,nx2)))
  chi<-stats::chisq.test(tab,correct = F,simulate.p.value = T)
  phi<-sqrt(chi$statistic/n)

  # Based on Hopkins Magnitude-Based Inference
  positive<-round(100*(1-stats::pnorm(.1,mean = (.5*log((1+phi)/(1-phi))),sd = (1/sqrt(n-3)))),digits = 1)
  negative<-round(100*(stats::pnorm(-.1,mean = (.5*log((1+phi)/(1-phi))),sd = (1/sqrt(n-3)))),digits = 1)
  trivial<-round(100-positive-negative,digits = 1)

  lower<-ifelse(negative<.5,"Most Unlikely",
                ifelse(negative<5,"Very Unlikely",
                       ifelse(negative<25,"Unlikely",
                              ifelse(negative<75,"Possibly",
                                     ifelse(negative<95,"Likely",
                                            ifelse(negative<99,"Most Likely",
                                                   ifelse(negative>=99,"Almost Certainly")))))))
  trivial2<-ifelse(trivial<.5,"Most Unlikely",
                   ifelse(trivial<5,"Very Unlikely",
                          ifelse(trivial<25,"Unlikely",
                                 ifelse(trivial<75,"Possibly",
                                        ifelse(trivial<95,"Likely",
                                               ifelse(negative<99,"Most Likely",
                                                      ifelse(negative>=99,"Almost Certainly")))))))
  higher<-ifelse(positive<.5,"Most Unlikely",
                 ifelse(positive<5,"Very Unlikely",
                        ifelse(positive<25,"Unlikely",
                               ifelse(positive<75,"Possibly",
                                      ifelse(positive<95,"Likely",
                                             ifelse(negative<99,"Most Likely",
                                                    ifelse(negative>=99,"Almost Certainly")))))))
  ############################################

  level<-paste(as.character(100*conf.int),"%",sep = "")
  trivadd<-negative+trivial
  cat("   Test of Two Proportions:\n")
  cat("   diff = ",round(diff,2),"\n",sep = "")
  cat("   chi-sq = ",unname(round(chi$statistic,digits = 2)),"\n",sep = "")
  cat("   ",level," CI ","[",round(diff.LL,digits = 2),", ",round(diff.UL,digits = 2),"]\n",sep = "")
  cat("   p value = ",unname(round(chi$p.value,digits = 2)),"\n\n",sep = "")
  cat("   phi = ",unname(round(phi,digits = 2)),"\n",sep = "")
  table<-matrix(c("Trivial","Substantial",trivadd,positive),nrow = 2,byrow = T)
  rownames(table)<-c(" ","MBI (%)")
  colnames(table)<-c(trivial2,higher)
  title<-("\n   Magnitude-Based Inference")
  cat(title,"\n\n")
  print(table)
  cat("\n")
  infer<-which.max(table[2,])
  infer3<-ifelse(infer == 1,lower,higher)
  mag<-ifelse(abs(phi) < 0.1,"Trivial",ifelse(abs(phi) < 0.3, "Small",ifelse(abs(phi) < 0.5,"Moderate",ifelse(abs(phi) < 0.7,"Large",ifelse(abs(phi) < 0.9,"Very Large",ifelse(abs(phi) >= .9,"Very Large"))))))
  if(abs(positive) >= 5 && abs(negative) > 5 || infer == 1){cat("Inference: Unclear Effect Size.")}
  else {cat("Inference:",infer3,mag,"Effect Size.",sep = " ")}

}


