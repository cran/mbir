#' Correlation Coefficient Test
#'
#'Provides magnitude-based inferences for the association between given data vectors. Evaluates normality assumption, performs either Pearson or Spearman correlation and subsequently estimates magnitude-based inferences.
#'
#'@param x,y numeric vectors of data values
#'@param conf.int (optional) confidence level of the interval. Defaults to \code{0.90}
#'@return Associated effect size measure, \emph{r}, and respective confidence intervals.
#'@details Refer to vignette for further information.
#'@examples a <- rnorm(25, 80, 35)
#'@examples b <- rnorm(25, 100, 35)
#'
#'@examples corr_test(a, b, 0.95)
#'@export

corr_test<-function(x,y,conf.int){

  if(length(x) != length(y) || sum(is.na(x)) > 0 || sum(is.na(y)) > 0){
    error<-"Sorry, data must be same length and complete cases."
    stop(error)
  }

  if(is.character(x) == TRUE || is.factor(x) == TRUE || is.character(y) == TRUE || is.factor(y) == TRUE){
    error<-"Sorry, data must be numeric or integer values."
    stop(error)
  }

  if(length(x) < 4 || length(y) < 4){
    error<-"Sorry, not enough data."
    stop(error)
  }

  if(missing(conf.int)){
    conf.int<-.9
  }

  x<-stats::na.omit(x)
  y<-stats::na.omit(y)
  full<-append(x,y)
  normal<-stats::shapiro.test(full)
  x.z<-(x-mean(x))/stats::sd(x)
  y.z<-(y-mean(y))/stats::sd(y)


  if(normal$p.value < .05 || max(x.z) > 3 || max(y.z) > 3){
    method<-"spearman"} else {method<-"pearson"}

  cor<-stats::cor.test(x,y,method = method,exact = F,na.action=na.omit,conf.level = conf.int)
  if(method == "spearman"){
    LL<-(exp(2*((.5*log((1+cor$estimate)/(1-cor$estimate)))+(stats::qnorm(((100-(100*conf.int))/100/2))/sqrt(length(x)-3))))-1)/(exp(2*((.5*log((1+cor$estimate)/(1-cor$estimate)))+(stats::qnorm(((100-(100*conf.int))/100/2))/sqrt(length(x)-3))))+1)
    UL<-(exp(2*((.5*log((1+cor$estimate)/(1-cor$estimate)))-(stats::qnorm(((100-(100*conf.int))/100/2))/sqrt(length(x)-3))))-1)/(exp(2*((.5*log((1+cor$estimate)/(1-cor$estimate)))-(stats::qnorm(((100-(100*conf.int))/100/2))/sqrt(length(x)-3))))+1)}

  # Based on Hopkins Magnitude-Based Inference
  positive<-round(100*(1-stats::pnorm(.1,mean = ((.5*log((1+cor$estimate)))/(1-cor$estimate)),sd = (1/sqrt(length(x)-3)))),digits = 1)
  negative<-round(100*(stats::pnorm(-.1,mean = ((.5*log((1+cor$estimate)))/(1-cor$estimate)),sd = (1/sqrt(length(x)-3)))),digits = 1)
  trivial<-round(100-positive-negative,digits = 1)

  ############################################

  level<-paste(as.character(100*conf.int),"%",sep = "")
  norm<-ifelse(method == "pearson","   Normality Observed, No Outliers Detected","   Skewness Observed or Outliers Detected")
  type<-ifelse(method == "pearson", "Pearson","Spearman")
  type2<-ifelse(method == "pearson","r = ","rho = ")
  cat(norm,"\n")
  cat("   Method: ",type,"\n\n",sep = " ")
  cat("   ",type2,round(cor$estimate,digits = 2),"\n",sep = "")

  if(method == "pearson"){cat("   ",level," CI ","[",round(cor$conf.int[1],digits = 2),", ",round(cor$conf.int[2],digits = 2),"]\n\n",sep = "")}
  else {cat("   ",level," CI ","[",round(LL,digits = 2),", ",round(UL,digits = 2),"]\n\n",sep = "")}

  table<-matrix(c("Negative","Trivial","Positive",negative,trivial,positive),nrow = 2,byrow = T)
  rownames(table)<-c(" ","MBI (%)")

  # Based on Hopkins Magnitude-Based Inference
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

  colnames(table)<-c(lower,trivial2,higher)
  title<-("   Magnitude-Based Inference")
  cat(title,"\n\n")
  print(table)
  cat("\n")
  infer<-which.max(table[2,])
  infer2<-ifelse(cor$estimate < 0,"Negative","Positive")
  infer3<-ifelse(infer == 1,lower,ifelse(infer == 2,trivial2,ifelse(infer == 3,higher)))
  mag<-ifelse(abs(cor$estimate) < 0.1 || infer == 2,"Trivial",ifelse(abs(cor$estimate) < 0.3, "Small",ifelse(abs(cor$estimate) < 0.5,"Moderate",ifelse(abs(cor$estimate) < 0.7,"Large",ifelse(abs(cor$estimate) < 0.9,"Very Large",ifelse(abs(cor$estimate) >= .9,"Very Large"))))))
  if(abs(positive) >= 5 && abs(negative) > 5){cat("Inference: Unclear Association.")}
  else {cat("Inference:",infer3,mag,infer2,"Correlation.",sep = " ")}
  r.stat<-cor$estimate
  if(method == "pearson"){
    r.LL<-cor$conf.int[1]
    r.UL<-cor$conf.int[2]
  } else {
      r.LL<-LL
      r.UL<-UL
      }
  rval<-list(r.stat=r.stat, r.LL=r.LL, r.UL=r.UL)


}
