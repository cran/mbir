#' Standardized Mean Difference
#'
#'Provides magnitude-based inferences upon given \emph{d}, \emph{p}-value, and degrees of freedom. Based upon WG Hopkins Microsoft Excel spreadsheet.
#'
#'@param es effect size measure (Cohen's \emph{d})
#'@param p associated \emph{p}-value from t-statistic
#'@param df associated degrees of freedom from t-statistic
#'@param conf.int (optional) confidence level of the interval. Defaults to \code{0.90}
#'@details Refer to vignette for further information.
#'@references Hopkins WG. (2007). A spreadsheet for deriving a confidence interval, mechanistic inference and clinical inference from a \emph{p} value. \emph{Sportscience} 11, 16-20. sportsci.org/2007/wghinf.htm
#'@examples smd(.75, 0.06, 20, 0.95)
#'@export

smd<-function(es,p,df,conf.int){

  if(is.character(es) == TRUE || is.factor(es) == TRUE || is.character(p) == TRUE || is.factor(p) == TRUE || is.character(df) == TRUE || is.factor(df) == TRUE){
    error<-"Sorry, data must be numeric or integer values."
    stop(error)
  }

  if(length(es) > 1){
    error<-"Please enter only one effect size."
    stop(error)
  }

  if(missing(conf.int)){
    conf.int<-.9
  }

  # Based on Hopkins Magnitude-Based Inference
  negative<-round(100*(ifelse((es--.5) > 0,stats::pt((es--.5)/abs(es)*abs(stats::qt(p/2,df)),df,lower.tail=F),(1-stats::pt((-.5-es)/abs(es)*abs(stats::qt(p/2,df)),df,lower.tail=F)))),digits = 1)
  positive<-round(100*(ifelse((es-.5) > 0,(1-stats::pt((es-.5)/abs(es)*abs(stats::qt(p/2,df)),df,lower.tail = F)),stats::pt((.5-es)/abs(es)*abs(stats::qt(p/2,df)),df,lower.tail = F))),digits = 1)
  trivial<-round((100-positive-negative),digits = 1)
  LL<-es-(stats::qt(((100-(100*conf.int))/100)/2,df))*abs(es)/stats::qt(p/2,df)
  UL<-es+(stats::qt(((100-(100*conf.int))/100)/2,df))*abs(es)/stats::qt(p/2,df)

  ############################################

  cat("   Standardized Mean Difference:\n")
  level<-paste(as.character(100*conf.int),"%",sep = "")
  cat("   es = ",es,"\n",sep = "")
  cat("   p value = ",p,"\n",sep = "")
  cat("   ",level," CI ","[",round(LL,digits = 2),", ",round(UL,digits = 2),"]\n\n",sep = "")
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
  infer2<-ifelse(infer == 1,lower,ifelse(infer == 2,trivial2,ifelse(infer == 3,higher)))
  mag<-ifelse(abs(es) < 0.2 || infer == 2,"Trivial",ifelse(abs(es) < 0.6, "Small",ifelse(abs(es) < 1.2,"Moderate",ifelse(abs(es) < 2.0,"Large",ifelse(abs(es) >= 2.0,"Very Large")))))
  dir<-ifelse(infer == 1,"Decrease.",ifelse(infer == 2, "Difference.",ifelse(infer == 3,"Increase.")))
  if(abs(positive) >= 5 && abs(negative) > 5){cat("Inference: Unclear Difference.")}
  else {cat("Inference:",infer2,mag,dir,sep = " ")}

}
