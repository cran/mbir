#' Sample Size Estimation: Standardized Mean Difference
#'
#'Estimates magnitude-based inferences upon planned sample size and \emph{d} value. Based upon WG Hopkins Microsoft Excel spreadsheet.
#'
#'@param exp planned sample size of experimental group
#'@param con planned sample size of control group
#'@param es planned Cohen's \emph{d}
#'@details Refer to vignette for further information.
#'@references Hopkins WG. (2006). Estimating sample size for magnitude-based inferences. \emph{Sportscience} 10, 63-70. sportsci.org/2006/wghss.htm
#'@examples ss_smd(exp = 20, con = 15, es = 0.6)
#'@export

ss_smd<-function(exp, con, es){

  if(is.character(exp) == TRUE || is.factor(exp) == TRUE || is.character(con) == TRUE || is.factor(con) == TRUE || is.character(es) == TRUE || is.factor(es) == TRUE){
    error<-"Sorry, data must be numeric or integer values."
    stop(error)
  }

  if(length(es) > 1){
    error<-"Please enter only one effect size."
    stop(error)
  }

  n<-(exp+con)
  percent<-100*(exp/n)
  factor<-1/(percent/100*(1-percent/100))
  positive<-round(100-(100*stats::pt(-abs((es-.2)/sqrt(factor/n)),(n-2))),digits=1)
  negative<-round(100*(1-stats::pt(abs((es+.2)/sqrt(factor/n)),(n-2))),digits = 1)
  trivial<-round(100-positive-negative,digits = 1)

  cat("   n = ",n,"\n",sep = "")
  cat("   es = ",es,"\n\n",sep = "")
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
  title<-("   Estimated:\n   Magnitude-Based Inference")
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
