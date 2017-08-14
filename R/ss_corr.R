#' Sample Size Estimation: Correlation Coefficient
#'
#'Estimates magnitude-based inferences upon planned sample size and \emph{r} value. Based upon WG Hopkins Microsoft Excel spreadsheet.
#'
#'@param n planned sample size
#'@param r planned correlation coefficient
#'@details Refer to vignette for further information.
#'@references Hopkins WG. (2006). Estimating sample size for magnitude-based inferences. \emph{Sportscience} 10, 63-70. sportsci.org/2006/wghss.htm
#'@examples ss_corr(n = 20, r = 0.2)
#'@export

ss_corr<-function(n, r){

  if(is.character(n) == TRUE || is.factor(n) == TRUE || is.character(r) == TRUE || is.factor(r) == TRUE){
    error<-"Sorry, data must be numeric or integer values."
    stop(error)
  }

  if(length(n) > 1){
    error<-"Please enter only one sample size."
    stop(error)
  }

  if(length(r) > 1){
    error<-"Please enter only one effect size."
    stop(error)
  }

  if(abs(r) > 1){
    error<-"Please double check. Correlation cannot surpass 1."
    stop(error)
  }

  positive<-round(100*stats::pnorm(1*(atanh(r)-atanh(.1))/sqrt(1/(n-3))),digits = 1)
  negative<-round(100*stats::pnorm(1*(atanh(-.1)-atanh(r))/sqrt(1/(n-3))),digits = 1)
  trivial<-round((100-positive-negative),digits = 1)
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
  cat("   n: ", n,"\n   ","r: ",r,"\n\n",sep = "")
  cat(title,"\n\n")
  print(table)
  cat("\n")
  infer<-which.max(table[2,])
  infer2<-ifelse(r < 0,"Negative","Positive")
  infer3<-ifelse(infer == 1,lower,ifelse(infer == 2,trivial2,ifelse(infer == 3,higher)))
  mag<-ifelse(abs(r) < 0.1 || infer == 2,"Trivial",ifelse(abs(r) < 0.3, "Small",ifelse(abs(r) < 0.5,"Moderate",ifelse(abs(r) < 0.7,"Large",ifelse(abs(r) < 0.9,"Very Large",ifelse(abs(r) >= .9,"Very Large"))))))
  if(abs(positive) >= 5 && abs(negative) > 5){cat("Inference: Unclear Association.")}
  else {cat("Inference:",infer3,mag,infer2,"Correlation.",sep = " ")}
}

