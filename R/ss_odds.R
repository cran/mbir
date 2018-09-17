#' Sample Size Estimation: Odds Ratio
#'
#'Estimates magnitude-based inferences upon planned sample size and odds ratio. Based upon WG Hopkins Microsoft Excel spreadsheet.
#'
#'@param exp planned sample size of experimental group
#'@param con planned sample size of control group
#'@param or planned odds ratio
#'@details Refer to vignette for further information.
#'@references Hopkins WG. (2006). Estimating sample size for magnitude-based inferences. \emph{Sportscience} 10, 63-70. sportsci.org/2006/wghss.htm
#'@examples ss_odds(exp = 15, con = 18, or = 3.25)
#'@export

ss_odds<-function(exp, con, or){

  if(is.character(exp) == TRUE || is.factor(exp) == TRUE || is.character(con) == TRUE || is.factor(con) == TRUE || is.character(or) == TRUE || is.factor(or) == TRUE){
    error<-"Sorry, data must be numeric or integer values."
    stop(error)
  }

  if(length(or) > 1){
    error<-"Please enter only one effect size."
    stop(error)
  }

  n<-(exp+con)
  percent<-100*(exp/n)
  factor<-(1/(percent/100))*(2/(1.11*50/(100-50))+1+1.11*50/(100-50))+(1/(1-percent/100))*(1/(50/100)+1/(1-50/100))
  positive<-round(100*(stats::pnorm((log(or)-log(1.11))/sqrt(factor/n))), digits = 1)
  negative<-round(100*(1-stats::pnorm((log(or)-log(0.9))/sqrt(factor/n))), digits = 1)
  trivial<-round((100-positive-negative), digits = 1)

  cat("   n = ",n,"\n",sep = "")
  cat("   OR = ",or,"\n\n",sep = "")
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
  dir<-ifelse(or < 1,"Decrease","Increase")
  infer2<-ifelse(infer == 1,lower,ifelse(infer == 2,trivial2,ifelse(infer == 3,higher)))
  if(or < 1){mag<-ifelse(or > (1/1.5) || infer == 2,"Trivial",ifelse(or > (1/3.5),"Small",ifelse(or > (1/9),"Moderate",ifelse(or > (1/32),"Large",ifelse(or <= (1/32),"Very Large")))))}
  else{mag<-ifelse(abs(or) < 1.5 || infer == 2,"Trivial",ifelse(abs(or) < 3.5, "Small",ifelse(abs(or) < 9,"Moderate",ifelse(abs(or) < 32,"Large",ifelse(abs(or) >= 32,"Very Large")))))}
  if(abs(positive) >= 5 && abs(negative) > 5){cat("Inference: Unclear Difference.")}
  else {cat("Inference:",infer2,mag,dir,sep = " ")}

  Inference <-  ifelse(abs(positive) >= 5 && abs(negative) > 5,
                       paste("Inference: Unclear Difference."),
                       paste("Inference:", infer2, mag, dir, sep = " "))

  invisible(list(
    mbiNegative = negative,
    mbiTrivial = trivial,
    mbiPositive = positive,
    Inference = Inference)
  )
}

