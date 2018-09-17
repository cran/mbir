#' Effect Size Converter
#'
#'Converts between equivalent effect size measures: \emph{d}, \emph{r}, odds ratio.
#'
#'@param x numeric value
#'@param from (character) current effect size of \code{x}
#'@param to (character) effect size measure to convert to
#'@details Refer to vignette for further information.
#'@references Rosenthal R. (1994). Parametric measures of effect size. In H. Cooper & LV. Hedges (Eds.), \emph{The Handbook of Research Synthesis}. New York, NY: Sage.
#'@references Borenstein M, Hedges LV, Higgins JPT & Rothstein HR. (2009). \emph{Introduction to Meta-Analysis}. Chichester, West Sussex, UK: Wiley.
#'@examples
#'# Odds ratio to Cohen's d
#'es_convert(1.25, from = "or", to = "d")
#'@export

es_convert<-function(x,from = c("d","or","r"),to = c("d","or","r")){

  if(is.character(x) == TRUE || is.factor(x) == TRUE){
    error<-"Sorry, data must be numeric or integer values."
    stop(error)
  }

  if(missing(from) || missing(to)){
    error<-"Missing Argument: Please specify type of conversion.\n  Options: d-or,d-r,or-d,or-r,r-d,r-or"
    stop(error)
  }

  if(from == to){
    error<-"Incorrect Argument: Please specify type of conversion.\n  Options: d-or,d-r,or-d,or-r,r-d,r-or"
    stop(error)
  }

  if(from == "d" || from == "or" || from == "r"){
    correct1<-from
  } else {
    error<-"Incorrect Argument. Conversion types: d, or, r"
    stop(error)
  }

  if(to == "d" || to == "or" || to == "r"){
    correct2<-to
  } else {
    error<-"Incorrect Argument. Conversion types: d, or, r"
    stop(error)
  }

  if(from == "r" && abs(x) > 1){
    error<-"Please double check. Correlation cannot surpass 1."
    stop(error)
  }

  if(from == "or" && x <= 0){
    error<-"Please double check. Odds Ratio must be positive."
    stop(error)
  }

  ####################################################################

  if(correct1 == "d" && correct2 == "or"){
    convert<-round(exp(x*pi/sqrt(3)),digits = 2)
  }
  else if(correct1 == "d" && correct2 == "r"){
    convert<-round(sign(x)*sqrt((x^2)/((x^2)+4)),digits = 2)
  }
  else if(correct1 == "or" && correct2 == "d"){
    convert<-round(log(x)*(sqrt(3)/pi),digits = 2)
  }
  else if(correct1 == "or" && correct2 == "r"){
    d<-log(x)*(sqrt(3)/pi)
    convert<-round(sign(x)*sqrt((d^2)/((d^2)+4)),digits = 2)
  }
  else if(correct1 == "r" && correct2 == "d"){
    convert<-round((2*x)/sqrt(1-(x^2)),digits = 2)
  }
  else if(correct1 == "r" && correct2 == "or"){
    d<-(2*x)/sqrt(1-(x^2))
    convert<-round(exp((d*pi)/sqrt(3)),digits = 2)
  }

  cat("Effect Size Conversion\n")
  cat(from,"=",x,"\n",sep = " ")
  cat(to,"=", convert,sep = " ")

  rval<-list(convert=convert)
}
