#' Odds Ratio
#'
#'Provides magnitude-based inferences upon given odds ratio and \emph{p}-value. Based upon WG Hopkins Microsoft Excel spreadsheet.
#'
#'@param or odds ratio
#'@param p associated \emph{p}-value
#'@param conf.int (optional) confidence level of the interval. Defaults to \code{0.90}
#'@details Refer to vignette for further information.
#'@references Hopkins WG. (2007). A spreadsheet for deriving a confidence interval, mechanistic inference and clinical inference from a \emph{p} value. \emph{Sportscience} 11, 16-20. sportsci.org/2007/wghinf.htm
#'@examples odds(1.25, 0.06, 0.95)
#'@export

odds <- function (or, p, conf.int=0.9)
{
  if (is.character(or) == TRUE || is.factor(or) == TRUE ||
      is.character(p) == TRUE || is.factor(p) == TRUE) {
    error <- "Sorry, data must be numeric or integer values."
    stop(error)
  }
  if (length(or) > 1) {
    error <- "Please enter only one effect size."
    stop(error)
  }
  #New addition here; may avoid more errors by users in the long run
  if (conf.int>= 1 || conf.int <=0) {
    error <- "Confidence interval preference must be between 0 and 1"
    stop(error)
  }

  negative <- round(100 * (stats::pnorm((log(0.9) - log(or))/abs(log(or)/stats::qnorm(p/2)))),
                    digits = 1)
  positive <- round(100 * (1 - stats::pnorm((log(1.11) - log(or))/abs(log(or)/stats::qnorm(p/2)))),
                    digits = 1)
  trivial <- round((100 - positive - negative), digits = 1)
  LL <- exp(log(or) + stats::qnorm((100 - (100 * conf.int))/100/2) *
              abs(log(or)/stats::qnorm(p/2)))
  UL <- exp(log(or) - stats::qnorm((100 - (100 * conf.int))/100/2) *
              abs(log(or)/stats::qnorm(p/2)))
  cat("   Odds Ratio:\n")
  level <- paste(as.character(100 * conf.int), "%", sep = "")
  cat("   or = ", or, "\n", sep = "")
  cat("   p value = ", p, "\n", sep = "")
  cat("   ", level, " CI ", "[", round(LL, digits = 2), ", ",
      round(UL, digits = 2), "]\n\n", sep = "")
  table <- matrix(c("Negative", "Trivial", "Positive", negative,
                    trivial, positive), nrow = 2, byrow = T)
  rownames(table) <- c(" ", "MBI (%)")
  lower <- ifelse(negative < 0.5, "Most Unlikely", ifelse(negative <
                                                            5, "Very Unlikely", ifelse(negative < 25, "Unlikely",
                                                                                       ifelse(negative < 75, "Possibly", ifelse(negative < 95,
                                                                                                                                "Likely", ifelse(negative < 99, "Most Likely", ifelse(negative >=
                                                                                                                                                                                        99, "Almost Certainly")))))))
  trivial2 <- ifelse(trivial < 0.5, "Most Unlikely", ifelse(trivial <
                                                              5, "Very Unlikely", ifelse(trivial < 25, "Unlikely",
                                                                                         ifelse(trivial < 75, "Possibly", ifelse(trivial < 95,
                                                                                                                                 "Likely", ifelse(trivial < 99, "Most Likely", ifelse(trivial >=
                                                                                                                                                                                        99, "Almost Certainly")))))))
  higher <- ifelse(positive < 0.5, "Most Unlikely", ifelse(positive <
                                                             5, "Very Unlikely", ifelse(positive < 25, "Unlikely",
                                                                                        ifelse(positive < 75, "Possibly", ifelse(positive < 95,
                                                                                                                                 "Likely", ifelse(positive < 99, "Most Likely", ifelse(positive >=
                                                                                                                                                                                         99, "Almost Certainly")))))))
  colnames(table) <- c(lower, trivial2, higher)
  title <- ("   Magnitude-Based Inference")
  cat(title, "\n\n")
  print(table)
  cat("\n")
  infer <- which.max(table[2, ])
  dir <- ifelse(or < 1, "Decrease", "Increase")
  infer2 <- ifelse(infer == 1, lower, ifelse(infer == 2, trivial2,
                                             ifelse(infer == 3, higher)))
  if (or < 1) {
    mag <- ifelse(or > (1/1.5) || infer == 2, "Trivial",
                  ifelse(or > (1/3.5), "Small", ifelse(or > (1/9),
                                                       "Moderate", ifelse(or > (1/32), "Large", ifelse(or <=
                                                                                                         (1/32), "Very Large")))))
  }
  else {
    mag <- ifelse(abs(or) < 1.5 || infer == 2, "Trivial",
                  ifelse(abs(or) < 3.5, "Small", ifelse(abs(or) < 9,
                                                        "Moderate", ifelse(abs(or) < 32, "Large", ifelse(abs(or) >=
                                                                                                           32, "Very Large")))))
  }

  inference <- ifelse(abs(positive) >= 5 && abs(negative) > 5,
                      paste("Inference: Unclear Difference."),
                      paste("Inference:", infer2, mag, dir, sep = " "))
  cat(inference)

  rval <- list(or=or, or.LL=round(LL,2),or.UL=round(UL,2),p.value=p,
               mbiPositive=positive, mbiTrivial=trivial, mbiNegative=negative,
               inference)

}
