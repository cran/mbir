#' Standardized Mean Difference
#'
#'Provides magnitude-based inferences upon given \emph{d}, \emph{p}-value, and degrees of freedom. Based upon WG Hopkins Microsoft Excel spreadsheet.
#'
#'@param es effect size measure (Cohen's \emph{d})
#'@param p associated \emph{p}-value from t-statistic
#'@param df associated degrees of freedom from t-statistic
#'@param conf.int (optional) confidence level of the interval. Defaults to \code{0.90}
#'@param swc (optional) number indicating smallest worthwhile change. Defaults to \code{0.5}
#'@param plot (optional) logical indicator specifying to print associated plot. Defaults to \code{FALSE}
#'@details Refer to vignette for further information.
#'@references Hopkins WG. (2007). A spreadsheet for deriving a confidence interval, mechanistic inference and clinical inference from a \emph{p} value. \emph{Sportscience} 11, 16-20. sportsci.org/2007/wghinf.htm
#'@examples smd(.75, 0.06, 20, 0.95)
#'@export

smd <- function (es, p, df, conf.int=0.9, swc=0.5, plot=FALSE){

  if (is.character(es) == TRUE || is.factor(es) == TRUE ||
      is.character(p) == TRUE || is.factor(p) == TRUE || is.character(df) ==
      TRUE || is.factor(df) == TRUE) {
    error <- "Sorry, data must be numeric or integer values."
    stop(error)
  }
  if (length(es) > 1) {
    error <- "Please enter only one effect size."
    stop(error)
  }
  #Returns error if p value is entered incorrectly
  if (p >= 1 || p<= 0 ) {
    error <- "Please enter the exact p-value in decimal form (e.g., p=0.049)"
    stop(error)
  }

  ###error messages with swc is incorrectly entered
  if (swc <= 0 ) {
    error <- "Sorry, the smallest effect size of interest (swc) must be a positive number"
    stop(error)
  }

  negative <- round(100 * (ifelse((es - -swc) > 0, stats::pt((es -
                                                                -swc)/abs(es) * abs(stats::qt(p/2, df)), df, lower.tail = F),
                                  (1 - stats::pt((-swc - es)/abs(es) * abs(stats::qt(p/2,
                                                                                     df)), df, lower.tail = F)))), digits = 1)
  positive <- round(100 * (ifelse((es - swc) > 0, (1 - stats::pt((es -
                                                                    swc)/abs(es) * abs(stats::qt(p/2, df)), df, lower.tail = F)),
                                  stats::pt((swc - es)/abs(es) * abs(stats::qt(p/2, df)),
                                            df, lower.tail = F))), digits = 1)
  trivial <- round((100 - positive - negative), digits = 1)
  LL <- es - (stats::qt(((100 - (100 * conf.int))/100)/2, df)) *
    abs(es)/stats::qt(p/2, df)
  UL <- es + (stats::qt(((100 - (100 * conf.int))/100)/2, df)) *
    abs(es)/stats::qt(p/2, df)
  cat("   Standardized Mean Difference:\n")
  level <- paste(as.character(100 * conf.int), "%", sep = "")
  cat("   es = ", es, "\n", sep = "")
  cat("   p value = ", p, "\n", sep = "")
  cat("   ", level, " CI ", "[", round(LL, digits = 2), ", ",
      round(UL, digits = 2), "]\n\n", sep = "")
  table <- matrix(c("Negative", "Trivial", "Positive", negative,
                    trivial, positive), nrow = 2, byrow = T)
  rownames(table) <- c(" ", "MBI (%)")
  lower <- ifelse(negative < swc, "Most Unlikely", ifelse(negative <
                                                            5, "Very Unlikely", ifelse(negative < 25, "Unlikely",
                                                                                       ifelse(negative < 75, "Possibly", ifelse(negative < 95,
                                                                                                                                "Likely", ifelse(negative < 99, "Most Likely", ifelse(negative >=
                                                                                                                                                                                        99, "Almost Certainly")))))))
  trivial2 <- ifelse(trivial < swc, "Most Unlikely", ifelse(trivial <
                                                              5, "Very Unlikely", ifelse(trivial < 25, "Unlikely",
                                                                                         ifelse(trivial < 75, "Possibly", ifelse(trivial < 95,
                                                                                                                                 "Likely", ifelse(negative < 99, "Most Likely", ifelse(negative >=
                                                                                                                                                                                         99, "Almost Certainly")))))))
  higher <- ifelse(positive < swc, "Most Unlikely", ifelse(positive <
                                                             5, "Very Unlikely", ifelse(positive < 25, "Unlikely",
                                                                                        ifelse(positive < 75, "Possibly", ifelse(positive < 95,
                                                                                                                                 "Likely", ifelse(negative < 99, "Most Likely", ifelse(negative >=
                                                                                                                                                                                         99, "Almost Certainly")))))))
  colnames(table) <- c(lower, trivial2, higher)
  title <- ("   Magnitude-Based Inference")
  cat(title, "\n\n")
  print(table)
  cat("\n")
  infer <- which.max(table[2, ])
  infer2 <- ifelse(infer == 1, lower, ifelse(infer == 2, trivial2,
                                             ifelse(infer == 3, higher)))
  mag <- ifelse(abs(swc) < 0.2 || infer == 2, "Trivial", ifelse(abs(swc) <
                                                                  0.6, "Small", ifelse(abs(swc) < 1.2, "Moderate", ifelse(abs(swc) <
                                                                                                                            2, "Large", ifelse(abs(swc) >= 2, "Very Large")))))
  dir <- ifelse(infer == 1, "Decrease.", ifelse(infer == 2,
                                                "Difference.", ifelse(infer == 3, "Increase.")))
  inference <- ifelse(abs(positive) >= 5 && abs(negative) > 5,
                      paste("Inference: Unclear Difference."),
                      paste("Inference:", infer2, mag, dir, sep = " "))
  cat(inference)

  if (plot == TRUE) {
    plot(NA, ylim = c(0, 1), xlim = c(min(LL, -swc) -
                                        max(UL - LL, swc - -swc)/10,
                                      max(UL, swc) + max(UL - LL, swc -
                                                           -swc)/10), bty = "l", yaxt = "n", ylab = "",
         xlab = "Effect Size")
    graphics::points(x = es, y = 0.5, pch = 15, cex = 2)
    graphics::abline(v = swc, lty = 2)
    graphics::abline(v = -swc, lty = 2)
    graphics::abline(v = 0, lty = 2, col = "grey")
    graphics::segments(LL, 0.5, UL, 0.5, lwd = 3)
    graphics::title(main = paste(
      "Cohen's d = ", round(es, digits = 3), " \n  ",
      100 * (conf.int), "% CI [", round(LL, digits = 3),
      ";", round(UL, digits = 3), "] ", " \n  ", "Inference: ", infer2, " ", mag," ", dir,
      sep = ""), cex.main = 1)
  }

  rval <- list(es=es, es.LL=LL, es.UL=UL,
               p.value=p, conf.int=conf.int, swc=swc,
               mbiPositive=positive, mbiTrivial=trivial, mbiNegative=negative, inference=inference)
}


