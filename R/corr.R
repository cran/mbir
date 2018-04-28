#' Correlation Coefficient
#'
#'Provides magnitude-based inferences upon given \emph{r} value and sample size. Based upon WG Hopkins Microsoft Excel spreadsheet.
#'
#'@param r correlation coefficient
#'@param n sample size
#'@param conf.int (optional) confidence level of the interval. Defaults to \code{0.90}
#'@param swc (optional) number indicating smallest worthwhile change. Defaults to \code{0.1}
#'@param plot (optional) logical indicator specifying to print associated plot. Defaults to \code{FALSE}
#'@details Refer to vignette for further information.
#'@references Hopkins WG. (2007). A spreadsheet for deriving a confidence interval, mechanistic inference and clinical inference from a \emph{p} value. \emph{Sportscience} 11, 16-20. sportsci.org/2007/wghinf.htm
#'@examples corr(.40, 25, 0.95)
#'@export

corr <- function (r, n, conf.int=0.9, swc=0.1, plot=FALSE) {
  if (is.character(r) == TRUE || is.factor(r) == TRUE || is.character(n) ==
      TRUE || is.factor(n) == TRUE) {
    error <- "Sorry, data must be numeric or integer values."
    stop(error)
  }
  if (n < 4) {
    error <- "Sorry, not enough data."
    stop(error)
  }
  if (length(r) > 1) {
    error <- "Please enter only one effect size."
    stop(error)
  }

  if (abs(r) > 1) {
    error <- "Please double check. Correlation cannot surpass 1."
    stop(error)
  }

  if (abs(swc) >= 1) {
    error <- "Please double check. Smallest effect size of interest cannot surpass 1."
    stop(error)
  }

  if (swc <= 0 ) {
    error <- "Sorry, the smallest effect size of interest (swc) must be a positive number"
    stop(error)
  }

  #Z-transformation for the swc correlation
  threshold <- (0.5*log((1+swc)/(1-swc)))

  positive <- round(100 * (1 - stats::pnorm(threshold, mean = (0.5 *
                                                                 log((1 + r)/(1 - r))), sd = (1/sqrt(n - 3)))), digits = 1)
  negative <- round(100 * (stats::pnorm(-threshold, mean = (0.5 *
                                                              log((1 + r)/(1 - r))), sd = (1/sqrt(n - 3)))), digits = 1)
  trivial <- round(100 - positive - negative, digits = 1)

  LL <- (exp(2 * ((0.5 * log((1 + r)/(1 - r))) + (stats::qnorm(((100 -
                                                                   (100 * conf.int))/100/2))/sqrt(n - 3)))) - 1)/(exp(2 *
                                                                                                                        ((0.5 * log((1 + r)/(1 - r))) + (stats::qnorm(((100 -
                                                                                                                                                                          (100 * conf.int))/100/2))/sqrt(n - 3)))) + 1)

  UL <- (exp(2 * ((0.5 * log((1 + r)/(1 - r))) - (stats::qnorm(((100 -
                                                                   (100 * conf.int))/100/2))/sqrt(n - 3)))) - 1)/(exp(2 *
                                                                                                                        ((0.5 * log((1 + r)/(1 - r))) - (stats::qnorm(((100 -
                                                                                                                                                                          (100 * conf.int))/100/2))/sqrt(n - 3)))) + 1)
  cat("   Correlation Coefficient:\n")
  level <- paste(as.character(100 * conf.int), "%", sep = "")
  cat("   r = ", r, "\n", sep = "")
  cat("   n = ", n, "\n", sep = "")
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
                                                                                                                                 "Likely", ifelse(negative < 99, "Most Likely", ifelse(negative >=
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
  infer2 <- ifelse(r < 0, "Negative", "Positive")
  infer3 <- ifelse(infer == 1, lower, ifelse(infer == 2, trivial2,
                                             ifelse(infer == 3, higher)))
  mag <- ifelse(abs(r) < 0.1 || infer == 2, "Trivial", ifelse(abs(r) <
                                                                0.3, "Small", ifelse(abs(r) < 0.5, "Moderate", ifelse(abs(r) <
                                                                                                                        0.7, "Large", ifelse(abs(r) < 0.9, "Very Large", ifelse(abs(r) >=
                                                                                                                                                                                  0.9, "Very Large"))))))
  inference <- ifelse(abs(positive) >= 5 && abs(negative) > 5,
                      paste("Inference: Unclear Association."),
                      paste("Inference:", infer3, mag, infer2, "Correlation.",
                            sep = " "))
  cat(inference)

  #Creates plots of MBI *Note will not print if normal=FALSE
  if (plot == TRUE) {
    plot(NA, ylim = c(0, 1), xlim = c(min(LL, -swc) -
                                        max(UL - LL, swc - -swc)/10,
                                      max(UL, swc) + max(UL - LL, swc -
                                                             -swc)/10), bty = "l", yaxt = "n", ylab = "",
         xlab = "Correlation")
    graphics::points(x = r, y = 0.5, pch = 15, cex = 2)
    graphics::abline(v = swc, lty = 2)
    graphics::abline(v = -swc, lty = 2)
    graphics::abline(v = 0, lty = 2, col = "grey")
    graphics::segments(LL, 0.5, UL, 0.5, lwd = 3)
    graphics::title(main = paste(
      "r = ", round(r, digits = 3), " \n  ",
      100 * (conf.int), "% CI [", round(LL, digits = 3),
      ";", round(UL, digits = 3), "] ", " \n  ", inference,
      sep = ""), cex.main = 1)
  }

  rval <- list(r=r, r.LL=LL, r.UL=UL,
               mbiPositive=positive, mbiTrivial=trivial, mbiNegative=negative,
               inference=inference, swc=swc, conf.int=conf.int)
}
