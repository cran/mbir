#'Test of Two Correlations
#'
#'Provides statistical inference upon the difference between two independent correlations.
#'
#'@param r1 correlation of group 1
#'@param n1 sample size of group 1
#'@param r2 correlation of group 2
#'@param n2 sample size of group 2
#'@param conf.int (optional) confidence level of the interval. Defaults to \code{0.90}
#'@param plot (optional) logical indicator specifying to print associated plot. Defaults to \code{FALSE}
#'@details Refer to vignette for further information.
#'@references Zou GY. (2007). Toward using confidence intervals to compare correlations. \emph{Psychological Methods}, 12, 399-413.
#'@examples corr_diff(r1 = 0.20, n1 = 71, r2 = 0.55, n2 = 46)
#'@export

corr_diff <- function (r1, n1, r2, n2, conf.int=0.9, plot=FALSE)
{
  if (is.character(r1) == TRUE || is.factor(r1) == TRUE ||
      is.character(n1) == TRUE || is.factor(n1) == TRUE) {
    error <- "Sorry, data must be numeric or integer values."
    stop(error)
  }
  if (is.character(r2) == TRUE || is.factor(r2) == TRUE ||
      is.character(n2) == TRUE || is.factor(n2) == TRUE) {
    error <- "Sorry, data must be numeric or integer values."
    stop(error)
  }
  if (length(r1) > 1 || length(n1) > 1 || length(r2) > 1 ||
      length(n2) > 1) {
    error <- "Please enter only one effect size."
    stop(error)
  }

  diff <- r2 - r1
  zcrit <- abs(stats::qnorm((1 - conf.int)/2))
  r1.z <- 0.5 * log((1 + r1)/(1 - r1))
  r1.sd <- 1/sqrt(n1 - 3)
  r1.ll <- r1.z - zcrit * r1.sd
  r1.ul <- r1.z + zcrit * r1.sd
  r2.z <- 0.5 * log((1 + r2)/(1 - r2))
  r2.sd <- 1/sqrt(n2 - 3)
  r2.ll <- r2.z - zcrit * r2.sd
  r2.ul <- r2.z + zcrit * r2.sd
  diff.UL <- diff + sqrt((r2.ul - r2)^2 + (r1 - r1.ll)^2)
  diff.LL <- diff - sqrt((r2 - r2.ll)^2 + (r1.ul - r1)^2)
  z.diff <- abs(r1.z - r2.z)
  z.diff.sd <- sqrt(1/(n1 - 3) + 1/(n2 - 3))
  z <- z.diff/z.diff.sd
  p <- 2 * (1 - stats::pnorm(z))
  dir <- ifelse(r2 > r1, ">", "<")
  level <- paste(as.character(100 * conf.int), "%", sep = "")
  cat("   Test of Two Correlations:\n")
  cat("   diff = ", diff, "\n", sep = "")
  cat("   ", level, " CI ", "[", round(diff.LL, digits = 2),
      ", ", round(diff.UL, digits = 2), "]\n", sep = "")
  cat("   p value = ", round(p, digits = 2), "\n\n", sep = "")
  #Save inference for list
  inference <- ifelse(diff.LL < 0 && diff.UL > 0,
                      paste("Inference: Lacking Evidence, r2 = r1, (CI contains 0).", sep = ""),
                      paste("Inference: Evidence Present, r2 ", dir, " r1, (CI does not contain 0).", sep = ""))

  #Creates plots of MBI *Note will not print if normal=FALSE
  if (plot == TRUE) {
    plot(NA, ylim = c(0, 1), xlim = c(diff.LL -
                                        (diff.UL - diff.LL)/10,
                                      (diff.UL) + (diff.UL - diff.LL)/10), bty = "l", yaxt = "n", ylab = "",
         xlab = "Difference in Correlations")
    graphics::points(x = diff, y = 0.5, pch = 15, cex = 2)
    graphics::abline(v = 0, lty = 2, col = "grey")
    graphics::segments(diff.LL, 0.5, diff.UL, 0.5, lwd = 3)
    graphics::title(main = paste(
      "difference = ", round(diff, digits = 2), " \n  ",
      100 * (conf.int), "% CI [", round(diff.LL, digits = 2),
      ";", round(diff.UL, digits = 2), "] ", " \n  ", inference,
      sep = ""), cex.main = 1)
  }


  #Save list of output
  rval <- list(diff=diff, diff.LL=diff.LL, diff.UL=diff.UL,
               n1=n1, r1=r1, r1.ll=r1.ll, r1.ul=r1.ul,
               n2=n2, r2=r2, r2.ll=r2.ll, r2.ul=r2.ul,
               p.value=p, inference=inference)
}

