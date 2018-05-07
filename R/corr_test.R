#' Correlation Coefficient Test
#'
#'Provides magnitude-based inferences for the association between given data vectors. Evaluates normality assumption, performs either Pearson or Spearman correlation and subsequently estimates magnitude-based inferences.
#'
#'@param x,y numeric vectors of data values
#'@param conf.int (optional) confidence level of the interval. Defaults to \code{0.90}
#'@param auto (character) logical indicator specifying if user wants function to programmatically detect statistical procedures. Defaults to \code{TRUE}
#'@param method (character) if \code{auto = F}, logical indicator specifying which correlation to execute (\code{pearson, spearman, kendall}). Defaults to \code{"pearson"}.
#'@param swc (optional) number indicating smallest worthwhile change. Defaults to \code{0.1}
#'@param plot (optional) logical indicator specifying to print associated plot. Defaults to \code{FALSE}
#'@return Associated effect size measure, \emph{r}, and respective confidence intervals.
#'@details Refer to vignette for further information.
#'@examples a <- rnorm(25, 80, 35)
#'@examples b <- rnorm(25, 100, 35)
#'
#'@examples corr_test(a, b, 0.95)
#'@export

corr_test <- function(x, y, conf.int=0.9, auto=TRUE, method="pearson", swc=0.1, plot=FALSE) {

  if (length(x) != length(y) || sum(is.na(x)) > 0 || sum(is.na(y)) > 0) {
    error <- "Sorry, data must be same length and complete cases."
    stop(error)
  }
  if (is.character(x) == TRUE || is.factor(x) == TRUE || is.character(y) == TRUE || is.factor(y) == TRUE) {
    error <- "Sorry, data must be numeric or integer values."
    stop(error)
  }
  if (length(x) < 4 || length(y) < 4) {
    error <- "Sorry, not enough data."
    stop(error)
  }
  if (missing(conf.int)) {
    conf.int <- 0.9
  }

  if (abs(swc) >= 1 ) {
    error <- "Sorry, the smallest effect size of interest (swc) must be less than 1"
    stop(error)
  }
  if (swc <= 0 ) {
    error <- "Sorry, the smallest effect size of interest (swc) must be a positive number"
    stop(error)
  }
  #Variance and normality checks were removed.
  x <- stats::na.omit(x)
  y <- stats::na.omit(y)
  full <- append(x, y)
  threshold <- (0.5*log((1+swc)/(1-swc)))

  #automated function
  if (auto==TRUE) {
    x.z <- (x - mean(x))/stats::sd(x)
    y.z <- (y - mean(y))/stats::sd(y)
    normal <- stats::shapiro.test(full)
    normlabel<-ifelse(normal$p.value<0.05 || max(x.z) > 3 || max(y.z) > 3,"   Normality Observed, No Outliers Detected","   Skewness Observed or Outliers Detected")
    method <- ifelse(normal$p.value<0.05 || max(x.z) > 3 || max(y.z) > 3, "spearman", "pearson")
  }

  #Pearson
  if (method == "pearson") {
    cor <- stats::cor.test(x, y, method = method, exact = F,
                           na.action = na.omit, conf.level = conf.int)
    corZ <- (0.5*log((1+cor$estimate)/(1-cor$estimate)))
  }


  if (method == "spearman") {
    cor <- stats::cor.test(x, y, method = method, exact = F,
                           na.action = na.omit, conf.level = conf.int)
    corZ <- (0.5*log((1+cor$estimate)/(1-cor$estimate)))

    LL <- (exp(2 * ((0.5 * log((1 + cor$estimate)/(1 - cor$estimate))) +
                      (stats::qnorm(((100 - (100 * conf.int))/100/2))/sqrt(length(x) -
                                                                             3)))) - 1)/(exp(2 * ((0.5 * log((1 + cor$estimate)/(1 -
                                                                                                                                   cor$estimate))) + (stats::qnorm(((100 - (100 * conf.int))/100/2))/sqrt(length(x) -
                                                                                                                                                                                                            3)))) + 1)
    UL <- (exp(2 * ((0.5 * log((1 + cor$estimate)/(1 - cor$estimate))) -
                      (stats::qnorm(((100 - (100 * conf.int))/100/2))/sqrt(length(x) -
                                                                             3)))) - 1)/(exp(2 * ((0.5 * log((1 + cor$estimate)/(1 -
                                                                                                                                   cor$estimate))) - (stats::qnorm(((100 - (100 * conf.int))/100/2))/sqrt(length(x) -
                                                                                                                                                                                                            3)))) + 1)
  }

  if (method == "kendall") {
    cor <- stats::cor.test(x, y, method = method, exact = F,
                           na.action = na.omit, conf.level = conf.int)
    corZ <- (0.5*log((1+cor$estimate)/(1-cor$estimate)))

    LL <- (exp(2 * ((0.5 * log((1 + cor$estimate)/(1 - cor$estimate))) +
                      (stats::qnorm(((100 - (100 * conf.int))/100/2))/sqrt(length(x) -
                                                                             3)))) - 1)/(exp(2 * ((0.5 * log((1 + cor$estimate)/(1 -
                                                                                                                                   cor$estimate))) + (stats::qnorm(((100 - (100 * conf.int))/100/2))/sqrt(length(x) -
                                                                                                                                                                                                            3)))) + 1)
    UL <- (exp(2 * ((0.5 * log((1 + cor$estimate)/(1 - cor$estimate))) -
                      (stats::qnorm(((100 - (100 * conf.int))/100/2))/sqrt(length(x) -
                                                                             3)))) - 1)/(exp(2 * ((0.5 * log((1 + cor$estimate)/(1 -
                                                                                                                                   cor$estimate))) - (stats::qnorm(((100 - (100 * conf.int))/100/2))/sqrt(length(x) -
                                                                                                                                                                                                            3)))) + 1)
  }
  #Originally an error see underlying code for  mbir::corr
  positive <- round(100 * (1 - stats::pnorm(threshold, mean = (0.5 *
                                                                 log((1 + cor$estimate)/(1 - cor$estimate))), sd = (1/sqrt(length(x) -
                                                                                                                             3)))), digits = 1)
  negative <- round(100 * (stats::pnorm(-threshold, mean = (0.5 *
                                                              log((1 + cor$estimate)/(1 - cor$estimate))), sd = (1/sqrt(length(x) -
                                                                                                                          3)))), digits = 1)
  trivial <- round(100 - positive - negative, digits = 1)
  level <- paste(as.character(100 * conf.int), "%", sep = "")
  type <- ifelse(method == "pearson", "Pearson", ifelse(method == "spearman", "Spearman", "Kendall"))
  type2 <- ifelse(method == "pearson", "r = ", "rho = ")
  if(auto==TRUE){cat(normlabel, "\n")}else{cat("\n")}
  cat("   Method: ", type, "\n\n", sep = " ")
  cat("   ", type2, round(cor$estimate, digits = 2), "\n",
      sep = "")
  if (method == "pearson") {
    cat("   ", level, " CI ", "[", round(cor$conf.int[1],
                                         digits = 2), ", ", round(cor$conf.int[2], digits = 2),
        "]\n\n", sep = "")
  }
  else {
    cat("   ", level, " CI ", "[", round(LL, digits = 2),
        ", ", round(UL, digits = 2), "]\n\n", sep = "")
  }
  table1 <- matrix(c("Negative", "Trivial", "Positive", negative,
                     trivial, positive), nrow = 2, byrow = T)
  rownames(table1) <- c(" ", "MBI (%)")
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
  colnames(table1) <- c(lower, trivial2, higher)
  title <- ("   Magnitude-Based Inference")
  cat(title, "\n\n")
  print(table1)
  cat("\n")
  infer <- which.max(table1[2, ])
  infer2 <- ifelse(cor$estimate < 0, "Negative", "Positive")
  infer3 <- ifelse(infer == 1, lower, ifelse(infer == 2, trivial2,
                                             ifelse(infer == 3, higher)))
  mag <- ifelse(abs(cor$estimate) < 0.1 || infer == 2, "Trivial",
                ifelse(abs(cor$estimate) < 0.3, "Small", ifelse(abs(cor$estimate) <
                                                                  0.5, "Moderate", ifelse(abs(cor$estimate) < 0.7,
                                                                                          "Large", ifelse(abs(cor$estimate) < 0.9, "Very Large",
                                                                                                          ifelse(abs(cor$estimate) >= 0.9, "Very Large"))))))
  #Inference code added to so the inference text can be saved
  inference <- ifelse(abs(positive) >= 5 && abs(negative) > 5,
                      paste("Inference: Unclear Association."),
                      paste("Inference:", infer3, mag, infer2, "Correlation.",
                            sep = " "))
  cat(inference)
  r.stat <- cor$estimate
  if (method == "pearson") {
    r.LL <- cor$conf.int[1]
    r.UL <- cor$conf.int[2]
  }
  else {
    r.LL <- LL
    r.UL <- UL
  }

  #Creates plots of MBI *Note will not print if normal=FALSE
  if (plot == TRUE) {
    plot(NA, ylim = c(0, 1), xlim = c(min(r.LL, -swc) -
                                        max(r.UL - r.LL, swc - -swc)/10,
                                      max(r.UL, swc) + max(r.UL - r.LL, swc -
                                                               -swc)/10), bty = "l", yaxt = "n", ylab = "",
         xlab = "Correlation")
    graphics::points(x = r.stat, y = 0.5, pch = 15, cex = 2)
    graphics::abline(v = swc, lty = 2)
    graphics::abline(v = -swc, lty = 2)
    graphics::abline(v = 0, lty = 2, col = "grey")
    graphics::segments(r.LL, 0.5, r.UL, 0.5, lwd = 3)
    graphics::title(main = paste(
      type2, round(r.stat, digits = 3), " \n  ",
      100 * (conf.int), "% CI [", round(r.LL, digits = 3),
      ";", round(r.UL, digits = 3), "] ", " \n  ", inference,
      sep = ""), cex.main = 1)
  }

  #Save List of values
  rval <- list(mean1 = round(mean(x, na.rm = T), digits = 3), sd1 = round(sd(x, na.rm = T), digits = 3),
               mean2 = round(mean(y, na.rm = T), digits = 3), sd2 = round(sd(y, na.rm = T), digits = 3),
               N = length(full), swc = swc,
               corr.stat = r.stat[[1]], z=corZ[[1]], z.LL = r.LL, z.UL = r.UL,
               norm=norm, type=type, inference=inference,
               mbiPositive=positive, mbiTrivial=trivial, mbiNegative=negative)
}
