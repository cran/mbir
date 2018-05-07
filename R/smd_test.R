#' Standardized Mean Difference Test
#'
#'Performs two-sample difference of means analysis to produce magnitude-based inferences. Evaluates both normality and homogeneity, performs either t-test or wilcoxon test, computes effect sizes and estimates magnitude-based inferences. Allows both independent and paired designs.
#'
#'@param x,y numeric vectors of data values
#'@param paired (character) logical indicator specifying if \code{x} and \code{y} are paired \code{(TRUE)} or independent \code{(FALSE)}
#'@param auto (character) logical indicator specifying if user wants function to programmatically detect statistical procedures. Defaults to \code{TRUE}
#'@param var (optional) if \code{auto = F}, logical indicator specifying if homogeneity of variance assumed. Defaults to \code{TRUE}
#'@param normal (optional) if \code{auto = F}, logical indicator specifying if normality assumed. Defaults to \code{TRUE}
#'@param conf.int (optional) confidence level of the interval. Defaults to \code{0.90}
#'@param mu (optional) number indicating true difference in means to test against. Defaults to zero.
#'@param swc (optional) number indicating smallest worthwhile change. Defaults to \code{0.5}
#'@param plot (optional) logical indicator specifying to print associated plot. Defaults to \code{FALSE}
#'@return Associated effect size measures (\emph{d}, \emph{r}, odds ratio) and respective confidence intervals based upon which statistical test(s) performed.
#'@details Refer to vignette for further information.
#'@examples a <- rnorm(25, 80, 35)
#'@examples b <- rnorm(25, 100, 50)
#'
#'@examples smd_test(a, b, paired = FALSE, conf.int=0.95)
#'@export

smd_test<-function(x, y, paired = c(TRUE, FALSE), auto=TRUE, var = TRUE, normal = TRUE, conf.int = 0.9, mu = 0, swc = 0.5, plot=FALSE) {

  if (is.character(x) == TRUE ||
      is.factor(x) == TRUE || is.character(y) ==
      TRUE || is.factor(y) == TRUE) {
    error <- "Sorry, data must be numeric or integer values."
    stop(error)
  }
  if (length(x) < 3 || length(y) < 3) {
    error <- "Sorry, not enough data."
    stop(error)
  }
  if (missing(paired)) {
    error <-
      "Missing Argument: Please specify if data is paired = TRUE/FALSE."
    stop(error)
  }

  if (swc <= 0 ) {
    error <- "Sorry, the smallest effect size of interest be a positive number"
    stop(error)
  }

  if (length(x) != length(y)) {
    max.len <- max(length(x), length(y))
    x <- c(x, rep(NA, max.len - length(x)))
    y <- c(y, rep(NA, max.len - length(y)))
  }

  if (auto==TRUE){
    normal.x<-stats::shapiro.test(x)
    normal.y<-stats::shapiro.test(y)
    Normal<-ifelse(normal.x$p.value < .05 || normal.y$p.value < .05, FALSE, TRUE)
    normlabel<-ifelse(normal.x$p.value < .05 || normal.y$p.value < .05,"   Skewness Observed","   Normality Observed")
    equal<-stats::var.test(x,y)
    variance<-ifelse(equal$p.value < .05, FALSE,TRUE)
    variance2<-ifelse(equal$p.value < .05, "Unequal Variance","Equal Variance")
  }

  if(auto==FALSE){
    NormTest <- FALSE
    VarTest <- FALSE
    Normal <- normal
    equal <- var
    variance <- ifelse(var == FALSE, FALSE, TRUE)
    variance2 <- ifelse(var == FALSE, "Unequal Variance","Equal Variance")
    normlabel<-ifelse(normal.x$p.value < .05 || normal.y$p.value < .05,"   Skewness Assumed","   Normality Assumed")
  }


  length.x <- length(x) - sum(is.na(x))
  length.y <- length(y) - sum(is.na(y))
  n2 <- length.x + length.y
  n <- length.x
  ind.stdr <- sqrt(1 / length.x + 1 / length.y)
  pair.stdr <- sqrt((sd(x, na.rm = T) ^ 2 + sd(y, na.rm = T) ^ 2) / 2)

  if (Normal == FALSE) {
    if (0 %in% x || 0 %in% y) {
      test <- stats::t.test(
        log(abs(x + 1)),
        log(abs(y +
                  1)),
        var.equal = variance,
        paired = paired,
        conf.level = conf.int,
        na.action = na.omit,
        mu = mu
      )
    }
    else {
      test <-
        stats::t.test(
          log(abs(x)),
          log(abs(y)),
          var.equal = variance,
          paired = paired,
          conf.level = conf.int,
          na.action = na.omit,
          mu = mu
        )
    }

    estimate <- ifelse(
      paired == T,
      unname(test$estimate),
      (unname(test$estimate[1]) - unname(test$estimate[2]))
    )

    OR <- round(exp(estimate), digits = 2)
    LL <- round(exp(test$conf.int[1]), digits = 2)
    UL <- round(exp(test$conf.int[2]), digits = 2)



    rank <- stats::wilcox.test(
      x,
      y,
      paired = paired,
      conf.int = T,
      conf.level = 0.9,
      na.action = na.omit,
      correct = F,
      exact = F
    )
    suppressWarnings(warning(rank))
    r <- ifelse(OR < 1,
                stats::qnorm(rank$p.value / 2) / sqrt(n2),
                abs(stats::qnorm(rank$p.value / 2) / sqrt(n2)))
    r.LL <-
      (exp(2 * ((swc * log((1 + r) / (1 - r)
      )) + (
        stats::qnorm(((
          100 -
            (100 * conf.int)
        ) / 100 / 2)) / sqrt(n - 3)
      ))) - 1) / (exp(2 *
                        ((swc * log((1 + r) /
                                        (1 - r)
                        )) + (
                          stats::qnorm(((
                            100 -
                              (100 * conf.int)
                          ) / 100 / 2)) / sqrt(n - 3)
                        ))) + 1)
    r.UL <-
      (exp(2 * ((swc * log((1 + r) / (1 - r)
      )) - (
        stats::qnorm(((
          100 -
            (100 * conf.int)
        ) / 100 / 2)) / sqrt(n - 3)
      ))) - 1) / (exp(2 *
                        ((swc * log((1 + r) /
                                        (1 - r)
                        )) - (
                          stats::qnorm(((
                            100 -
                              (100 * conf.int)
                          ) / 100 / 2)) / sqrt(n - 3)
                        ))) + 1)
    level <- paste(as.character(100 * conf.int), "%", sep = "")
    negative <-
      round(100 * (stats::pnorm((log(0.9) - log(OR)) / abs(log(OR) / stats::qnorm(test$p.value /
                                                                                    2))
      )),
      digits = 1)
    positive <- round(100 * (1 - stats::pnorm((log(1.11) -
                                                 log(OR)) / abs(log(OR) /
                                                                  stats::qnorm(test$p.value / 2))
    )),
    digits = 1)
    trivial <- round((100 - positive - negative), digits = 1)
    table1 <- matrix(c(
      round(test$statistic, digits = 5),
      round(test$parameter, digits = 5),
      round(test$p.value,
            digits = 5)
    ),
    ncol = 1,
    byrow = T)
    rownames(table1) <- c("t", "df", "P value")
    colnames(table1) <- c("Log-Transformed")
    table2 <- matrix(c(
      round(rank$statistic, digits = 3),
      round(rank$p.value, digits = 3),
      round(rank$conf.int[1],
            digits = 3),
      round(rank$conf.int[2], digits = 3)
    ),
    ncol = 1,
    byrow = T)
    if (paired == TRUE) {
      rownames(table2) <- c("V",
                            "P value",
                            paste(level,
                                  "CI LL", sep = " "),
                            paste(level, "CI UL", sep = " "))
    }
    else {
      rownames(table2) <- c("W",
                            "P value",
                            paste(level,
                                  "CI LL", sep = " "),
                            paste(level, "CI UL", sep = " "))
    }
    colnames(table2) <- c("Rank-Based")
    title <- paste(normlabel, variance2, sep = ", ")
    title2 <- paste("   Method:", test$method, sep = " ")
    title3 <- paste("   Method:", rank$method, sep = " ")
    cat("   Mean of ", deparse(substitute(x)), " = ", round(mean(x,
                                                                 na.rm = T), digits = 2), "; ", "Mean of ", deparse(substitute(y)),
        " = ", round(mean(y, na.rm = T), digits = 2), "\n",
        sep = "")
    cat(title, title2, sep = "\n")
    cat("\n")
    print(table1)
    cat("\n")
    cat("   ", "Factor", " = ", round(OR, digits = 2), "\n",
        sep = "")
    cat("   ", level, " CI ", "[", round(LL, digits = 2),
        ", ", round(UL, digits = 2), "]\n", sep = "")
    table3 <- matrix(
      c("Lower", "Trivial", "Higher", negative,
        trivial, positive),
      nrow = 2,
      byrow = T
    )
    rownames(table3) <- c(" ", "MBI (%)")
    lower <-
      ifelse(negative < 0.5,
             "Most Unlikely",
             ifelse(
               negative <
                 5,
               "Very Unlikely",
               ifelse(
                 negative < 25,
                 "Unlikely",
                 ifelse(
                   negative < 75,
                   "Possibly",
                   ifelse(
                     negative <
                       95,
                     "Likely",
                     ifelse(
                       negative < 99,
                       "Most Likely",
                       ifelse(negative >= 99, "Almost Certainly")
                     )
                   )
                 )
               )
             ))
    trivial2 <-
      ifelse(trivial < 0.5,
             "Most Unlikely",
             ifelse(
               trivial <
                 5,
               "Very Unlikely",
               ifelse(
                 trivial < 25,
                 "Unlikely",
                 ifelse(
                   trivial < 75,
                   "Possibly",
                   ifelse(
                     trivial <
                       95,
                     "Likely",
                     ifelse(
                       negative < 99,
                       "Most Likely",
                       ifelse(negative >= 99, "Almost Certainly")
                     )
                   )
                 )
               )
             ))
    higher <-
      ifelse(positive < 0.5,
             "Most Unlikely",
             ifelse(
               positive <
                 5,
               "Very Unlikely",
               ifelse(
                 positive < 25,
                 "Unlikely",
                 ifelse(
                   positive < 75,
                   "Possibly",
                   ifelse(
                     positive <
                       95,
                     "Likely",
                     ifelse(
                       positive < 99,
                       "Most Likely",
                       ifelse(positive >= 99, "Almost Certainly")
                     )
                   )
                 )
               )
             ))
    colnames(table3) <- c(lower, trivial2, higher)
    title4 <- ("   Magnitude-Based Inference\n")
    cat("\n")
    cat(title4, "\n")
    print(table3)
    infer <- which.max(table3[2,])
    dir <- ifelse(OR < 1, "Decrease.", "Increase.")
    infer2 <- ifelse(infer == 1, lower, ifelse(infer == 2,
                                               trivial2, ifelse(infer == 3, higher)))
    if (OR < 1) {
      mag <- ifelse(OR > (1 / 1.5) || infer == 2,
                    "Trivial",
                    ifelse(OR > (1 / 3.5), "Small", ifelse(
                      OR > (1 / 9),
                      "Moderate", ifelse(OR > (1 /
                                                 32), "Large", ifelse(OR <=
                                                                        (1 /
                                                                           32), "Very Large"))
                    )))
    }
    else {
      mag <- ifelse(abs(OR) < 1.5 || infer == 2,
                    "Trivial",
                    ifelse(
                      abs(OR) < 3.5,
                      "Small",
                      ifelse(abs(OR) <
                               9, "Moderate", ifelse(
                                 abs(OR) < 32, "Large",
                                 ifelse(abs(OR) >= 32, "Very Large")
                               ))
                    ))
    }

    LogInference <-  ifelse(abs(positive) >= 5 && abs(negative) > 5,
                            paste("Inference: Unclear Difference."),
                            paste("Inference:", infer2, mag, dir,sep = " "))
    cat(LogInference)

    cat("\n-----------------------------------------------")
    cat("\n   Nonparametric Inference", title3, sep = "\n")
    cat("\n")
    print(table2)
    cat("\n")
    cat("   ", "r", " = ", round(r, digits = 2), "\n", sep = "")
    cat("   ", level, " CI ", "[", round(r.LL, digits = 2),
        ", ", round(r.UL, digits = 2), "]\n\n", sep = "")
    positive.r <-
      round(100 * (1 - stats::pnorm(
        0.1, mean = (swc *
                       log((1 + r) /
                             (1 - r))), sd = (1 / sqrt(n - 3))
      )), digits = 1)
    negative.r <- round(100 * (stats::pnorm(
      -0.1, mean = (swc *
                      log((1 + r) /
                            (1 - r))), sd = (1 / sqrt(n - 3))
    )), digits = 1)
    trivial.r <- round(100 - positive.r - negative.r, digits = 1)
    lower.r <- ifelse(negative.r < 0.5,
                      "Most Unlikely",
                      ifelse(
                        negative.r < 5,
                        "Very Unlikely",
                        ifelse(
                          negative.r <
                            25,
                          "Unlikely",
                          ifelse(
                            negative.r < 75,
                            "Possibly",
                            ifelse(
                              negative.r < 95,
                              "Likely",
                              ifelse(
                                negative.r <
                                  99,
                                "Most Likely",
                                ifelse(negative.r >= 99,
                                       "Almost Certainly")
                              )
                            )
                          )
                        )
                      ))
    trivial2.r <- ifelse(trivial.r < 0.5,
                         "Most Unlikely",
                         ifelse(
                           trivial.r < 5,
                           "Very Unlikely",
                           ifelse(
                             trivial.r <
                               25,
                             "Unlikely",
                             ifelse(
                               trivial.r < 75,
                               "Possibly",
                               ifelse(
                                 trivial.r < 95,
                                 "Likely",
                                 ifelse(
                                   trivial.r <
                                     99,
                                   "Most Likely",
                                   ifelse(trivial.r >= 99,
                                          "Almost Certainly")
                                 )
                               )
                             )
                           )
                         ))
    higher.r <- ifelse(positive.r < 0.5,
                       "Most Unlikely",
                       ifelse(
                         positive.r < 5,
                         "Very Unlikely",
                         ifelse(
                           positive.r <
                             25,
                           "Unlikely",
                           ifelse(
                             positive.r < 75,
                             "Possibly",
                             ifelse(
                               positive.r < 95,
                               "Likely",
                               ifelse(
                                 positive.r <
                                   99,
                                 "Most Likely",
                                 ifelse(positive.r >= 99,
                                        "Almost Certainly")
                               )
                             )
                           )
                         )
                       ))
    table.r <- matrix(
      c(
        "Negative",
        "Trivial",
        "Positive",
        negative.r,
        trivial.r,
        positive.r
      ),
      nrow = 2,
      byrow = T
    )
    rownames(table.r) <- c(" ", "MBI (%)")
    colnames(table.r) <- c(lower.r, trivial2.r, higher.r)
    title.r <- ("   Magnitude-Based Inference")
    cat(title.r, "\n\n")
    print(table.r)
    cat("\n")
    infer.r <- which.max(table.r[2,])
    infer3.r <- ifelse(infer.r == 1,
                       lower.r,
                       ifelse(infer.r ==
                                2, trivial2.r, ifelse(infer.r == 3, higher.r)))
    mag.r <- ifelse(abs(r) < 0.1 || infer.r == 2,
                    "Trivial",
                    ifelse(abs(r) < 0.3, "Small", ifelse(
                      abs(r) < 0.5,
                      "Moderate", ifelse(
                        abs(r) < 0.7,
                        "Large",
                        ifelse(abs(r) <
                                 0.9, "Very Large", ifelse(abs(r) >= 0.9, "Very Large"))
                      )
                    )))

    NonParametricInference <-  ifelse(abs(positive.r) >= 5 && abs(negative.r) > 5,
                                      paste("Inference: Unclear Difference."),
                                      paste("Inference:", infer3.r, mag.r, "Effect Size",sep = " "))
    cat(NonParametricInference)

    invisible(list(or.stat = OR, or.LL = LL, or.UL = UL, r.stat = r, r.LL = r.LL, r.UL = r.UL,
                   t.value = test$statistic[[1]], df = test$parameter[[1]], p.value = test$p.value,
                   mean1 = round(mean(x, na.rm = T), digits = 3), sd1 = round(sd(x, na.rm = T), digits = 3),
                   mean2 = round(mean(y, na.rm = T), digits = 3), sd2 = round(sd(y, na.rm = T), digits = 3),
                   LogmbiNegative = negative, LogmbiTrivial = trivial, LogmbiPositive = positive,
                   mbiNegative.r = negative.r, mbiTrivial.r = trivial.r, mbiPositive.r = positive.r,
                   LogInference=LogInference, NonParametricInference=NonParametricInference))
  }
  else {
    test <- stats::t.test(x,y, var.equal = variance, paired = paired, conf.level = conf.int,
                          na.action = na.omit, mu = mu)
    d <- ifelse(paired == T, ((mean(x, na.rm = T) - mean(y,
                                                         na.rm = T)) / pair.stdr), test$statistic * ind.stdr)
    LL <- d - (stats::qt(((
      100 - (100 * conf.int)
    ) / 100) / 2,
    test$parameter)) * abs(d) / stats::qt(test$p.value /
                                            2,
                                          test$parameter)
    UL <- d + (stats::qt(((
      100 - (100 * conf.int)
    ) / 100) / 2,
    test$parameter)) * abs(d) / stats::qt(test$p.value /
                                            2,
                                          test$parameter)
    if (test$parameter < 30) {
      d <- d * (1 - (3 / (4 * (
        test$parameter - 1
      ))))
      LL <- LL * (1 - (3 / (4 * (
        test$parameter - 1
      ))))
      UL <- UL * (1 - (3 / (4 * (
        test$parameter - 1
      ))))
    }
    parameter <- ifelse(test$parameter < 30, "Cohen d Adj",
                        "Cohen d")
    negative <-
      round(100 * (ifelse((d--swc) > 0,
                          stats::pt((d--swc) / abs(d) * abs(test$statistic),
                                    test$parameter,
                                    lower.tail = F
                          ),
                          (1 - stats::pt((-swc - d) / abs(d) *
                                           abs(test$statistic),
                                         test$parameter,
                                         lower.tail = F
                          ))
      )),
      digits = 1)
    positive <-
      round(100 * (ifelse((d - swc) > 0,
                          (1 - stats::pt((d -
                                            swc) /
                                           abs(d) * abs(test$statistic),
                                         test$parameter,
                                         lower.tail = F
                          )),
                          stats::pt((swc - d) / abs(d) * abs(test$statistic),
                                    test$parameter,
                                    lower.tail = F
                          )
      )), digits = 1)
    trivial <- round((100 - positive - negative), digits = 1)
    table1 <- matrix(c(
      round(test$statistic, digits = 3),
      round(test$parameter, digits = 0),
      round(test$p.value,
            digits = 3),
      round(test$conf.int[1], digits = 3),
      round(test$conf.int[2], digits = 3)
    ),
    ncol = 1,
    byrow = T)
    level <- paste(as.character(100 * conf.int), "%", sep = "")
    rownames(table1) <- c("t",
                          "df",
                          "P value",
                          paste(level,
                                "CI LL", sep = " "),
                          paste(level, "CI UL", sep = " "))
    colnames(table1) <- c("Raw Scale")
    table2 <- matrix(
      c("Lower", "Trivial", "Higher", negative,
        trivial, positive),
      nrow = 2,
      byrow = T
    )
    rownames(table2) <- c(" ", "MBI (%)")
    lower <-
      ifelse(negative < 0.5,
             "Most Unlikely",
             ifelse(
               negative <
                 5,
               "Very Unlikely",
               ifelse(
                 negative < 25,
                 "Unlikely",
                 ifelse(
                   negative < 75,
                   "Possibly",
                   ifelse(
                     negative <
                       95,
                     "Likely",
                     ifelse(
                       negative < 99,
                       "Most Likely",
                       ifelse(negative >= 99, "Almost Certainly")
                     )
                   )
                 )
               )
             ))
    trivial2 <-
      ifelse(trivial < 0.5,
             "Most Unlikely",
             ifelse(
               trivial <
                 5,
               "Very Unlikely",
               ifelse(
                 trivial < 25,
                 "Unlikely",
                 ifelse(
                   trivial < 75,
                   "Possibly",
                   ifelse(
                     trivial <
                       95,
                     "Likely",
                     ifelse(
                       negative < 99,
                       "Most Likely",
                       ifelse(negative >= 99, "Almost Certainly")
                     )
                   )
                 )
               )
             ))
    higher <-
      ifelse(positive < 0.5,
             "Most Unlikely",
             ifelse(
               positive <
                 5,
               "Very Unlikely",
               ifelse(
                 positive < 25,
                 "Unlikely",
                 ifelse(
                   positive < 75,
                   "Possibly",
                   ifelse(
                     positive <
                       95,
                     "Likely",
                     ifelse(
                       positive < 99,
                       "Most Likely",
                       ifelse(positive >= 99, "Almost Certainly")
                     )
                   )
                 )
               )
             ))
    colnames(table2) <- c(lower, trivial2, higher)
    title3 <- ("   Magnitude-Based Inference")
    title <- paste(normlabel, variance2, sep = ", ")
    title2 <- paste("   Method:", test$method, sep = " ")
    cat("   Mean of ", deparse(substitute(x)), " = ", round(mean(x,
                                                                 na.rm = T), digits = 2), "; ", "Mean of ", deparse(substitute(y)),
        " = ", round(mean(y, na.rm = T), digits = 2), "\n",
        sep = "")
    cat(title, title2, sep = "\n")
    cat("\n")
    print(table1)
    cat("\n")
    cat("   ", parameter, " = ", round(d, digits = 2), "\n",
        sep = "")
    cat("   ", level, " CI ", "[", round(LL, digits = 2),
        ", ", round(UL, digits = 2), "]\n\n", sep = "")
    cat(title3, "\n\n")
    print(table2)
    cat("\n")
    infer <- which.max(table2[2,])
    infer2 <- ifelse(infer == 1, lower, ifelse(infer == 2,
                                               trivial2, ifelse(infer == 3, higher)))
    mag <- ifelse(abs(d) < 0.2 || infer == 2,
                  "Trivial",
                  ifelse(abs(d) < 0.6, "Small", ifelse(
                    abs(d) < 1.2,
                    "Moderate", ifelse(abs(d) < 2, "Large", ifelse(abs(d) >=
                                                                     2, "Very Large"))
                  )))
    dir <- ifelse(d < 0, "Decrease.", "Increase.")


    Inference <-  ifelse(abs(positive) >= 5 && abs(negative) > 5,
                         paste("Inference: Unclear Difference."),
                         paste("Inference:", infer2, mag, dir, sep = " "))
    cat(Inference)

    if (plot == TRUE) {
      plot(NA, ylim = c(0, 1), xlim = c(min(LL, -swc) -
                                          max(UL - LL, swc - -swc)/10,
                                        max(UL, swc) + max(UL - LL, swc -
                                                               -swc)/10), bty = "l", yaxt = "n", ylab = "",
           xlab = "Effect Size")
      graphics::points(x = d, y = 0.5, pch = 15, cex = 2)
      graphics::abline(v = swc, lty = 2)
      graphics::abline(v = -swc, lty = 2)
      graphics::abline(v = 0, lty = 2, col = "grey")
      graphics::segments(LL, 0.5, UL, 0.5, lwd = 3)
      graphics::title(main = paste(
        "Cohen's d = ", round(d, digits = 3), " \n  ",
        100 * (conf.int), "% CI [", round(LL, digits = 3),
        ";", round(UL, digits = 3), "] ", " \n  ", "Inference: ", infer2, " ", mag," ", dir,
        sep = ""), cex.main = 1)
    }

    invisible(list(
      d.stat = d,
      d.LL = LL,
      d.UL = UL,
      t.value = test$statistic[[1]],
      df = test$parameter[[1]],
      p.value = test$p.value,
      conf.int=conf.int,
      mean1 = round(mean(x, na.rm = T), digits = 3),
      sd1 = round(sd(x, na.rm = T), digits = 3),
      mean2 = round(mean(y, na.rm = T), digits = 3),
      sd2 = round(sd(y, na.rm = T), digits = 3),
      mbiNegative = negative,
      mbiTrivial = trivial,
      mbiPositive = positive,
      Inference = Inference)
    )
  }
}
