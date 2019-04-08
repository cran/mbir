
# This file is a generated template, your changes will not be overwritten

dataMBIpairedClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "dataMBIpairedClass",
    inherit = dataMBIpairedBase,
    private = list(
        .init = function(){
          
          table <- self$results$ttest
          EffSize <- self$results$effect
          MBI <- self$results$mbi
          desc <- self$results$desc
          

          
          for (pair in self$options$pairs) {
            table$setRow(rowKey=pair,  list(i1=pair[[1]], i2=pair[[2]]))
            EffSize$setRow(rowKey=pair, list(i1=pair[[1]], i2=pair[[2]]))
            MBI$setRow(rowKey=pair, list(i1=pair[[1]], i2=pair[[2]]))
            desc$setRow(rowKey=pair, list(`name[1]`=pair[[1]], `name[2]`=pair[[2]]))
          }
          
          EffSize$getColumn('diffLL')$setSuperTitle(jmvcore::format('{}% Confidence Interval', self$options$confint))
          EffSize$getColumn('diffUL')$setSuperTitle(jmvcore::format('{}% Confidence Interval', self$options$confint)) 
          EffSize$getColumn('dLL')$setSuperTitle(jmvcore::format('{}% Confidence Interval', self$options$confint))
          EffSize$getColumn('dUL')$setSuperTitle(jmvcore::format('{}% Confidence Interval', self$options$confint))          
          
          
        },
        .run = function() {

          
          table <- self$results$ttest
          EffSize <- self$results$effect
          MBI <- self$results$mbi
          desc <- self$results$desc
          
          conf.int <- (self$options$confint)/100
          SWC <- self$options$SWC
          
          for (pair in self$options$pairs) {
            
            if (is.null(pair[[1]]))
              next()
            if (is.null(pair[[2]]))
              next()
            
            #Gather descriptives
            i1 <- jmvcore::toNumeric(self$data[[pair[[1]] ]])
            i2 <- jmvcore::toNumeric(self$data[[pair[[2]] ]])
            data <- data.frame(i1=i1, i2=i2)
            data <- na.omit(data)
            n <- nrow(data)
            i1 <- data$i1
            i2 <- data$i2
            m1 <- base::mean(i1)
            m2 <- base::mean(i2)
            med1 <- stats::median(i1)
            med2 <- stats::median(i2)
            sd1  <- stats::sd(i1)
            sd2  <- stats::sd(i2)
            se1  <- sd1/sqrt(n)
            se2  <- sd2/sqrt(n)
            
            result <- t.test(i1, i2, paired=TRUE)
            corr <- cor.test(x=i1,y=i2)
            corr <- unname(corr$estimate)
            
            #Mean Difference
            mdiff <- m1-m2
            diffSE <- sqrt(((sd1^2/n)+(sd2^2/n))-(2*corr*(sd1/sqrt(n))*(sd2/sqrt(n))))
            diffLL <- mdiff+diffSE*qt((1-conf.int)/2,n)
            diffUL <- mdiff-diffSE*qt((1-conf.int)/2,n)
            
            #Effect Size
            d <- (result$statistic/sqrt(n))
            
            LL <- d - (stats::qt(((
              100 - (100 * conf.int)
            ) / 100) / 2,
            result$parameter)) * abs(d) / stats::qt(result$p.value /
                                                      2,
                                                    result$parameter)
            
            UL <- d + (stats::qt(((
              100 - (100 * conf.int)
            ) / 100) / 2,
            result$parameter)) * abs(d) / stats::qt(result$p.value /
                                                      2,
                                                    result$parameter)
            type <- paste("Cohen's d")
            
            #Hedges Correlation
            if (result$parameter < 30) {
              type <- paste("Hedge's g")
              d <- d * (1 - (3 / (4 * (
                result$parameter - 1
              ))))
              LL <- LL * (1 - (3 / (4 * (
                result$parameter - 1
              ))))
              UL <- UL * (1 - (3 / (4 * (
                result$parameter - 1
              ))))
            }
            

            
            #Correct for correlation
            if (self$options$correl==TRUE) {
              
              d <- d * sqrt(2*(1 - corr))
              
              LL <- d - (stats::qt(((
                100 - (100 * conf.int)
              ) / 100) / 2,
              result$parameter)) * abs(d) / stats::qt(result$p.value /
                                                        2,
                                                      result$parameter)
              UL <- d + (stats::qt(((
                100 - (100 * conf.int)
              ) / 100) / 2,
              result$parameter)) * abs(d) / stats::qt(result$p.value /
                                                        2,
                                                      result$parameter)
            }
            
            
            ###MBI
            negative <- round(100 * (ifelse((d--SWC) > 0,
                                            stats::pt((d--SWC) / abs(d) * abs(result$statistic),
                                                      result$parameter,
                                                      lower.tail = F
                                            ),
                                            (1 - stats::pt((-SWC - d) / abs(d) *
                                                             abs(result$statistic),
                                                           result$parameter,
                                                           lower.tail = F
                                            ))
            )),
            digits = 1)
            positive <- round(100 * (ifelse((d - SWC) > 0,
                                            (1 - stats::pt((d -
                                                              SWC) /
                                                             abs(d) * abs(result$statistic),
                                                           result$parameter,
                                                           lower.tail = F
                                            )),
                                            stats::pt((SWC - d) / abs(d) * abs(result$statistic),
                                                      result$parameter,
                                                      lower.tail = F
                                            )
            )), digits = 1)
            trivial <- round((100 - positive - negative), digits = 1)
            
            lower <- ifelse(negative < 0.5,
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
            trivial2 <- ifelse(trivial < 0.5,
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
                                         trivial < 99,
                                         "Most Likely",
                                         ifelse(trivial >= 99, "Almost Certainly")
                                       )
                                     )
                                   )
                                 )
                               ))
            higher <- ifelse(positive < 0.5,
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
            
            table2 <- matrix(
              c("Lower", "Trivial", "Higher", negative,
                trivial, positive),
              nrow = 2,
              byrow = T
            )
            
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
            
            infer3 <-  ifelse(abs(positive) >= 5 && abs(negative) > 5,
                              paste("Unclear Difference."),
                              paste(infer2, mag, dir, sep = " "))
            
            #Add MBI results
            
            #Results
            #T-test and effect size
            table$setRow(rowKey=pair, list(  
              t=result$statistic,
              df=result$parameter,
              p=result$p.value))
            
            #effect size
            EffSize$setRow(rowKey=pair, list(
              mdiff=mdiff,
              diffLL=diffLL,
              diffUL=diffUL,
              type=type,
              d=d, 
              dLL=LL, 
              dUL=UL))
            
            #MBI results
            MBI$setRow(rowKey=pair, list(
              inference=infer3,
              negative=negative, 
              trivial=trivial, 
              positive=positive))
            
            #Descriptives results
            desc$setRow(rowKey=pair, list(
              `n[1]`=n, `m[1]`=m1, `med[1]`=med1, `sd[1]`=sd1, `se[1]`=se1,
              `n[2]`=n, `m[2]`=m2, `med[2]`=med2, `sd[2]`=sd2, `se[2]`=se2))
            
          }
            # `self$data` contains the data
            # `self$options` contains the options
            # `self$results` contains the results object (to populate)

        })
)
