
# This file is a generated template, your changes will not be overwritten

dataMBInppairClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "dataMBInppairClass",
    inherit = dataMBInppairBase,
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
        
        EffSize$getColumn('rLL')$setSuperTitle(jmvcore::format('{}% Confidence Interval', self$options$confint))
        EffSize$getColumn('rUL')$setSuperTitle(jmvcore::format('{}% Confidence Interval', self$options$confint))          
        
        
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
           
            rank <- stats::wilcox.test(x=i1, y=i2, 
                                       paired = TRUE,
                                       conf.int = T,
                                       conf.level = conf.int,
                                       na.action = na.omit,
                                       correct = F,
                                       exact = F)  
            suppressWarnings(warning(rank))
            diff <- rank$estimate
            
            r <- ifelse(diff < 0,
                        stats::qnorm(rank$p.value / 2) / sqrt(n),
                        abs(stats::qnorm(rank$p.value / 2) / sqrt(n)))
            
            r.LL <-( (exp(2 * ((SWC * log((1 + r) / (1 - r)
            )) + (
              stats::qnorm(((
                100 -
                  (100 * conf.int)
              ) / 100 / 2)) / sqrt(n - 3)
            ))) - 1) / (exp(2 *
                              ((SWC * log((1 + r) /
                                              (1 - r)
                              )) + (
                                stats::qnorm(((
                                  100 -
                                    (100 * conf.int)
                                ) / 100 / 2)) / sqrt(n - 3)
                              ))) + 1) )
            
            rLL <- as.numeric(r.LL)
            
            r.UL <- ((exp(2 * ((SWC * log((1 + r) / (1 - r)
            )) - (
              stats::qnorm(((
                100 -
                  (100 * conf.int)
              ) / 100 / 2)) / sqrt(n - 3)
            ))) - 1) / (exp(2 *
                              ((SWC * log((1 + r) /
                                              (1 - r)
                              )) - (
                                stats::qnorm(((
                                  100 -
                                    (100 * conf.int)
                                ) / 100 / 2)) / sqrt(n - 3)
                              ))) + 1))
            
            rUL <- as.numeric(r.UL)
            
            level <- paste(as.character(100 * conf.int), "%", sep = "")
            
            positive.r <-
              round(100 * (1 - stats::pnorm(
                0.1, mean = (SWC *
                               log((1 + r) /
                                     (1 - r))), sd = (1 / sqrt(n - 3))
              )), digits = 1)
            negative.r <- round(100 * (stats::pnorm(
              -0.1, mean = (SWC *
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
            table.r <- matrix(c("Negative", "Trivial","Positive",negative.r,trivial.r,positive.r),
                              nrow = 2,
                              byrow = T
            )
            
            
            
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
            
            Inference <-  ifelse(abs(positive.r) >= 5 && abs(negative.r) > 5,
                                 paste("Unclear Difference."),
                                 paste(infer3.r, mag.r, "Effect Size",sep = " "))
            
            table$setRow(rowKey=pair, list(  
              statistic=rank$statistic,
              p=rank$p.value))
            
            #effect size
            EffSize$setRow(rowKey=pair, list(
              r=r, 
              rLL=rLL, 
              rUL=rUL))
            
            MBI$setRow(rowKey=pair, list(
              inference=Inference,
              negative=negative.r, 
              trivial=trivial.r, 
              positive=positive.r))
            
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
