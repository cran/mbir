
# This file is a generated template, your changes will not be overwritten

dataMBInpindClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "dataMBInpindClass",
    inherit = dataMBInpindBase,
    private = list(
      .init = function() {
        
        
        table <- self$results$ttest
        EffSize <- self$results$effect
        MBI <- self$results$mbi
        desc <- self$results$desc
        
        
        EffSize$getColumn('rLL')$setSuperTitle(jmvcore::format('{}% Confidence Interval', self$options$confint))
        EffSize$getColumn('rUL')$setSuperTitle(jmvcore::format('{}% Confidence Interval', self$options$confint))
        
        groupName <- self$options$group
        groups <- NULL
        if ( ! is.null(groupName))
          groups <- base::levels(self$data[[groupName]])
        if (length(groups) != 2)
          groups <- c('Group 1', 'Group 2')
        
        desc <- self$results$desc
        for (key in desc$rowKeys) {
          desc$setRow(rowKey=key, values=list(
            `name[1]`=groups[1],
            `name[2]`=groups[2]))
        }
        
 
        
      },
        .run = function() {
          
          if (is.null(self$options$group) || length(self$options$deps) == 0)
            return()
          
          table <- self$results$ttest
          EffSize <- self$results$effect
          MBI <- self$results$mbi
          desc <- self$results$desc
          
          conf.int <- (self$options$confint)/100
          SWC <- self$options$SWC
          
          groupName <- self$options$group
          group <- self$data[[groupName]]
          group <- as.factor(group)
          group <- droplevels(group)
          
          groupLevels <- base::levels(group)
          if (length(groupLevels) != 2)
            jmvcore::reject("Grouping variable must have exactly 2 levels", code="grouping_var_must_have_2_levels")
          
          
          for (depName in self$options$deps) {
            
            
            dep <- self$data[[depName]]
            dep <- jmvcore::toNumeric(dep)
            dataTTest <- data.frame(dep=dep, group=group)
            dataTTest <- na.omit(dataTTest)
            
            n <- tapply(dataTTest$dep, dataTTest$group, length)
            n1 <- n[1]
            n2 <- n[2]
            n3 <- n1+n2
            
            v <- tapply(dataTTest$dep, dataTTest$group, function(x) jmvcore::tryNaN(var(x)))
            v1 <- v[1]
            v2 <- v[2]
            
            m <- tapply(dataTTest$dep, dataTTest$group, function(x) jmvcore::tryNaN(mean(x)))
            m1 <- m[1]
            m2 <- m[2]
            
            med <- tapply(dataTTest$dep, dataTTest$group, function(x) jmvcore::tryNaN(median(x)))
            med1 <- med[1]
            med2 <- med[2]
            
            se <- sqrt(v/n)
            se1 <- se[1]
            se2 <- se[2]
            
            sd <- sqrt(v)
            sd1 <- sd[1]
            sd2 <- sd[2]
            
            
            rank <- stats::wilcox.test(dep ~ group, dataTTest,
              paired = FALSE,
              conf.int = T,
              conf.level = conf.int,
              na.action = na.omit,
              correct = F,
              exact = F
            )
            
            suppressWarnings(warning(rank))
            diff <- rank$estimate
            
            r <- ifelse(diff < 0,
                        stats::qnorm(rank$p.value / 2) / sqrt(n3),
                        abs(stats::qnorm(rank$p.value / 2) / sqrt(n3)))
            
            r.LL <-( (exp(2 * ((SWC * log((1 + r) / (1 - r)
              )) + (
                stats::qnorm(((
                  100 -
                    (100 * conf.int)
                ) / 100 / 2)) / sqrt(n3 - 3)
              ))) - 1) / (exp(2 *
                                ((SWC * log((1 + r) /
                                                (1 - r)
                                )) + (
                                  stats::qnorm(((
                                    100 -
                                      (100 * conf.int)
                                  ) / 100 / 2)) / sqrt(n3 - 3)
                                ))) + 1) )
            
            rLL <- as.numeric(r.LL)
              
            r.UL <- ((exp(2 * ((SWC * log((1 + r) / (1 - r)
              )) - (
                stats::qnorm(((
                  100 -
                    (100 * conf.int)
                ) / 100 / 2)) / sqrt(n3 - 3)
              ))) - 1) / (exp(2 *
                                ((SWC * log((1 + r) /
                                                (1 - r)
                                )) - (
                                  stats::qnorm(((
                                    100 -
                                      (100 * conf.int)
                                  ) / 100 / 2)) / sqrt(n3 - 3)
                                ))) + 1))
            
            rUL <- as.numeric(r.UL)
            
            level <- paste(as.character(100 * conf.int), "%", sep = "")
            
            positive.r <-
              round(100 * (1 - stats::pnorm(
                0.1, mean = (SWC *
                               log((1 + r) /
                                     (1 - r))), sd = (1 / sqrt(n3 - 3))
              )), digits = 1)
            negative.r <- round(100 * (stats::pnorm(
              -0.1, mean = (SWC *
                              log((1 + r) /
                                    (1 - r))), sd = (1 / sqrt(n3 - 3))
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
            
            table$setRow(rowKey=depName, list(  
              statistic=rank$statistic,
              p=rank$p.value))
            
            #effect size
            EffSize$setRow(rowKey=depName, list(
              r=r, 
              rLL=rLL, 
              rUL=rUL))
            
            MBI$setRow(rowKey=depName, list(
              inference=Inference,
              negative=negative.r, 
              trivial=trivial.r, 
              positive=positive.r))
            
            desc$setRow(rowKey=depName, list(
              `n[1]`=n[1], `m[1]`=m[1], `med[1]`=med[1], `sd[1]`=sd[1], `se[1]`=se[1],
              `n[2]`=n[2], `m[2]`=m[2], `med[2]`=med[2], `sd[2]`=sd[2], `se[2]`=se[2]))
            
          }
            # `self$data` contains the data
            # `self$options` contains the options
            # `self$results` contains the results object (to populate)

        })
)
