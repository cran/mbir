
# This file is a generated template, your changes will not be overwritten

dataMBItwoClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "dataMBItwoClass",
    inherit = dataMBItwoBase,
    private = list(
      .init = function() {
        
        
        table <- self$results$ttest
        EffSize <- self$results$effect
        MBI <- self$results$mbi
        desc <- self$results$desc
        
        EffSize$getColumn('diffLL')$setSuperTitle(jmvcore::format('{}% Confidence Interval', self$options$confint))
        EffSize$getColumn('diffUL')$setSuperTitle(jmvcore::format('{}% Confidence Interval', self$options$confint))        
        EffSize$getColumn('dLL')$setSuperTitle(jmvcore::format('{}% Confidence Interval', self$options$confint))
        EffSize$getColumn('dUL')$setSuperTitle(jmvcore::format('{}% Confidence Interval', self$options$confint))
        
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
        
        if ( ! self$options$varEq)
          table$setNote('varEq', "Welch's, unequal variances, t-test")
        
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
          
          #formula <- jmvcore::constructFormula(dep, self$options$group)
          #formula <- as.formula(formula)
          
          #result <- t.test(formula, self$data, var.equal=self$options$varEq)
          
          dep <- self$data[[depName]]
          dep <- jmvcore::toNumeric(dep)
          dataTTest <- data.frame(dep=dep, group=group)
          dataTTest <- na.omit(dataTTest)
          
          n <- tapply(dataTTest$dep, dataTTest$group, length)
          n1 <- n[1]
          n2 <- n[2]
          
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
          
          result <- t.test(dep ~ group, dataTTest, var.equal=self$options$varEq)
          
          mdiff <- unname(result$estimate[1]-result$estimate[2])
          diffLL <- unname(result$conf.int[1])
          diffUL <- unname(result$conf.int[2])
          
          ind.stdr <- sqrt(1 / n[1] + 1 / n[2])
          
          d <- result$statistic * ind.stdr
          
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
          
          table$setRow(rowKey=depName, list(  
            t=result$statistic,
            df=result$parameter,
            p=result$p.value))
          
          #effect size
          EffSize$setRow(rowKey=depName, list(
            mdiff=mdiff,
            diffLL=diffLL,
            diffUL=diffUL,
            type=type,
            d=d, 
            dLL=LL, 
            dUL=UL))
          
          MBI$setRow(rowKey=depName, list(
            inference=infer3,
            negative=negative, 
            trivial=trivial, 
            positive=positive))
          
          desc$setRow(rowKey=depName, list(
            `n[1]`=n[1], `m[1]`=m[1], `med[1]`=med[1], `sd[1]`=sd[1], `se[1]`=se[1],
            `n[2]`=n[2], `m[2]`=m[2], `med[2]`=med[2], `sd[2]`=sd[2], `se[2]`=se[2]))
          #desc$setRow(rowKey=dep, values=list(
          # `n[1]`=n[1], m[1]=m[1], med[1]=med[1], sd[1]=sd[1], se[1]=se[1],
          #n[2]=n[2], m[2]=m[2], med[2]=med[2], sd[2]=sd[2], se[2]=se[2]))
          
        }
        
        
        
        
        
        # `self$data` contains the data
        # `self$options` contains the options
        # `self$results` contains the results object (to populate)
        
      }
      )
)
