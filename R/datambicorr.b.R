
# This file is a generated template, your changes will not be overwritten

dataMBIcorrClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "dataMBIcorrClass",
    inherit = dataMBIcorrBase,
    private = list(
      .init = function() {
        
        
        MBI <- self$results$mbi
        desc <- self$results$desc
        
        confint <- self$options$confint / 100
        
        
        for (pair in self$options$pairs) {
          MBI$setRow(rowKey=pair, list(i1=pair[[1]], i2=pair[[2]]))
          desc$setRow(rowKey=pair, list(`name[1]`=pair[[1]], `name[2]`=pair[[2]]))
        }
        
        

        
        
        
      },
      
        .run = function() {

          MBI <- self$results$mbi
          desc <- self$results$desc
          
          
          conf.int <- (self$options$confint)/100
          SWC <- self$options$SWC
          threshold <- (0.5*log((1+SWC)/(1-SWC)))         
          
          for (pair in self$options$pairs) {
            
            if (is.null(pair[[1]]))
              next()
            if (is.null(pair[[2]]))
              next()
            
            i1 <- jmvcore::toNumeric(self$data[[ pair[[1]] ]])
            i2 <- jmvcore::toNumeric(self$data[[ pair[[2]] ]])
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
            
            if (self$options$spear) {
              
            result <- stats::cor.test(i1, i2, method="spearman", exact=F, conf.level=conf.int)
            
            type_Spear <- paste("Spearman")
            type_Spear2 <- paste("rho =")
            p_Spear <- unname(result$p.value)
            correl_Spear <- unname(result$estimate)
            z_Spear <- (0.5 * 
                         log((1 + correl_Spear)/(1 - correl_Spear)))
            
            zSE_Spear <- (1/sqrt(n-3))
            
            Positive_Spear <- round(100 * (1 - stats::pnorm(threshold, mean = z_Spear, sd = (1/sqrt(n - 3)))), digits = 1)
            Negative_Spear <- round(100 * (stats::pnorm(-threshold, mean = z_Spear, sd = (1/sqrt(n - 3)))), digits = 1)
            Trivial_Spear <- round(100 - Positive_Spear - Negative_Spear, digits = 1)
            
            
            lower <- ifelse(Negative_Spear < 0.5, "Most Unlikely", ifelse(Negative_Spear < 
                                                                           5, "Very Unlikely", ifelse(Negative_Spear < 25, "Unlikely", 
                                                                                                      ifelse(Negative_Spear < 75, "Possibly", ifelse(Negative_Spear < 95, 
                                                                                                                                                    "Likely", ifelse(Negative_Spear < 99, "Most Likely", ifelse(Negative_Spear >= 
                                                                                                                                                                                                                 99, "Almost Certainly")))))))
            trivial2 <- ifelse(Trivial_Spear < 0.5, "Most Unlikely", ifelse(Trivial_Spear < 
                                                                             5, "Very Unlikely", ifelse(Trivial_Spear < 25, "Unlikely", 
                                                                                                        ifelse(Trivial_Spear < 75, "Possibly", ifelse(Trivial_Spear < 95, 
                                                                                                                                                     "Likely", ifelse(Trivial_Spear < 99, "Most Likely", ifelse(Trivial_Spear >= 
                                                                                                                                                                                                             99, "Almost Certainly")))))))
            higher <- ifelse(Positive_Spear < 0.5, "Most Unlikely", ifelse(Positive_Spear < 
                                                                            5, "Very Unlikely", ifelse(Positive_Spear < 25, "Unlikely", 
                                                                                                       ifelse(Positive_Spear < 75, "Possibly", ifelse(Positive_Spear < 95, 
                                                                                                                                                     "Likely", ifelse(Positive_Spear < 99, "Most Likely", ifelse(Positive_Spear >= 
                                                                                                                                                                                                                  99, "Almost Certainly")))))))
            
            
            
            
            table <- matrix(c("Negative", "Trivial", "Positive", Negative_Spear, 
                              Trivial_Spear, Positive_Spear), nrow = 2, byrow = T)
            infer <- which.max(table[2, ])
            infer2 <- ifelse(correl_Spear < 0, "Negative", "Positive")
            infer3 <- ifelse(infer == 1, lower, ifelse(infer == 2, trivial2, 
                                                       ifelse(infer == 3, higher)))
            
            mag <- ifelse(abs(correl_Spear) < 0.1 || infer == 2, "Trivial", ifelse(abs(correl_Spear) < 
                                                                                    0.3, "Small", ifelse(abs(correl_Spear) < 0.5, "Moderate", ifelse(abs(correl_Spear) < 
                                                                                                                                                      0.7, "Large", ifelse(abs(correl_Spear) < 0.9, "Very Large", ifelse(abs(correl_Spear) >= 
                                                                                                                                                                                                                          0.9, "Very Large"))))))
            
            Inference_Spear <- ifelse(abs(Positive_Spear) >= 5 && abs(Negative_Spear) > 5, 
                                     paste("Unclear Association."),
                                     paste(infer3, mag, infer2, "Correlation.", 
                                           sep = " "))                                                                                                                                                                               
            
           
            #Save result for Pearson
            MBI$setRow(rowKey=pair, list(
              `type[1]`=type_Spear, `p[1]`=p_Spear,
              `type_Spear2[1]`=type_Spear2, `correl[1]`=correl_Spear,  `z[1]`=z_Spear,
              `zSE[1]`=zSE_Spear, `Positive[1]`=Positive_Spear,
              `Trivial[1]`=Trivial_Spear, `Negative[1]`=Negative_Spear,
              `Inference[1]`=Inference_Spear))
            
            #MBI$setRow(rowKey=pair, list(
            #  `type_Spear[1]`=type_Spear, `p_Spear[1]`=p_Spear,
            #  `type_Spear2[1]`=type_Spear2, `correl_Spear[1]`=correl_Spear,  `z_Spear[1]`=z_Spear,
            #  `zSE_Spear[1]`=zSE_Spear, `Positive_Spear[1]`=Positive_Spear,
            #  `Trivial_Spear[1]`=Trivial_Spear, `Negative_Spear[1]`=Negative_Spear,
            #  `Inference_Spear[1]`=Inference_Spear))
            
            }
            
            if (self$options$pearson) {
              
              result <- stats::cor.test(i1, i2, method="pearson", exact=F, conf.level=conf.int)
              
              type_Pear <- paste("Pearson")
              type_Pear2 <- paste("r =")
              p_Pear <- unname(result$p.value)
              correl_Pear <- unname(result$estimate)
              z_Pear <- (0.5 * 
                           log((1 + correl_Pear)/(1 - correl_Pear)))
              
              zSE_Pear <- (1/sqrt(n-3))
              
              Positive_Pear <- round(100 * (1 - stats::pnorm(threshold, mean = z_Pear, sd = (1/sqrt(n - 3)))), digits = 1)
              Negative_Pear <- round(100 * (stats::pnorm(-threshold, mean = z_Pear, sd = (1/sqrt(n - 3)))), digits = 1)
              Trivial_Pear <- round(100 - Positive_Pear - Negative_Pear, digits = 1)
              
              
              lower <- ifelse(Negative_Pear < 0.5, "Most Unlikely", ifelse(Negative_Pear < 
                                                                              5, "Very Unlikely", ifelse(Negative_Pear < 25, "Unlikely", 
                                                                                                         ifelse(Negative_Pear < 75, "Possibly", ifelse(Negative_Pear < 95, 
                                                                                                                                                  "Likely", ifelse(Negative_Pear < 99, "Most Likely", ifelse(Negative_Pear >= 
                                                                                                                                                                                                          99, "Almost Certainly")))))))
              trivial2 <- ifelse(Trivial_Pear < 0.5, "Most Unlikely", ifelse(Trivial_Pear < 
                                                                                5, "Very Unlikely", ifelse(Trivial_Pear < 25, "Unlikely", 
                                                                                                           ifelse(Trivial_Pear < 75, "Possibly", ifelse(Trivial_Pear < 95, 
                                                                                                                                                   "Likely", ifelse(Trivial_Pear < 99, "Most Likely", ifelse(Trivial_Pear >= 
                                                                                                                                                                                                           99, "Almost Certainly")))))))
              higher <- ifelse(Positive_Pear < 0.5, "Most Unlikely", ifelse(Positive_Pear < 
                                                                               5, "Very Unlikely", ifelse(Positive_Pear < 25, "Unlikely", 
                                                                                                          ifelse(Positive_Pear < 75, "Possibly", ifelse(Positive_Pear < 95, 
                                                                                                                                                   "Likely", ifelse(Positive_Pear < 99, "Most Likely", ifelse(Positive_Pear >= 
                                                                                                                                                                                                           99, "Almost Certainly")))))))
              
              
              
              
              
              table <- matrix(c("Negative", "Trivial", "Positive", Negative_Pear, 
                                Trivial_Pear, Positive_Pear), nrow = 2, byrow = T)
              infer <- which.max(table[2, ])
              infer2 <- ifelse(correl_Pear < 0, "Negative", "Positive")
              infer3 <- ifelse(infer == 1, lower, ifelse(infer == 2, trivial2, 
                                                         ifelse(infer == 3, higher)))
              mag <- ifelse(abs(correl_Pear) < 0.1 || infer == 2, "Trivial", ifelse(abs(correl_Pear) < 
                                                                            0.3, "Small", ifelse(abs(correl_Pear) < 0.5, "Moderate", ifelse(abs(correl_Pear) < 
                                                                                                                                    0.7, "Large", ifelse(abs(correl_Pear) < 0.9, "Very Large", ifelse(abs(correl_Pear) >= 
                                                                                                                                                                                              0.9, "Very Large"))))))
              
              Inference_Pear <- ifelse(abs(Positive_Pear) >= 5 && abs(Negative_Pear) > 5, 
                                        paste("Unclear Association."),
                                        paste(infer3, mag, infer2, "Correlation.", 
                                              sep = " "))                                                                                                                                                                               
              
              

#Save result for Pearson
              
              MBI$setRow(rowKey=pair, list(
                `type[0]`=type_Pear, `p[0]`=p_Pear,
                `type_Pear2[0]`=type_Pear2, `correl[0]`=correl_Pear,  `z[0]`=z_Pear,
                `zSE[0]`=zSE_Pear, `Positive[0]`=Positive_Pear,
                `Trivial[0]`=Trivial_Pear, `Negative[0]`=Negative_Pear,
                `Inference[0]`=Inference_Pear))
              
              #MBI$setRow(rowKey=pair, list(
              #  type_Pear=type_Pear, p_Pear=unname(result$p.value),
              #  type_Pear2=type_Pear2, correl_Pear=correl_Pear,  z_Pear=z_Pear,
              #  zSE_Pear=zSE_Pear, Positive_Pear=Positive_Pear,
              #  Trivial_Pear=Trivial_Pear, Negative_Pear=Negative_Pear,
              #  Inference_Pear=Inference_Pear))
                            
            }

            if (self$options$kendall) {
              
              result <- stats::cor.test(i1, i2, method="kendall", exact=F, conf.level=conf.int)
              
              type_Ken <- paste("Kendall")
              type_Ken2 <- paste("tau =")
              p_Ken <- unname(result$p.value)
              correl_Ken <- unname(result$estimate)
              z_Ken <- (0.5 * 
                           log((1 + correl_Ken)/(1 - correl_Ken)))
              
              zSE_Ken <- (1/sqrt(n-3))
              
              Positive_Ken <- round(100 * (1 - stats::pnorm(threshold, mean = z_Ken, sd = (1/sqrt(n - 3)))), digits = 1)
              Negative_Ken <- round(100 * (stats::pnorm(-threshold, mean = z_Ken, sd = (1/sqrt(n - 3)))), digits = 1)
              Trivial_Ken <- round(100 - Positive_Ken - Negative_Ken, digits = 1)
              
              
              lower <- ifelse(Negative_Ken < 0.5, "Most Unlikely", ifelse(Negative_Ken < 
                                                                             5, "Very Unlikely", ifelse(Negative_Ken < 25, "Unlikely", 
                                                                                                        ifelse(Negative_Ken < 75, "Possibly", ifelse(Negative_Ken < 95, 
                                                                                                                                                      "Likely", ifelse(Negative_Ken < 99, "Most Likely", ifelse(Negative_Ken >= 
                                                                                                                                                                                                                   99, "Almost Certainly")))))))
              trivial2 <- ifelse(Trivial_Ken < 0.5, "Most Unlikely", ifelse(Trivial_Ken < 
                                                                               5, "Very Unlikely", ifelse(Trivial_Ken < 25, "Unlikely", 
                                                                                                          ifelse(Trivial_Ken < 75, "Possibly", ifelse(Trivial_Ken < 95, 
                                                                                                                                                       "Likely", ifelse(Trivial_Ken  < 99, "Most Likely", ifelse(Trivial_Ken  >= 
                                                                                                                                                                                                               99, "Almost Certainly")))))))
              higher <- ifelse(Positive_Ken < 0.5, "Most Unlikely", ifelse(Positive_Ken < 
                                                                              5, "Very Unlikely", ifelse(Positive_Ken < 25, "Unlikely", 
                                                                                                         ifelse(Positive_Ken < 75, "Possibly", ifelse(Positive_Ken < 95, 
                                                                                                                                                       "Likely", ifelse(Positive_Ken < 99, "Most Likely", ifelse(Positive_Ken >= 
                                                                                                                                                                                                                    99, "Almost Certainly")))))))
              
              
              
              
              table <- matrix(c("Negative", "Trivial", "Positive", Negative_Ken, 
                                Trivial_Ken, Positive_Ken), nrow = 2, byrow = T)
              infer <- which.max(table[2, ])
              infer2 <- ifelse(correl_Ken < 0, "Negative", "Positive")
              infer3 <- ifelse(infer == 1, lower, ifelse(infer == 2, trivial2, 
                                                         ifelse(infer == 3, higher)))
              mag <- ifelse(abs(correl_Ken) < 0.1 || infer == 2, "Trivial", ifelse(abs(correl_Ken) < 
                                                                                       0.3, "Small", ifelse(abs(correl_Ken) < 0.5, "Moderate", ifelse(abs(correl_Ken) < 
                                                                                                                                                          0.7, "Large", ifelse(abs(correl_Ken) < 0.9, "Very Large", ifelse(abs(correl_Ken) >= 
                                                                                                                                                                                                                               0.9, "Very Large"))))))
              
              Inference_Ken <- ifelse(abs(Positive_Ken) >= 5 && abs(Negative_Ken) > 5, 
                                        paste("Unclear Association."),
                                        paste(infer3, mag, infer2, "Correlation.", 
                                              sep = " "))                                                                                                                                                                               
              
              #Save result for Pearson
              MBI$setRow(rowKey=pair, list(
                `type[2]`=type_Ken, `p[2]`=p_Ken,
                `type_Ken2[2]`=type_Ken2, `correl[2]`=correl_Ken,  `z[2]`=z_Ken,
                `zSE[2]`=zSE_Ken, `Positive[2]`=Positive_Ken,
                `Trivial[2]`=Trivial_Ken, `Negative[2]`=Negative_Ken,
                `Inference[2]`=Inference_Ken))
              
             # MBI$setRow(rowKey=pair, list(
             #   type_Ken=type_Ken, p_Ken=unname(result$p.value),
             #   type_Ken2=type_Ken2, correl_Ken=correl_Ken,  z_Ken=z_Ken,
             #   zSE_Ken=zSE_Ken, Positive_Ken=Positive_Ken,
             #   Trivial_Ken=Trivial_Ken, Negative_Ken=Negative_Ken,
             #   Inference_Ken=Inference_Ken))             
                            
            }
            
            #Save result for Descriptives
            desc$setRow(rowKey=pair, list(
              `n[1]`=n, `m[1]`=m1, `med[1]`=med1, `sd[1]`=sd1, `se[1]`=se1,
              `n[2]`=n, `m[2]`=m2, `med[2]`=med2, `sd[2]`=sd2, `se[2]`=se2))
            

          }
          

            # `self$data` contains the data
            # `self$options` contains the options
            # `self$results` contains the results object (to populate)

        })
)
