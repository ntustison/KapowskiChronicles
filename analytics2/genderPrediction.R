# Load required libraries

library( MASS )
library( popbio )
library( gdata )
library( ROCR )
library( caret )
library( ggplot2 )

nPermutations <- 1000
trainingPortions <- c( 0.5 )

resultsData <- data.frame( Pipeline = character( 0 ), Correlation = numeric( 0 ), MeanErrors = numeric( 0 ) )
aucData <- data.frame( Pipeline = character( 0 ), FPR = numeric( 0 ), TPR = numeric( 0 ) )


for( trainingPortion in trainingPortions )
  {
  cat( "trainingPortion = ", trainingPortion, "\n", sep = '' );
  for( n in seq( 1, nPermutations, by = 1 ) )
    {
    cat( "  Permutation ", n, "\n", sep = '' )

    thicknessTypes = c( 'ANTs', 'FreeSurfer' )
    for( whichPipeline in thicknessTypes )
      {
      resultsIXI <- read.csv( paste0( 'labelresults', whichPipeline, 'I.csv' ) )
      resultsKirby <- read.csv( paste0( 'labelresults', whichPipeline, 'K.csv' ) )
      resultsNKI <- read.csv( paste0( 'labelresults', whichPipeline, 'N.csv' ) )
      resultsOasis <- read.csv( paste0( 'labelresults', whichPipeline, 'O.csv' ) )

      resultsCombined <- rbind( resultsIXI, resultsKirby, resultsNKI, resultsOasis );
      resultsCombined$SITE <- as.factor( resultsCombined$SITE )
      resultsCombined$SEX <- resultsCombined$SEX - 1

      corticalLabels <- tail( colnames( resultsCombined ), n = 62 )

      drops <- c( "ID", "SITE" )
      resultsCombined <- resultsCombined[, !( names( resultsCombined ) %in% drops )]

  #     trainingIndices <- createDataPartition( resultsCombined$SEX, p = trainingPortion, list = FALSE, times = 1 )
      trainingIndices <- sample( 1:nrow( resultsCombined ), size = 550, replace = FALSE )
      trainingData <- resultsCombined[trainingIndices,]
      testingData <- resultsCombined[-trainingIndices,]

      #
      # Do gender prediction
      #   Include all AGE/REGION interaction terms
      #

      regionalTerms <- paste( corticalLabels, collapse = " + " )
      myFormula <- as.formula( paste( "SEX ~ AGE + ", regionalTerms, " + VOLUME ", sep = '' ) )

      genderModel <- glm( myFormula, data = trainingData, family = "binomial" )
      genderPrediction <- predict( genderModel, newdata = testingData, type = 'response' )

      genderPredictionError <- abs( mean( genderPrediction - testingData$SEX, na.rm = TRUE ) )

      genderPredictionROC <- prediction( genderPrediction, testingData$SEX )
      genderPerformanceROC <- performance( genderPredictionROC, "tpr", "fpr" )
      auc <- performance( genderPredictionROC, "auc" )@y.values[[1]]

      oneData <- data.frame( Pipeline = whichPipeline, AUC = auc, MeanError = genderPredictionError )
      resultsData <- rbind( resultsData, oneData )

      fpr <- genderPerformanceROC@x.values[[1]]
      tpr <- genderPerformanceROC@y.values[[1]]
      whichPipelineIndices <- which( aucData$Pipeline == whichPipeline )

      if( n == 1 )
        {
        oneAucData <- data.frame( Pipeline = rep( whichPipeline, length( fpr ) ),
                                  FPR = fpr, TPR = tpr )
        aucData <- rbind( aucData, oneAucData )
        }
      else
        {
        myApprox <- approx( fpr, tpr, n = length( whichPipelineIndices ) )
        fpr <- myApprox$x
        tpr <- myApprox$y

        aucData$FPR[whichPipelineIndices] <- aucData$FPR[whichPipelineIndices] + ( fpr - aucData$FPR[whichPipelineIndices] ) / ( n + 1 )
        aucData$TPR[whichPipelineIndices] <- aucData$TPR[whichPipelineIndices] + ( tpr - aucData$TPR[whichPipelineIndices] ) / ( n + 1 )
        }
      }
    }

  aucPlot <- ggplot( resultsData, aes( x = AUC, fill = Pipeline ) ) +
                     scale_y_continuous( "Density" ) +
                     geom_density( alpha = 0.5 )
  ggsave( filename = paste( "~/Desktop/genderAuc", trainingPortion, ".pdf", sep = "" ), plot = aucPlot, width = 6, height = 6, units = 'in' )

  meanErrorPlot <- ggplot( resultsData, aes( x = MeanError, fill = Pipeline ) ) +
                     scale_y_continuous( "Density" ) +
                     geom_density( alpha = 0.5 )
  ggsave( filename = paste( "~/Desktop/genderMeanError", trainingPortion, ".pdf", sep = "" ), plot = meanErrorPlot, width = 6, height = 6, units = 'in' )

  for( whichPipeline in thicknessTypes )
    {
    genderPlot <- ggplot( aucData, aes( x = FPR, y = TPR, colour = Pipeline ) ) +
          geom_line( size = 1.0 ) +
          geom_abline ( intercept = 0, slope = 1, colour = "black", linetype = "dashed", size = 1.0 ) +
          scale_x_continuous( "False positive rate (1-specificity)" ) +
          scale_y_continuous( "True positive rate (sensitivity)" )
    ggsave( filename = paste( "~/Desktop/genderPlot", trainingPortion, ".pdf", sep = "" ),
            plot = genderPlot, width = 8, height = 6, units = 'in' )
    }

  aucFS <- resultsData$AUC[which( resultsData$Pipeline == 'FreeSurfer' )]
  aucANTs <- resultsData$AUC[which( resultsData$Pipeline == 'ANTs' )]

  cat( "Mean FS auc = ", mean( aucFS, na.rm = TRUE ), "\n", sep = '' );
  cat( "Mean ANTs auc = ", mean( aucANTs, na.rm = TRUE ), "\n", sep = '' );

  myTest <- t.test( x = aucFS, y = aucANTs, alternative = "two.sided", paired = TRUE )
  cat( "p.value = ", myTest$p.value, "\n", sep = '' );
  }
