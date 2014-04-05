library( kernlab )
library( caret )
library( randomForest )
library( ggplot2 )

nPermutations <- 1000

trainingPortions <- c( 0.5 )

for( p in trainingPortions )
  {
  trainingPortion <- p
  cat( "trainingPortion = ", trainingPortion, "\n", sep = '' )

  resultsData <- data.frame( Pipeline = character( 0 ), Correlation = numeric( 0 ), MeanErrors = numeric( 0 ) )

  for( n in seq( 1, nPermutations, by = 1 ) )
    {
    cat( "  Permutation ", n, "\n", sep = '' )

    thicknessTypes = c( 'ANTs', 'FreeSurfer' )
    for( whichPipeline in thicknessTypes )
      {
#       resultsIXI <- read.csv( paste0( 'labelresults', whichPipeline, 'I.csv' ) )
#       resultsKirby <- read.csv( paste0( 'labelresults', whichPipeline, 'K.csv' ) )
#       resultsNKI <- read.csv( paste0( 'labelresults', whichPipeline, 'N.csv' ) )
      resultsOasis <- read.csv( paste0( 'labelresults', whichPipeline, 'O.csv' ) )
      resultsCombined <- rbind( resultsOasis )

#       resultsCombined <- rbind( resultsIXI, resultsKirby, resultsNKI, resultsOasis )
      resultsCombined$SITE <- as.factor( resultsCombined$SITE )
      resultsCombined$SEX <- as.factor( resultsCombined$SEX )

      resultsCombined <- resultsCombined[which( resultsCombined$AGE >= 20 & resultsCombined$AGE <= 80 ),]

      corticalLabels <- tail( colnames( resultsCombined ), n = 62 )

      drops <- c( "ID", "SITE" )
      resultsCombined <- resultsCombined[, !( names( resultsCombined ) %in% drops )]

      trainingIndices <- createDataPartition( resultsCombined$SEX, p = trainingPortion, list = FALSE, times = 1 )
      trainingData <- resultsCombined[trainingIndices,]
      testingData <- resultsCombined[-trainingIndices,]

      brainAgeRF <- randomForest( AGE ~ ., data = trainingData,
                        na.action = na.omit, replace = FALSE, ntree = 200 )
      predictedAge <- predict( brainAgeRF, testingData )

#       regionalQuadraticTerms <- paste0( "I(", corticalLabels, collapse = "^2) * SEX + " )
#       myFormula <- as.formula( paste( "AGE ~ SEX + ", regionalTerms, " + ", regionalQuadraticTerms, "^2) + VOLUME ", sep = '' ) )
#       regionalTerms <- paste( corticalLabels, collapse = " * SEX + " )
#       myFormula <- as.formula( paste( "AGE ~ SEX + ", regionalTerms, " + VOLUME ", sep = '' ) )
#       brainAgeLM <- lm( myFormula, data = trainingData, na.action = na.omit )
#       predictedAge <- predict( brainAgeLM, testingData )

      rmse <- sqrt( mean( ( ( testingData$AGE - predictedAge )^2 ), na.rm = TRUE ) )

      oneData <- data.frame( Pipeline = whichPipeline, RMSE = rmse )
      resultsData <- rbind( resultsData, oneData )
      }
    }

  rmsePlot <- ggplot( resultsData, aes( x = RMSE, fill = Pipeline ) ) +
                      scale_y_continuous( "Density" ) +
                      scale_x_continuous( "RMSE", limits = c( 9, 15.5 ) ) +
                     geom_density( alpha = 0.5 )
  ggsave( filename = paste( "~/Desktop/rfRmse", p, ".pdf", sep = "" ), plot = rmsePlot, width = 6, height = 6, units = 'in' )
  cat( "Mean FS rmse = ", mean( resultsData$RMSE[which( resultsData$Pipeline == 'FreeSurfer' )], na.rm = TRUE ), "\n", sep = '' );
  cat( "Mean ANTs rmse = ", mean( resultsData$RMSE[which( resultsData$Pipeline == 'ANTs' )], na.rm = TRUE ), "\n", sep = '' );
  }

