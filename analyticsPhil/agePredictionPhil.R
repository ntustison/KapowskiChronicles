library( caret )
library( randomForest )
library( ggplot2 )

# Test of pull request
# https://github.com/ANTsX/ANTs/pull/514

nPermutations <- 1000

trainingPortions <- c( 0.5 )

resultsIXI <- read.csv( 'labelresultsIXI.csv' )
resultsKirby <- read.csv( 'labelresultsKirby.csv' )
resultsNKI <- read.csv( 'labelresultsNKI.csv' )

thicknessTypes <- c( 'Ants', 'Current', 'Phil' )

for( p in trainingPortions )
  {
  trainingPortion <- p
  cat( "trainingPortion = ", trainingPortion, "\n", sep = '' )

  resultsData <- data.frame( Pipeline = character( 0 ), Correlation = numeric( 0 ), MeanErrors = numeric( 0 ) )

  for( n in seq( 1, nPermutations, by = 1 ) )
    {
    cat( "  Permutation ", n, "\n", sep = '' )

    for( whichPipeline in thicknessTypes )
      {
      resultsIXI.pipeline <- resultsIXI[which( resultsIXI$PIPELINE == whichPipeline ),]
      resultsKirby.pipeline <- resultsKirby[which( resultsKirby$PIPELINE == whichPipeline ),]
      resultsNKI.pipeline <- resultsNKI[which( resultsNKI$PIPELINE == whichPipeline ),]
      
      resultsCombined <- rbind( resultsIXI.pipeline, resultsKirby.pipeline, resultsNKI.pipeline )
      resultsCombined$SITE <- as.factor( resultsCombined$SITE )
      resultsCombined$SEX <- as.factor( resultsCombined$SEX )

      resultsCombined <- resultsCombined[which( resultsCombined$AGE >= 20 & resultsCombined$AGE <= 80 ),]

      corticalLabels <- tail( colnames( resultsCombined ), n = 62 )

      drops <- c( "ID", "SITE", "PIPELINE" )
      resultsCombined <- resultsCombined[, !( names( resultsCombined ) %in% drops )]

      trainingIndices <- createDataPartition( resultsCombined$SEX, p = trainingPortion, list = FALSE, times = 1 )
      trainingData <- resultsCombined[trainingIndices,]
      testingData <- resultsCombined[-trainingIndices,]

      brainAgeRF <- randomForest( AGE ~ ., data = trainingData,
                        na.action = na.omit, replace = FALSE, ntree = 200 )
      predictedAge <- predict( brainAgeRF, testingData )

      rmse <- sqrt( mean( ( ( testingData$AGE - predictedAge )^2 ), na.rm = TRUE ) )

      oneData <- data.frame( Pipeline = whichPipeline, RMSE = rmse )
      resultsData <- rbind( resultsData, oneData )
      }
    }

  rmsePlot <- ggplot( resultsData, aes( x = RMSE, fill = Pipeline ) ) +
                      scale_y_continuous( "Density" ) +
                      scale_x_continuous( "RMSE", limits = c( 6, 14. ) ) +
                     geom_density( alpha = 0.5 )
  ggsave( filename = paste( "~/Desktop/rfRmse", p, ".pdf", sep = "" ), plot = rmsePlot, width = 6, height = 6, units = 'in' )
  cat( "Mean Ants rmse = ", mean( resultsData$RMSE[which( resultsData$Pipeline == 'Ants' )], na.rm = TRUE ), "\n", sep = '' );
  cat( "Mean Current rmse = ", mean( resultsData$RMSE[which( resultsData$Pipeline == 'Current' )], na.rm = TRUE ), "\n", sep = '' );
  cat( "Mean Phil rmse = ", mean( resultsData$RMSE[which( resultsData$Pipeline == 'Phil' )], na.rm = TRUE ), "\n", sep = '' );
  }

