library( kernlab )
library( caret )
library( randomForest )
library( ggplot2 )
library( grid )

nPermutations <- 500

trainingPortions <- c( 0.5 ) # seq( 0.05, 0.25, by = 0.01 ); #

antsRmse <- c()
fsRmse <- c()

count <- 1
for( p in trainingPortions )
  {
  trainingPortion <- p
  cat( "trainingPortion = ", trainingPortion, "\n", sep = '' )

  resultsData <- data.frame( Pipeline = character( 0 ), Correlation = numeric( 0 ), MeanErrors = numeric( 0 ) )

  forestImp <- list()
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
      resultsCombined$SEX <- as.factor( resultsCombined$SEX )

      corticalLabels <- tail( colnames( resultsCombined ), n = 62 )

      drops <- c( "ID", "SITE" )
      resultsCombined <- resultsCombined[, !( names( resultsCombined ) %in% drops )]

      trainingIndices <- createDataPartition( resultsCombined$SEX, p = trainingPortion, list = FALSE, times = 1 )
      trainingData <- resultsCombined[trainingIndices,]
      testingData <- resultsCombined[-trainingIndices,]

    #   brainAgeVM <- rvm( AGE ~ .,
    #                       data = trainingData, type = "regression",
    #                       kernel = "rbfdot",
    #                       verbosity = 1, tol = .Machine$double.eps,
    #                       minmaxdiff = 1e-3, cross = 0, fit = TRUE,
    #                       na.action = na.omit , iterations = 1000 )
    #   predictedAge <- predict( brainAgeVM, testingData )

       brainAgeRF <- randomForest( AGE ~ ., data = trainingData, importance = TRUE,
                    na.action = na.omit, replace = FALSE, ntree = 200 )
       if( n == 1 )
         {
         if( whichPipeline == 'ANTs' )
           {
           forestImp[[1]] <- importance( brainAgeRF, type = 1 )
           } else {
           forestImp[[2]] <- importance( brainAgeRF, type = 1 )
           }
         } else {
         if( whichPipeline == 'ANTs' )
           {
           forestImp[[1]] <- forestImp[[1]] + importance( brainAgeRF, type = 1 )
           } else {
           forestImp[[2]] <- forestImp[[2]] + importance( brainAgeRF, type = 1 )
           }
         }
       predictedAge <- predict( brainAgeRF, testingData )

#       regionalTerms <- paste( corticalLabels, collapse = " * SEX + " )
#       myFormula <- as.formula( paste( "AGE ~ SEX + ", regionalTerms, " + VOLUME ", sep = '' ) )
#       brainAgeLM <- lm( myFormula, data = trainingData, na.action = na.omit )
#       predictedAge <- predict( brainAgeLM, testingData )

      rmse <- sqrt( mean( ( ( testingData$AGE - predictedAge )^2 ), na.rm = TRUE ) )

      oneData <- data.frame( Pipeline = whichPipeline, RMSE = rmse )
      resultsData <- rbind( resultsData, oneData )
      }
    }

  for( n in 1:2 )
    {
#     forestImp[[n]] <- forestImp[[n]] / nPermutations

    forestImp.df <- data.frame( Statistic = names( forestImp[[n]][,1] ), Importance = as.numeric( forestImp[[n]][,1] )  )
    forestImp.df <- forestImp.df[order( forestImp.df$Importance ),]

    forestImp.df$Statistic <- factor( x = forestImp.df$Statistic, levels = forestImp.df$Statistic )
    levels( forestImp.df$Statistic )[which( levels( forestImp.df$Statistic ) == 'SEX' )] <- 'Gender'
    levels( forestImp.df$Statistic )[which( levels( forestImp.df$Statistic ) == 'VOLUME' )] <- 'Volume'

    vPlot <- ggplot( data = forestImp.df, aes( x = Importance, y = Statistic ) ) +
             geom_point( aes( color = Importance ) ) +
             ylab( "" ) +
             scale_x_continuous( "MeanDecreaseAccuracy", limits = c( 0, 16 ) ) +
             scale_color_continuous( low = "navyblue", high = "darkred" ) +
             theme( axis.text.y = element_text( size = 8 ) ) +
             theme( plot.margin = unit( c( 0.1, 0.1, 0.1, -0.5 ), "cm" ) ) +
             theme( axis.title = element_text( size = 9 ) ) +
             theme( legend.position = "none" )

    ggsave( file = paste( "~/Desktop/importanceCombined", thicknessTypes[n], p, ".pdf", sep = "" ), plot = vPlot, width = 3, height = 8 )
    }

  rmsePlot <- ggplot( resultsData, aes( x = RMSE, fill = Pipeline ) ) +
                      scale_y_continuous( "Density" ) +
                      scale_x_continuous( "RMSE" ) +
                     geom_density( alpha = 0.5 )
  ggsave( filename = paste( "~/Desktop/xrfRmse", p, ".pdf", sep = "" ), plot = rmsePlot, width = 6, height = 6, units = 'in' )

  cat( "Mean FS rmse = ", mean( resultsData$RMSE[which( resultsData$Pipeline == 'FreeSurfer' )], na.rm = TRUE ), "\n", sep = '' );
  cat( "Mean ANTs rmse = ", mean( resultsData$RMSE[which( resultsData$Pipeline == 'ANTs' )], na.rm = TRUE ), "\n", sep = '' );

  antsRmse[count] <- mean( resultsData$RMSE[which( resultsData$Pipeline == 'ANTs' )], na.rm = TRUE )
  fsRmse[count] <- mean( resultsData$RMSE[which( resultsData$Pipeline == 'FreeSurfer' )], na.rm = TRUE )
  count <- count + 1

# for( i in 1:length( forestImp[[1]] ) )
#   {
#   cat( rownames( forestImp[[1]] )[i], forestImp[[1]][i], forestImp[[2]][i],  sep = ',' )
#   cat( "\n" )
#   }

  }




# cat( fsRmse, "\n" )
# cat( antsRmse, "\n" )
#
# resultsData <- data.frame( TrainingPortions <- c( trainingPortions, trainingPortions ),
#                            RMSE = c( antsRmse, fsRmse ),
#                            Pipeline = c( rep( 'ANTs', length( antsRmse ) ),  rep( 'FreeSurfer', length( fsRmse ) ) ) )
# rmsePlot <- ggplot( resultsData, aes( x = TrainingPortions, y = RMSE, colour = Pipeline ) ) +
#                     geom_line( size = 1 ) +
#                     scale_y_continuous( "Mean RMSE" ) +
#                     scale_x_continuous( "Training portion" )
# ggsave( filename = paste( "~/Desktop/Plots/lmTotalRmse.pdf", sep = "" ), plot = rmsePlot, width = 6, height = 6, units = 'in' )







