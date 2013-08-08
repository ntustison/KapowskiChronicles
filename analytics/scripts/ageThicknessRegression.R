# Load required libraries

library( MASS )
library( ggplot2 )

resultsIXI <- read.csv( '../labelresultsI.csv' )
resultsKirby <- read.csv( '../labelresultsK.csv' )
resultsNKI <- read.csv( '../labelresultsN.csv' )
resultsOasis <- read.csv( '../labelresultsO.csv' )

resultsCombined <- rbind( resultsIXI, resultsKirby, resultsNKI, resultsOasis )
resultsCombined$SEX <- resultsCombined$SEX - 1

trainingPercentage <- 0.5
partitioning <- runif( nrow( resultsCombined ) )
trainingIndices <- which( partitioning < trainingPercentage )
testingIndices <- which( partitioning >= trainingPercentage )
trainingData <- resultsCombined[trainingIndices,];
testingData <- resultsCombined[testingIndices,];

#
# Do gender prediction
#   Include all AGE/REGION interaction terms
#
#

regionalTerms <- paste( names( resultsCombined )[6:37] , collapse = " + " )

myFormula <- as.formula( paste( "AGE ~ ", regionalTerms, " + SITE + VOLUME", sep = '' ) )
ageModel <- lm( myFormula, data = trainingData )
ageModelAIC <- stepAIC( ageModel, direction = "both", k = 4, trace = TRUE )
agePrediction <- predict( ageModelAIC, newdata = testingData, type = 'response' )

agePredictionError <- mean( round( agePrediction ) - testingData$AGE, na.rm = TRUE )

print( summary( ageModelAIC ) )
cat( "age prediction error = ", agePredictionError );

plotData <- data.frame( TrueAge = testingData$AGE, PredictedAge = agePrediction )

brainAgePlot <- ggplot( plotData, aes( x = TrueAge, y = PredictedAge ) ) +
                stat_smooth( colour = "navyblue", formula = y ~ 1 + x, method = "lm",
                  size = 1, n = 1000, level = 0.95, se = TRUE, fullrange = TRUE, fill = 'black', alpha = 0.25 ) +
                geom_point( colour = "darkred", size = 4, alpha = 0.75 ) +
                scale_x_continuous( "True age (years)", breaks = seq( 0, 100, by = 10 ), labels = seq( 0, 100, by = 10 ), limits = c( 0, 100 ) ) +
                scale_y_continuous( "Predicted age (years)", breaks = seq( 0, 100, by = 10 ), labels = seq( 0, 100, by = 10 ), limits = c( 0, 100 ) )
ggsave( filename = paste( "../ageRegressionPredict.pdf", sep = "" ), plot = brainAgePlot, width = 6, height = 6, units = 'in' )
