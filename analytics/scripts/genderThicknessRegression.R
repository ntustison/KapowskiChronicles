# Load required libraries

library( MASS )
library( popbio )
library( gdata )
library( ROCR )
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
# Do sex prediction
#   Include all AGE/REGION interaction terms
#
#

regionalTerms <- paste( names( resultsCombined )[6:37] , collapse = " * AGE + " )

myFormula <- as.formula( paste( "SEX ~ ", regionalTerms, " * AGE + VOLUME + (AGE^2) ", sep = '' ) )
sexModel <- glm( myFormula, data = trainingData, family = "binomial" )
sexModelAIC <- stepAIC( sexModel, direction = "both", k = 4, trace = TRUE )
sexPrediction <- predict( sexModelAIC, newdata = testingData, type = 'response' )

sexPredictionError <- mean( round( sexPrediction ) - testingData$SEX, na.rm = TRUE )

print( summary( sexModelAIC ) )
cat( "sex prediction error = ", sexPredictionError );

# plot results

# sexPredictionTF <- sexPrediction
#
# sexPredictionTF[which( round( sexPrediction ) == 1 )] <- TRUE  # MALE
# sexPredictionTF[which( round( sexPrediction ) == 0 )] <- FALSE # FEMALE
#
# pdf( "../figs/ageGenderThicknessRegression_brainVolume.pdf", width = 8, height = 6 )
# logi.hist.plot( testingData$VOLUME, sexPredictionTF, boxp = FALSE, type = "hist", col = "gray",
#   mainlabel = "", xlabel = expression( paste( "Brain Volume (mm"^"3", ")", sep = '' ) ),
#   ylabel = "Probability of Gender (Female)", ylabel2 = "Volume Frequency" )
# dev.off()
#
# sexCoefficients <- summary( sexModelAIC )$coefficients
#
# for( i in 1:length( rownames( sexCoefficients ) ) )
#   {
#   name <- rownames( sexCoefficients )[i]
#   if( startsWith( name, "LABEL_" ) )
#     {
#     index <- as.integer( gsub( "LABEL_", replacement = "", name ) )
#     pdf( paste0( "../figs/ageGenderThicknessRegression_", name, ".pdf" ), width = 8, height = 6 )
#     logi.hist.plot( testingData[,index+5], sexPredictionTF, boxp = FALSE, type = "hist", col = "gray",
#       mainlabel = corticalLabels[index], xlabel = expression( paste( "Thickness (mm"^"3", ")", sep = '' ) ),
#       ylabel = "Probability of Gender (Female)",ylabel2 = "Thickness Frequency" )
#     dev.off()
#     }
#   }

# plot ROC curve

sexPredictionROC <- prediction( sexPrediction, testingData$SEX )
sexPerformanceROC <- performance( sexPredictionROC, "tpr", "fpr" )

#
# http://www.rrandomness.com/r/simple-roc-plots-with-ggplot2-part-1/
# Hanley JA, McNeil BJ. The meaning and use of the area under a receiver operating
# characteristic (ROC) curve. Radiology 1982;143:29-36.
#
# http://gim.unmc.edu/dxtests/ROC3.htm
#
# The accuracy of the test depends on how well the test separates the group being
# tested into those with and without the disease in question. Accuracy is measured
# by the area under the ROC curve. An area of 1 represents a perfect test; an area
# of .5 represents a worthless test. A rough guide for classifying the accuracy of
# a diagnostic test is the traditional academic point system:
#
# .90-1 = excellent (A)
# .80-.90 = good (B)
# .70-.80 = fair (C)
# .60-.70 = poor (D)
# .50-.60 = fail (F)
#

auc <- performance( sexPredictionROC, "auc" )@y.values[[1]]
q1 <- auc / ( 2 - auc )
q2 <- ( 2 * auc^2 ) / ( 1 + auc )

neg <- sexPrediction[testingData$SEX == 1]
pos <- sexPrediction[testingData$SEX == 0]

se.auc <- sqrt( ( ( auc * ( 1 - auc ) ) + ( ( length( sexPrediction ) - 1 ) *
  ( q1 - auc^2) ) + ( ( length( neg ) - 1 ) *( q2 - auc^2 ) ) ) / ( length( pos ) * length( neg ) ) )
ci.upper <- auc + (se.auc * 0.96)
ci.lower <- auc - (se.auc * 0.96)

sexPlotData <- data.frame( x = sexPerformanceROC@x.values[[1]], y = sexPerformanceROC@y.values[[1]] )

annotation <- paste0( "AUC = ", signif( auc, 2 ), " (95%CI: ", signif( ci.lower, 2 ), " - ", signif( ci.upper, 2 ), ")" )

sexPlot <- ggplot( sexPlotData, aes( x = x, y = y ) ) +
      geom_line( aes( colour = "navyblue" ), size = 1.0 ) +
      geom_abline ( intercept = 0, slope = 1, colour = "darkred", linetype = "dashed", size = 1.0 ) +
      scale_x_continuous( "False positive rate (1-specificity)" ) +
      scale_y_continuous( "True positive rate (sensitivity)" ) +
      scale_colour_manual( labels = annotation, values = "navyblue" ) +
      theme( legend.justification = c( 1, 0 ), legend.position = c( 1, 0 ), legend.title = element_blank(), legend.background = element_blank(), legend.key = element_blank(), legend.text = element_text( colour = 'navyblue' ) ) +
      guides( colour = guide_legend( override.aes = list( shape = NA ) ), keywidth = 0 )
ggsave( filename = paste( "../figs/sexPlot.pdf", sep = "" ), plot = sexPlot, width = 6, height = 6, units = 'in' )

