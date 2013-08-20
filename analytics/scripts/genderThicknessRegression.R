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

myFormula <- as.formula( paste( "SEX ~ AGE + ", regionalTerms, " * AGE + VOLUME ", sep = '' ) )
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

# Call:
# glm(formula = SEX ~ AGE + LABEL_1 + LABEL_2 + LABEL_4 + LABEL_7 +
#     LABEL_8 + LABEL_10 + LABEL_12 + LABEL_14 + LABEL_15 + LABEL_19 +
#     LABEL_22 + LABEL_24 + LABEL_27 + LABEL_29 + LABEL_32 + VOLUME +
#     AGE:LABEL_7 + AGE:LABEL_24 + AGE:LABEL_27 + AGE:LABEL_29,
#     family = "binomial", data = trainingData)
#
# Deviance Residuals:
#     Min       1Q   Median       3Q      Max
# -2.4730  -0.6031   0.1456   0.5776   2.6735
#
# Coefficients:
#                Estimate Std. Error z value Pr(>|z|)
# (Intercept)   1.294e+01  4.088e+00   3.166 0.001546 **
# AGE           5.494e-02  6.356e-02   0.864 0.387410
# LABEL_1       1.566e+00  6.861e-01   2.282 0.022490 *
# LABEL_2      -1.966e+00  7.218e-01  -2.724 0.006455 **
# LABEL_4       1.448e+00  5.927e-01   2.444 0.014539 *
# LABEL_7      -8.254e-01  7.584e-01  -1.088 0.276453
# LABEL_8      -1.246e+00  4.948e-01  -2.518 0.011803 *
# LABEL_10     -1.718e+00  7.308e-01  -2.350 0.018749 *
# LABEL_12      1.962e+00  5.172e-01   3.794 0.000148 ***
# LABEL_14     -2.101e+00  5.163e-01  -4.070 4.71e-05 ***
# LABEL_15      1.539e+00  5.053e-01   3.046 0.002315 **
# LABEL_19     -1.473e+00  5.730e-01  -2.570 0.010172 *
# LABEL_22     -1.690e+00  6.659e-01  -2.538 0.011148 *
# LABEL_24      2.589e+00  1.297e+00   1.997 0.045879 *
# LABEL_27     -3.615e+00  1.649e+00  -2.193 0.028309 *
# LABEL_29      4.442e+00  1.333e+00   3.332 0.000862 ***
# LABEL_32      3.322e+00  8.151e-01   4.075 4.61e-05 ***
# VOLUME       -1.528e-05  1.432e-06 -10.671  < 2e-16 ***
# AGE:LABEL_7   3.005e-02  1.410e-02   2.130 0.033145 *
# AGE:LABEL_24 -4.548e-02  2.217e-02  -2.052 0.040212 *
# AGE:LABEL_27  1.012e-01  3.241e-02   3.123 0.001792 **
# AGE:LABEL_29 -9.585e-02  2.791e-02  -3.434 0.000594 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# (Dispersion parameter for binomial family taken to be 1)
#
#     Null deviance: 719.91  on 520  degrees of freedom
# Residual deviance: 410.44  on 499  degrees of freedom
# AIC: 454.44
#
# Number of Fisher Scoring iterations: 6
#
# sex prediction error =  0.05347594
