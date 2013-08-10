# Load required libraries

library( MASS )
library( ggplot2 )

resultsIXI <- read.csv( '../labelresultsI.csv' )
resultsKirby <- read.csv( '../labelresultsK.csv' )
resultsNKI <- read.csv( '../labelresultsN.csv' )
resultsOasis <- read.csv( '../labelresultsO.csv' )

resultsCombined <- rbind( resultsIXI, resultsKirby, resultsNKI, resultsOasis )

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

myFormula <- as.formula( paste( "AGE ~ ", regionalTerms, " + SEX + SITE + VOLUME", sep = '' ) )
ageModel <- lm( myFormula, data = trainingData )
ageModelAIC <- stepAIC( ageModel, direction = "both", k = 4, trace = TRUE )
predictedAge <- predict( ageModelAIC, newdata = testingData, type = 'response' )

correlation <- cor( testingData$AGE, predictedAge )
meanerr <- mean( abs(  testingData$AGE - predictedAge ) )

annotation <- paste0( "r = ", signif( correlation, 2 ), ", mean error = ", signif( meanerr, 2 ), " years" )

print( summary( ageModelAIC ) )

plotData <- data.frame( TrueAge = testingData$AGE, PredictedAge = predictedAge )

brainAgePlot <- ggplot( plotData, aes( x = TrueAge, y = PredictedAge ) ) +
                stat_smooth( colour = "navyblue", formula = y ~ 1 + x, method = "lm",
                  size = 1, n = 1000, level = 0.95, se = TRUE, fullrange = TRUE, fill = 'black', alpha = 0.25 ) +
                geom_point( aes( colour = "darkred" ), size = 4, alpha = 0.5 ) +
                scale_x_continuous( "True age (years)", breaks = seq( 0, 100, by = 10 ), labels = seq( 0, 100, by = 10 ), limits = c( 0, 100 ) ) +
                scale_y_continuous( "Predicted age (years)", breaks = seq( 0, 100, by = 10 ), labels = seq( 0, 100, by = 10 ), limits = c( 0, 100 ) ) +
                scale_colour_manual( labels = annotation, values = "darkred" ) +
                theme( legend.justification = c( 1, 0 ), legend.position = c( 1, 0 ), legend.title = element_blank(), legend.background = element_blank(), legend.key = element_blank(), legend.text = element_text( colour = 'navyblue' ) ) +
                guides( colour = guide_legend( override.aes = list( shape = NA ) ), keywidth = 0 )
ggsave( filename = paste( "../ageRegressionPredict.pdf", sep = "" ), plot = brainAgePlot, width = 6, height = 6, units = 'in' )



# Call:
# lm(formula = AGE ~ LABEL_2 + LABEL_3 + LABEL_5 + LABEL_7 + LABEL_9 +
#     LABEL_11 + LABEL_12 + LABEL_13 + LABEL_15 + LABEL_16 + LABEL_21 +
#     LABEL_22 + LABEL_23 + LABEL_24 + LABEL_25 + LABEL_26 + LABEL_28 +
#     LABEL_31 + SITE + VOLUME, data = trainingData)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -43.663  -8.444   0.641   8.309  40.089
#
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)
# (Intercept)  9.198e+01  7.381e+00  12.463  < 2e-16 ***
# LABEL_2     -5.778e+00  2.132e+00  -2.711 0.006939 **
# LABEL_3     -4.746e+00  2.416e+00  -1.965 0.049967 *
# LABEL_5      4.226e+00  1.373e+00   3.079 0.002188 **
# LABEL_7      3.495e+00  1.509e+00   2.316 0.020947 *
# LABEL_9     -1.966e+01  3.004e+00  -6.545 1.43e-10 ***
# LABEL_11     5.998e+00  2.970e+00   2.019 0.043951 *
# LABEL_12     9.984e+00  2.862e+00   3.488 0.000527 ***
# LABEL_13    -6.731e+00  1.926e+00  -3.495 0.000514 ***
# LABEL_15     6.110e+00  2.741e+00   2.229 0.026231 *
# LABEL_16     6.378e+00  2.482e+00   2.570 0.010454 *
# LABEL_21    -1.505e+01  3.838e+00  -3.921 9.98e-05 ***
# LABEL_22    -8.432e+00  3.679e+00  -2.292 0.022322 *
# LABEL_23     1.469e+01  3.391e+00   4.332 1.77e-05 ***
# LABEL_24    -8.680e+00  3.367e+00  -2.578 0.010212 *
# LABEL_25    -1.634e+01  4.959e+00  -3.296 0.001048 **
# LABEL_26    -1.498e+01  5.047e+00  -2.968 0.003138 **
# LABEL_28     1.077e+01  3.656e+00   2.946 0.003361 **
# LABEL_31     1.361e+01  4.099e+00   3.320 0.000962 ***
# SITEHH      -1.036e+00  1.896e+00  -0.546 0.585114
# SITEIOP     -7.608e-01  2.844e+00  -0.267 0.789203
# SITEKirby   -1.121e+01  3.912e+00  -2.866 0.004328 **
# SITENKI     -3.149e+00  2.138e+00  -1.473 0.141377
# SITEOASIS    8.384e+00  1.970e+00   4.256 2.46e-05 ***
# VOLUME      -1.301e-05  4.463e-06  -2.915 0.003714 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 12.45 on 521 degrees of freedom
# Multiple R-squared:  0.6518,	Adjusted R-squared:  0.6358
# F-statistic: 40.64 on 24 and 521 DF,  p-value: < 2.2e-16
#
# age prediction error =  -0.870453
