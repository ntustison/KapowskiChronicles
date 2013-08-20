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
# lm(formula = AGE ~ LABEL_2 + LABEL_3 + LABEL_7 + LABEL_9 + LABEL_12 +
#     LABEL_13 + LABEL_15 + LABEL_16 + LABEL_21 + LABEL_23 + LABEL_25 +
#     LABEL_26 + LABEL_27 + LABEL_28 + LABEL_31 + SITE, data = trainingData)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -28.402  -7.868   0.209   7.055  32.831
#
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept) 109.8142     6.8401  16.055  < 2e-16 ***
# LABEL_2     -12.5064     2.2169  -5.641 2.82e-08 ***
# LABEL_3     -22.2376     2.5452  -8.737  < 2e-16 ***
# LABEL_7       3.8367     1.5585   2.462  0.01416 *
# LABEL_9     -12.6620     2.2789  -5.556 4.48e-08 ***
# LABEL_12      7.1609     2.1789   3.287  0.00109 **
# LABEL_13     -5.5380     1.9021  -2.912  0.00376 **
# LABEL_15      6.1994     2.4185   2.563  0.01066 *
# LABEL_16      7.0774     2.1533   3.287  0.00108 **
# LABEL_21    -16.8904     2.7239  -6.201 1.17e-09 ***
# LABEL_23     14.0643     2.7166   5.177 3.26e-07 ***
# LABEL_25     -8.9672     4.1883  -2.141  0.03275 *
# LABEL_26    -11.4948     4.3091  -2.668  0.00789 **
# LABEL_27     17.5986     3.8935   4.520 7.72e-06 ***
# LABEL_28     -8.5283     3.9651  -2.151  0.03196 *
# LABEL_31      7.7227     3.3973   2.273  0.02343 *
# SITEHH       -0.5685     1.6813  -0.338  0.73539
# SITEIOP       6.0230     2.7156   2.218  0.02701 *
# SITEKirby     1.2137     4.8197   0.252  0.80129
# SITENKI       6.0722     2.1591   2.812  0.00511 **
# SITEOASIS    20.5323     2.2535   9.111  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 11.65 on 503 degrees of freedom
# Multiple R-squared:  0.6714,	Adjusted R-squared:  0.6584
# F-statistic: 51.39 on 20 and 503 DF,  p-value: < 2.2e-16
#
# Warning messages:
# 1: Removed 2 rows containing missing values (stat_smooth).
# 2: Removed 2 rows containing missing values (geom_point).
# [ntustison@Nietzschean-Numerics Mon Aug 19 13:46:40] $
