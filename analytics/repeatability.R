library( psych )
library( ggplot2 )

kirby <- read.csv( "labelResultsK_pairwise.csv" )
oasis <- read.csv( "labelResultsO_pairwise.csv" )

corticalLabels <- c( "L occipital", "R occipital",
                     "L cingulate", "R cingulate",
                     "L insula",    "R insula",
                     "L temporal pole",   "R temporal pole",
                     "L superior temporal", "R superior temporal",
                     "L infero temporal", "R infero temporal",
                     "L parahippocampal", "R parahippocampal",
                     "L frontal pole",    "R frontal pole",
                     "L superior frontal","R superior frontal",
                     "L middle frontal",  "R middle frontal",
                     "L inferior",        "R inferior",
                     "L orbital frontal", "R orbital frontal",
                     "L precentral",      "R precentral",
                     "L superior parietal", "R superior parietal",
                     "L inferior parietal", "R inferior parietal",
                     "L postcentral",       "R postcentral" )

#
# The repeat scan is listed at every other row, i.e.
#
#    subject 1, scan 1, ...
#    subject 1, scan 2, ...
#    subject 2, scan 1, ...
#    subject 2, scan 2, ...
#

subjectID <- paste0( "S", c( 1:( 0.5 * ( nrow( kirby ) + nrow( oasis ) ) ) ) )
subjectID <- rbind( subjectID, subjectID )
subjectID <- as.vector( tmp )

scanOrder <- rep( c( "First", "Repeat" ), 0.5 * ( nrow( kirby ) + nrow( oasis ) ) )

repeatabilityDataFrame <- data.frame( rbind( kirby, oasis ) )
repeatabilityDataFrame$ID <- as.factor( subjectID )
repeatabilityDataFrame$ScanOrder <- as.factor( scanOrder )

iccData <- data.frame( First = as.vector( data.matrix( repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "First" ), 6:37] ) ),
                       Repeat = as.vector( data.matrix( repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "Repeat" ), 6:37] ) ) )

iccResults <- ICC( iccData )

# > ICC( iccData)
# Call: ICC(x = iccData)
#
# Intraclass correlation coefficients
#                          type  ICC  F  df1  df2 p lower bound upper bound
# Single_raters_absolute   ICC1 0.95 40 1279 1280 0        0.95        0.96
# Single_random_raters     ICC2 0.95 41 1279 1279 0        0.94        0.96
# Single_fixed_raters      ICC3 0.95 41 1279 1279 0        0.95        0.96
# Average_raters_absolute ICC1k 0.98 40 1279 1280 0        0.97        0.98
# Average_random_raters   ICC2k 0.98 41 1279 1279 0        0.97        0.98
# Average_fixed_raters    ICC3k 0.98 41 1279 1279 0        0.97        0.98
#
#  Number of subjects = 1280     Number of Judges =  2
#
# > iccResults$summary
#               Df Sum Sq Mean Sq F value   Pr(>F)
# subs        1279 1534.0  1.1994   41.45  < 2e-16 ***
# ind            1    1.2  1.2388   42.81 8.71e-11 ***
# Residuals   1279   37.0  0.0289

iccData$Difference <- iccData$First - iccData$Repeat
iccData$Average <- 0.5 * ( iccData$First + iccData$Repeat )

upperLineIntercept <- mean( iccData$Difference ) + 1.96 * sd( iccData$Difference )
lowerLineIntercept <- mean( iccData$Difference ) - 1.96 * sd( iccData$Difference )

repeatabilityPlot <- ggplot( iccData, aes( x = Average, y = Difference ) ) +
                geom_point( colour = "darkred", size = 4, alpha = 0.75 ) +
                geom_hline( aes( yintercept = upperLineIntercept ), colour = "navyblue", linetype = "dashed" ) +
                geom_hline( aes( yintercept = lowerLineIntercept ), colour = "navyblue", linetype = "dashed" ) +
                scale_x_continuous( "Average First/Repeat thickness (mm)" ) +
                scale_y_continuous( "First - Repeated thickness (mm)", breaks = seq( -1.5, 1.5, by = 0.5 ), labels = seq( -1.5, 1.5, by = 0.5 ), limits = c( -1.5, 1.5 ) ) +
                ggtitle( paste0( "Reproducibility (ICC = ", round( iccResults$results[2,2], digits = 2 ), ")" ) )
ggsave( filename = paste( "reproducibilityBA.pdf", sep = "" ), plot = repeatabilityPlot, width = 8, height = 6, units = 'in' )


repeatabilityError <- 100 * abs( repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "First" ), 6:37] -
  repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "Repeat" ), 6:37] ) /
  ( 0.5 * (  repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "First" ), 6:37] +
    repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "Repeat" ), 6:37] ) )

errorMeans <- c()
errorStds <- c()

for( i in 1:32 )
  {
  errorMeans[i] <- mean( repeatabilityError[,i] )
  errorStds[i] <- sd( repeatabilityError[,i] )
  }

for( i in 1:16 )
  {
  cat( corticalLabels[2 * i - 1], " & $",
    round( errorMeans[2 * i - 1], digits = 2 ), " \\pm ", round( errorStds[2 * i - 1], digits = 2 ), "$ & $",
    round( errorMeans[2 * i], digits = 2 ), " \\pm ", round( errorStds[2 * i], digits = 2 ), "$\\\\\n",
    sep = '' );
  }



# qvalues <- c()
# meanDifferences <- c()
# for( i in 6:37 )
#   {
#   myTest <- t.test( repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "First" ), i],
#                     repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "Repeat" ), i],
#                     paired = TRUE, conf.int = TRUE )
#   qvalues[i-5] <- myTest$p.value
#   meanDifferences[i-5] <- myTest$estimate
#   }
# qvalues <- p.adjust( qvalues, method = "fdr" )
#
# for( i in 1:32 )
#   {
#   cat( corticalLabels[i], ": mean of the differences = ", meanDifferences[i],
#        " (q-value = ", qvalues[i], ")\n", sep = '' )
#   }

