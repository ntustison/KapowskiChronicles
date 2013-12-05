library( psych )
library( ggplot2 )
library( gridExtra )

oasis <- read.csv( "../labelresultsPairwiseFreesurferO.csv" )

#
# The repeat scan is listed at every other row, i.e.
#
#    subject 1, scan 1, ...
#    subject 1, scan 2, ...
#    subject 2, scan 1, ...
#    subject 2, scan 2, ...
#

subjectID <- paste0( "S", c( 1:( 0.5 * ( nrow( oasis ) ) ) ) )
scanOrder <- rep( c( "First", "Repeat" ), 0.5 * ( nrow( oasis ) ) )

repeatabilityDataFrame <- data.frame( oasis )
repeatabilityDataFrame$ID <- as.factor( subjectID )
repeatabilityDataFrame$ScanOrder <- as.factor( scanOrder )
repeatabilityDataFrame$SEX <- as.factor( repeatabilityDataFrame$SEX )
repeatabilityDataFrame$SITE <- as.factor( repeatabilityDataFrame$SITE )

numberOfRegions <- ncol( oasis ) - 5;

iccData <- data.frame( First = as.vector( data.matrix( repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "First" ), 6:ncol( oasis )] ) ),
                       Repeat = as.vector( data.matrix( repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "Repeat" ), 6:ncol( oasis )] ) ) )

iccData <- na.omit( iccData )
iccResults <- ICC( iccData )

# Intraclass correlation coefficients
#                          type  ICC  F  df1  df2 p lower bound
# Single_raters_absolute   ICC1 0.83 11 1323 1324 0        0.82
# Single_random_raters     ICC2 0.84 12 1323 1323 0        0.80
# Single_fixed_raters      ICC3 0.85 12 1323 1323 0        0.83
# Average_raters_absolute ICC1k 0.91 11 1323 1324 0        0.90
# Average_random_raters   ICC2k 0.91 12 1323 1323 0        0.89
# Average_fixed_raters    ICC3k 0.92 12 1323 1323 0        0.91
#                         upper bound
# Single_raters_absolute         0.85
# Single_random_raters           0.87
# Single_fixed_raters            0.86
# Average_raters_absolute        0.92
# Average_random_raters          0.93
# Average_fixed_raters           0.93
#
#  Number of subjects = 1324     Number of Judges =  2


iccData$Difference <- iccData$First - iccData$Repeat
iccData$Average <- 0.5 * ( iccData$First + iccData$Repeat )

upperLineIntercept <- mean( iccData$Difference ) + 1.96 * sd( iccData$Difference )
lowerLineIntercept <- mean( iccData$Difference ) - 1.96 * sd( iccData$Difference )

repeatabilityPlot <- ggplot( iccData, aes( x = Average, y = Difference ) ) +
                geom_point( colour = "darkred", size = 4, alpha = 0.5 ) +
                geom_hline( aes( yintercept = upperLineIntercept ), colour = "navyblue", linetype = "dashed" ) +
                geom_hline( aes( yintercept = lowerLineIntercept ), colour = "navyblue", linetype = "dashed" ) +
                scale_x_continuous( "Average First/Repeat thickness (mm)" ) +
                scale_y_continuous( "First - Repeated thickness (mm)", breaks = seq( -1.5, 1.5, by = 0.5 ), labels = seq( -1.5, 1.5, by = 0.5 ), limits = c( -1.5, 1.5 ) ) +
                ggtitle( paste0( "Reproducibility (ICC = ", round( iccResults$results[2,2], digits = 2 ), ")" ) )
ggsave( filename = paste( "reproducibilityBA.pdf", sep = "" ), plot = repeatabilityPlot, width = 8, height = 6, units = 'in' )

################################# do individual subject plots



# plotData <- data.frame( Subject = factor( rep.int( subjectID, numberOfRegions ), levels = subjectID ),
#                         First = as.vector( data.matrix( repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "First" ), 6:ncol( oasis )] ) ),
#                         Repeat = as.vector( data.matrix( repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "Repeat" ), 6:ncol( oasis )] ) ) )
#
# myPlot <- ggplot( plotData, aes( x = First, y = Repeat ) ) +
#           geom_point( aes( colour = Subject ), alpha = 0.5, size = 3 )

colors <- rainbow( length( subjectID ) )


subjectPlots <- list()
for( i in 1:length( subjectID ) )
  {
		plotData <- data.frame( First = as.vector( data.matrix( repeatabilityDataFrame[2*i-1,6:ncol( oasis )] ) ),
																										Repeat = as.vector( data.matrix( repeatabilityDataFrame[2*i,6:ncol( oasis )] ) ) )

		subjectPlots[[i]] <- ggplot( plotData, aes( x = First, y = Repeat ) ) +
												             geom_point( colour = colors[i], alpha = 1.0, size = 5 )
#   ggsave( filename = paste0( "~/Desktop/reproducibility", i, ".pdf" ), plot = myPlot, width = 8, height = 8, units = 'in' )
  }

pdf( paste0( "allSubjects.pdf" ), width = 20, height = 16 )
do.call("grid.arrange", c(subjectPlots, ncol=5))
dev.off()



