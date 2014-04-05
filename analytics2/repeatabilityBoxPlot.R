library( psych )
library( ggplot2 )

kirby <- list()
oasis <- list()

kirby[[1]] <- read.csv( "labelresultsPairwiseK.csv" )
oasis[[1]] <- read.csv( "labelresultsPairwiseO.csv" )
kirby[[2]] <- read.csv( "labelresultsPairwiseFreesurferK.csv" )
oasis[[2]] <- read.csv( "labelresultsPairwiseFreesurferO.csv" )

whichThickness <- c( 'ANTs', 'FreeSurfer' )
sites <- c( 'Kirby', 'Oasis' )

#
# The repeat scan is listed at every other row, i.e.
#
#    subject 1, scan 1, ...
#    subject 1, scan 2, ...
#    subject 2, scan 1, ...
#    subject 2, scan 2, ...
#

boxPlotDataFrame <- data.frame( Site = character( 0 ), Pipeline = character( 0 ),
                                Hemisphere = character( 0 ), CorticalRegion = character( 0 ),
                                Error = numeric( 0 ) )

for( i in 1:length( whichThickness ) )
  {
  for( j in 1:2 )
    {
    whichSite <- kirby[[i]]
    if( j == 2 )
      {
      whichSite <- oasis[[i]]
      }

    subjectID <- paste0( "S", c( 1:( 0.5 * ( nrow( whichSite ) ) ) ) )
    repeatabilityDataFrame <- data.frame( whichSite )
    scanOrder <- rep( c( "First", "Repeat" ), 0.5 * ( nrow( whichSite ) ) )

    subjectID <- rbind( subjectID, subjectID )


    corticalLabels <- tail( colnames( repeatabilityDataFrame ), n = 62 )
    numberOfLabels <- length( corticalLabels )
    beginIndex <- length( colnames( repeatabilityDataFrame ) ) - 62 + 1
    endIndex <- length( colnames( repeatabilityDataFrame ) )
    indices <- beginIndex:endIndex

    repeatabilityDataFrame$ID <- as.factor( subjectID )
    repeatabilityDataFrame$ScanOrder <- as.factor( scanOrder )

    repeatabilityError <- 100 * abs( repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "First" ), indices] -
      repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "Repeat" ), indices] ) /
      ( 0.5 * (  repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "First" ), indices] +
        repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "Repeat" ), indices] ) )

    for( k in 1:length( corticalLabels ) )
      {
      r <- repeatabilityError[,k]
      nR <- length( r )

      hemisphere <- rep( 'Right', nR )
      if( k < 32 )
        {
        hemisphere <- rep( 'Left', nR )
        }
      site <- rep( 'Oasis', nR )
      if( j == 1 )
        {
        site <- rep( 'Kirby', nR )
        }
      pipeline <- rep( whichThickness[i], nR )

      roiVector <- strsplit( colnames( repeatabilityError )[k], ".", fixed = TRUE )
      roi <- paste( roiVector[2:length( roiVector )], collapse = ' ' )

      boxPlotDataFrame <- rbind( boxPlotDataFrame, data.frame( Site = site, Pipeline = pipeline,
                                 Hemisphere = hemisphere, CorticalRegion = roi, Error = r ) )
      }
    cat( whichThickness[i], " (", sites[j], ") ",
      "mean repeatability error:  ", mean( as.matrix( repeatabilityError ) ),
      " (", sd( as.matrix( repeatabilityError ) ), ")",
      "\n", sep = ""  )
    }
  }


boxPlotLeftDataFrame <- boxPlotDataFrame[which( boxPlotDataFrame$Hemisphere == 'Left' ),]
boxPlotRightDataFrame <- boxPlotDataFrame[which( boxPlotDataFrame$Hemisphere == 'Right' ),]

myBoxPlotLeft <- ggplot( boxPlotLeftDataFrame, aes( x = CorticalRegion, y = Error ) ) +
                   geom_boxplot( aes( fill = Pipeline ) ) +
                   scale_x_discrete( "Cortical region (left hemisphere)", labels = 1:32 ) +
                   scale_y_continuous( "Error (% variability)", limits = c( 0, 20 ) )
ggsave( filename = paste( "~/Desktop/repeatabilityLeft.pdf", sep = "" ), plot = myBoxPlotLeft, width = 12, height = 4, units = 'in' )

myBoxPlotRight <- ggplot( boxPlotRightDataFrame, aes( x = CorticalRegion, y = Error ) ) +
                   geom_boxplot( aes( fill = Pipeline ) ) +
                   scale_x_discrete( "Cortical region (right hemisphere)", labels = 1:32 ) +
                   scale_y_continuous( "Error (% variability)", limits = c( 0, 20 ) )
ggsave( filename = paste( "~/Desktop/repeatabilityRight.pdf", sep = "" ), plot = myBoxPlotRight, width = 12, height = 4, units = 'in' )
