library( ggplot2 )

whichThickness <- c( 'ANTs', 'FreeSurfer' )
sites <- c( 'IXI', 'Kirby', 'NKI', 'Oasis' )

boxPlotDataFrame <- data.frame( Site = character( 0 ), Pipeline = character( 0 ), ElapsedTime = numeric( 0 ) )

for( i in 1:length( sites ) )
  {
  for( j in 1:length( whichThickness ) )
    {
    times <- read.csv( paste0( 'Time', whichThickness[j], sites[i], '.csv' ) )

    site <- rep( sites[i], length( times$Id ) )
    pipeline <- rep( whichThickness[j], length( times$Id ) )

    boxPlotDataFrame <- rbind( boxPlotDataFrame, data.frame( Site = site, Pipeline = pipeline,
                               ElapsedTime = times$ElapsedTime ) )
    }
  }

myBoxPlot <- ggplot( boxPlotDataFrame, aes( x = Site, y = ElapsedTime ) ) +
                   geom_boxplot( aes( fill = Pipeline ) ) +
                   scale_x_discrete( "Cohort ", labels = sites ) +
                   scale_y_continuous( "Elapsed time (hours)" )
ggsave( filename = paste( "~/Desktop/Times.pdf", sep = "" ), plot = myBoxPlot, width = 6, height = 6, units = 'in' )
