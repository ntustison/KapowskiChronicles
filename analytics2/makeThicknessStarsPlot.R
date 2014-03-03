# Load required libraries

library( graphics )
library( nortest )

whichPipeline <- c( "ANTs", "FreeSurfer" )

resultsANTsIXI <- read.csv( "labelresultsANTsI.csv" )
resultsANTsKirby <- read.csv( "labelresultsANTsK.csv" )
resultsANTsNKI <- read.csv( "labelresultsANTsN.csv" )
resultsANTsOasis <- read.csv( "labelresultsANTsO.csv" )
resultsANTsCombined <- rbind( resultsANTsIXI, resultsANTsKirby, resultsANTsNKI, resultsANTsOasis )

tag <- c()
for( i in 1:nrow( resultsANTsCombined ) )
  {
  tag <- append( tag, paste0( resultsANTsCombined$ID[i], resultsANTsCombined$SITE[i] ) );
  }
resultsANTsCombined$TAG <- tag;

resultsFreeSurferIXI <- read.csv( "labelresultsFreeSurferI.csv" )
resultsFreeSurferKirby <- read.csv( "labelresultsFreeSurferK.csv" )
resultsFreeSurferNKI <- read.csv( "labelresultsFreeSurferN.csv" )
resultsFreeSurferOasis <- read.csv( "labelresultsFreeSurferO.csv" )
resultsFreeSurferCombined <- rbind( resultsFreeSurferIXI, resultsFreeSurferKirby, resultsFreeSurferNKI, resultsFreeSurferOasis )

tag <- c()
for( i in 1:nrow( resultsFreeSurferCombined ) )
  {
  tag <- append( tag, paste0( resultsFreeSurferCombined$ID[i], resultsFreeSurferCombined$SITE[i] ) );
  }
resultsFreeSurferCombined$TAG <- tag;

tagsCombined <- c( resultsANTsCombined$TAG, resultsFreeSurferCombined$TAG )

duplicateTags <- tagsCombined[duplicated( tagsCombined )]

for( h in 1:length( whichPipeline ) )
  {
  cat( whichPipeline[h], "\n" )

  resultsIXI <- read.csv( paste0( "labelresults", whichPipeline[h], "I.csv" ) )
  resultsKirby <- read.csv( paste0( "labelresults", whichPipeline[h], "K.csv" ) )
  resultsNKI <- read.csv( paste0( "labelresults", whichPipeline[h], "N.csv" ) )
  resultsOasis <- read.csv( paste0( "labelresults", whichPipeline[h], "O.csv" ) )

  resultsCombined <- rbind( resultsIXI, resultsKirby, resultsNKI, resultsOasis )

  resultsCombined$COLOR <- c( rep( 'red', nrow( resultsIXI ) ),
                              rep( 'green', nrow( resultsKirby ) ),
                              rep( 'blue', nrow( resultsNKI ) ),
                              rep( 'orange', nrow( resultsOasis ) ) )
  resultsCombined <- resultsCombined[order( resultsCombined$AGE ),]
  resultsCombined$SITE <- as.factor( resultsCombined$SITE )

  keepIndices <- c();
  for( i in 1:nrow( resultsCombined ) )
    {
    tag <- paste0( resultsCombined$ID[i], resultsCombined$SITE[i] )
    if( tag %in% duplicateTags )
      {
      keepIndices <- append( keepIndices, i );
      }
    }
  resultsCombined <- resultsCombined[keepIndices,];

  ################################################
  #
  # Plot by individuals and check for normality on the residuals
  #
  #################################################

  res <- residuals( lm( as.matrix( resultsCombined[,6:67] ) ~ resultsCombined$VOLUME + resultsCombined$SITE + resultsCombined$AGE ) )
  res.qvalues <- rep( NA, nrow( res ) )

  for( i in 1:ncol( res ) )
    {
    res[,i] <- ( res[,i] - min( res[,i] ) ) / ( max( res[,i] - min( res[,i] ) ) )
    }
  for( i in 1:length( res.qvalues ) )
    {
    res.qvalues[i] <- sf.test( res[i,] )$p.value
    }
  res.qvalues <- p.adjust( res.qvalues, "fdr" )

  ageIDLabels <- c();
  for( i in 1:nrow( resultsCombined ) )
    {
    tmp <- strsplit( as.character( resultsCombined$ID[i] ), "_" )
    id <- tmp[[1]][1]
    if( length( tmp[[1]] ) > 1 )
      {
      id <- tmp[[1]][2]
      }
    if( res.qvalues[i] <= 0.05 )
      {
      ageIDLabels <- append( ageIDLabels, paste0( id, "* (", resultsCombined$AGE[i], ")" ) )
#       cat( paste0( id, "* (", resultsCombined$AGE[i], ")\n" ) )
      } else {
      ageIDLabels <- append( ageIDLabels, paste0( id, " (", resultsCombined$AGE[i], ")\n" ) )
      }
    }

  pdf( paste0( "~/Desktop/BrainConstellationMap", whichPipeline[h], ".pdf" ) )
  stars( res,
         labels = ageIDLabels, cex = 0.2, scale = TRUE, radius = FALSE, full = TRUE, flip.labels = FALSE,
         mar = c( 0, 0, 2, 0 ),
         col.lines = resultsCombined$COLOR, main = paste0( "Brain Constellation Map (", whichPipeline[h], ")" ) )
  dev.off()
  }



