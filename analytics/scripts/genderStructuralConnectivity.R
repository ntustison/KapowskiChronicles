# load required libraries
library( boot )
library( ggplot2 )
library( scales )
library( igraph )
library( reshape )
library( ANTsR )
library( ape )

#################################################################
##
##   function definitions
##
##     * calculateCorrelationMatrix
##
#################################################################

calculateCorrelationMatrix <- function( mat, weights, nuis )
  {
  if( missing( mat ) | missing( weights ) )
    {
    print( args( calculateCorrelationMatrix ) )
    return( 1 )
    }
  correlationMatrix <- matrix( rep( NA, ncol( mat ) * ncol( mat ) ), ncol = ncol( mat ) )
  for ( x in 1:ncol( mat ) )
    {
    for ( y in 1:ncol( mat ) )
      {
      correlationMatrix[x,y] <- corr( cbind( mat[,x], mat[,y] ), w = weights / max( weights ) )
      }
    }
  return( correlationMatrix )
  }

#################################################################
##
##   main routine
##
#################################################################

maximumNumberOfPermutations <- 1000
sigma <- 5
ages <- seq( 10, 90, by = 5 )

resultsIXI <- read.csv( '../labelresultsI.csv' )
resultsKirby <- read.csv( '../labelresultsK.csv' )
resultsNKI <- read.csv( '../labelresultsN.csv' )
resultsOasis <- read.csv( '../labelresultsO.csv' )

resultsCombined <- rbind( resultsIXI, resultsKirby, resultsNKI, resultsOasis );
resultsCombined$SITE <- as.factor( resultsCombined$SITE )

# resultsCombined <- resultsCombined[which( resultsCombined$AGE >= 10 & resultsCombined$AGE <= 80 ),];

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

numberOfAges <- length( ages )
numberOfLabels <- length( corticalLabels )
thicknessColumns <- 6:ncol( resultsCombined )

tstatisticMatrix <- matrix( rep( NA, numberOfLabels * numberOfAges ), ncol = numberOfAges )
pvalueMatrix <- matrix( rep( NA, numberOfLabels * numberOfAges ), ncol = numberOfAges )
networkDifferenceMatrix <- matrix( rep( NA, numberOfLabels * numberOfAges ), ncol = numberOfAges )
networkMales <- matrix( rep( NA, numberOfLabels * numberOfAges ), ncol = numberOfAges )
networkFemales <- matrix( rep( NA, numberOfLabels * numberOfAges ), ncol = numberOfAges )
networkAll <- matrix( rep( NA, numberOfLabels * numberOfAges ), ncol = numberOfAges )

corrs <- rep( NA, numberOfAges )
weightedAges <- rep( NA, numberOfAges )

rownames( networkAll ) <- corticalLabels
colnames( networkAll ) <- ages


############################################
#
# Section 3.5. Gender Structural Connectivity Across Age Using Cortical Thickness
# permutation testing on sex differences across age
#
############################################

#
#  begin rendering code
#

# rgl.bg( color = "white" )

template <- antsImageRead( '../glasshead_female.nii.gz', 3 )
brain <- antsImageRead( '../glassbrain_female.nii.gz', 3 )
leftright <- antsImageRead( '../leftright_female.nii.gz', 3 )
nirepLabels <- antsImageRead( '../nirep_female.nii.gz', 3 )
centroids <- LabelImageCentroids( nirepLabels, physical = TRUE )

# template <- maskImage( template, leftright, 2 )
# brain <- maskImage( brain, leftright, 2 )

# id <- rotationMatrix( 0, 0, 1, 0)
# lateralLeft <- rotationMatrix( pi/2, 0, -1, 0 ) %*% rotationMatrix( pi/2, -1, 0, 0 )
# frontal <- rotationMatrix( pi*3/2, 1,  0, 0 )  # %*% rotationMatrix( pi/2, -1, 0, 0 )
# lateralRight <- rotationMatrix( pi/2, 0, 1, 0 ) %*% rotationMatrix( pi/2, -1, 0, 0 )
# par3d( userMatrix = id, windowRect = c( 0, 0, 512, 512 ), zoom = 0.7 )
# mysurf <- renderSurfaceFunction( list( template,brain ), alphasurf = c( 0.3, 0.3 ),
#   surfval = 0.5, smoothsval = 1.5, alphafunc = 1, mycol = "cadetblue1" )
#
# par3d( userMatrix = lateralLeft, windowRect = c( 0, 0, 512, 512 ), zoom = 0.7 )
# par3d( userMatrix = id, windowRect = c( 0, 0, 512, 512 ), zoom = 0.7 )
# par3d( userMatrix = lateralRight, windowRect = c( 0, 0, 512, 512 ), zoom = 0.7 )

#
# end rendering code
#

count <- 1
for( age in ages )
  {
  ageDifference <- ( resultsCombined$AGE - age )
  resultsSubsetBasedOnAgeDifference <- subset( resultsCombined, abs( ageDifference ) < 5000 * sigma )
  thicknessValues <- resultsSubsetBasedOnAgeDifference[,thicknessColumns]

  ageDifference <- ( resultsSubsetBasedOnAgeDifference$AGE - age )
  cweights <- exp( -1.0 * ageDifference^2 / sigma^2 )
  cweights <- cweights / sum( cweights )

  weightedAges[count] <- sum( resultsSubsetBasedOnAgeDifference$AGE * cweights )

#   thicknessResiduals <- residuals( lm( as.matrix( thicknessValues ) ~ resultsSubsetBasedOnAgeDifference$SITE ) )
  thicknessResiduals <- residuals( lm( as.matrix( thicknessValues ) ~ SITE + VOLUME, data = resultsSubsetBasedOnAgeDifference ) )

  initialNetworkDifference <- c()

  pb <- txtProgressBar( min = 0, max = maximumNumberOfPermutations, style = 3 )
  permutationCount <- rep( 0, numberOfLabels )
  for( permutation in 0:( maximumNumberOfPermutations - 1 ) )
    {
    samplesSex <- resultsSubsetBasedOnAgeDifference$SEX
    if( permutation > 0 & permutation < maximumNumberOfPermutations )
      {
      samplesSex <- sample( resultsSubsetBasedOnAgeDifference$SEX )
      }

    gdensity <- 0.25
    temp <- calculateCorrelationMatrix( thicknessResiduals, cweights )
    networkAll[,count] <- makeGraph( temp, gdensity )$localtransitivity

    if( permutation == 0 )
      {
      subnet0 <- reduceNetwork( temp, N = 200 )
      corrs[count] <- mean( subnet0$network[subnet0$network > 0] )
      }

    males <- samplesSex == 1
    temp <- calculateCorrelationMatrix( thicknessResiduals[males,], cweights[males] )
    networkMales[,count] <- makeGraph( temp, gdensity )$localtransitivity

    females <- samplesSex == 2
    temp <- calculateCorrelationMatrix( thicknessResiduals[females,], cweights[females] )
    networkFemales[,count] <- makeGraph( temp, gdensity )$localtransitivity

    networkDifference <- ( ( networkFemales[,count] ) - ( networkMales[,count] ) )

    if( permutation == 0 & TRUE )
      {
      locations <- list( vertices = centroids$vertices )

      ########## first males ############

#       gender <- samplesSex == 1
#       temp <- calculateCorrelationMatrix( thicknessResiduals[gender,], cweights[gender] )
#       maleGraph <- makeGraph( temp, gdensity )
#
#       renderNetwork( maleGraph$adjacencyMatrix, locations )
#
#       filename <- paste0( '../figs/temp_male_community_' , age, '.pdf' )
#       pdf( filename )
#       plot( maleGraph$walktrapcomm, maleGraph$mygraph, layout = layout.fruchterman.reingold,
#             vertex.size = 15,
#             vertex.label = corticalLabels, vertex.label.cex = 0.75, vertex.label.font = 2 )
#       dev.off()
#
#       filename <- paste0( '../figs/temp_male_community_X_' , age, '.pdf' )
#       pdf( filename )
#       V( maleGraph$mygraph )$name <- corticalLabels
#       dendPlot( walktrap.community( maleGraph$mygraph ), mode = "phylo", type = "fan", edge.width = 2,
#         show.tip.label = TRUE, x.lim = c( -30, 30 ), no.margin = TRUE, font = 2, label.offset = 1
#          )
# #       axisPhylo()
#       dev.off()
#
#       filename <- paste0( '../figs/temp_male_network_F_', age, '.png' )
#       par3d( userMatrix = id, windowRect = c( 0, 0, 512, 512 ), zoom = 0.7 )
#       par3d( userMatrix = frontal, windowRect = c( 0, 0, 512, 512 ), zoom = 0.7 )
#       rgl.snapshot( filename )
#
#       filename <- paste0( '../figs/temp_male_network_L_', age, '.png' )
#       par3d( userMatrix = id, windowRect = c( 0, 0, 512, 512 ), zoom = 0.7 )
#       par3d( userMatrix = lateralLeft, windowRect = c( 0, 0, 512, 512 ), zoom = 0.7 )
#       rgl.snapshot( filename )
#       filename <- paste0( '../figs/temp_male_network_R_', age, '.png' )
#
#       par3d( userMatrix = id, windowRect = c( 0, 0, 512, 512 ), zoom = 0.7 )
#       par3d( userMatrix = lateralRight, windowRect = c( 0, 0, 512, 512 ), zoom = 0.7 )
#       rgl.snapshot( filename )
#       rgl.pop()

      ######## now females ###########

#       gender <- samplesSex == 2
#       temp <- calculateCorrelationMatrix( thicknessResiduals[gender,], cweights[gender] )
#       femaleGraph <- makeGraph( temp, gdensity )
#
#       renderNetwork( femaleGraph$adjacencyMatrix, locations )
#
#       filename <- paste0( '../figs/temp_female_community_' , age, '.pdf' )
#       pdf( filename )
#       plot( femaleGraph$walktrapcomm, femaleGraph$mygraph, layout = layout.fruchterman.reingold,
#             vertex.size = 15,
#             vertex.label = corticalLabels, vertex.label.cex = 0.75, vertex.label.font = 2 )
#       dev.off()
#
#       filename <- paste0( '../figs/temp_female_community_X_' , age, '.pdf' )
#       pdf( filename )
#       V( femaleGraph$mygraph )$name <- corticalLabels
#       dendPlot( walktrap.community( femaleGraph$mygraph ), mode = "phylo", type = "fan", edge.width = 2,
#         show.tip.label = TRUE, x.lim = c( -30, 30 ), no.margin = TRUE, font = 2, label.offset = 1
#          )
# #       axisPhylo()
#       dev.off()
#
#       filename <- paste0( '../figs/temp_female_network_F_', age, '.png' )
#       par3d( userMatrix = id, windowRect = c( 0, 0, 512, 512 ), zoom = 0.7 )
#       par3d( userMatrix = frontal, windowRect = c( 0, 0, 512, 512 ), zoom = 0.7 )
#       rgl.snapshot( filename )
#
#       filename <- paste0( '../figs/temp_female_network_L_', age, '.png' )
#       par3d( userMatrix = id, windowRect = c( 0, 0, 512, 512 ), zoom = 0.7 )
#       par3d( userMatrix = lateralLeft, windowRect = c( 0, 0, 512, 512 ), zoom = 0.7 )
#       rgl.snapshot( filename )
#       filename <- paste0( '../figs/temp_female_network_R_', age, '.png' )
#
#       par3d( userMatrix = id, windowRect = c( 0, 0, 512, 512 ), zoom = 0.7 )
#       par3d( userMatrix = lateralRight, windowRect = c( 0, 0, 512, 512 ), zoom = 0.7 )
#       rgl.snapshot( filename )
#       rgl.pop()
      }

    if( permutation == 0 )
      {
      initialNetworkDifference <- networkDifference
      mysign <- as.numeric( initialNetworkDifference > 0 )
      mysign[mysign == 0] <- -1
      } else {
      permutationCount <- permutationCount + as.numeric( ( networkDifference * mysign ) >= initialNetworkDifference )
      }
    setTxtProgressBar( pb, permutation )
    }

  pvalueMatrix[,count] <- permutationCount / maximumNumberOfPermutations
  networkDifferenceMatrix[,count] <- initialNetworkDifference
  cat( "Age: ", age, "\n", sep = '' )
  for( ff in 1:length( mysign ) )
    {
    if( pvalueMatrix[ff, count] < 0.05 )
      {
      print( paste( corticalLabels[ff], initialNetworkDifference[ff], pvalueMatrix[ff, count] ) )
      #  cat( "  p-value (permutation testing) =", permutationCount / maximumNumberOfPermutations, "\n", sep = ' ' )
      }
    }
  count <- count+1
  }

#
#  Create plots
#

# (weighted) age vs. average thickness network plot

corrsPlotData <- data.frame( weightedAges = weightedAges, correlationValues = corrs )
corrsPlot <- ggplot( corrsPlotData, aes( x = weightedAges, y = correlationValues ) ) +
             geom_point( colour = "darkred", size = 4, alpha = 0.5 ) +
             stat_smooth( colour = "navyblue", formula = y ~ 1 + x + I(x^2) + I(x^3) + I(x^4), method = "lm",
                          size = 1, n = 1000, level = 0.95, se = TRUE, fullrange = TRUE, fill = 'black', alpha = 0.25 ) +
             scale_x_continuous( "Age", breaks = seq( 10, 80, by = 10 ), labels = seq( 10, 80, by = 10 ), limits = c( 10, 80 ) ) +
             scale_y_continuous( "Correlation", breaks = seq( 0.7, 1, by = 0.05 ), labels = seq( 0.7, 1, by = 0.05 ), limits = c( 0.685, 1 ) ) +
             ggtitle( "Average thickness network correlation vs. age" )
ggsave( filename = paste( "averageThicknessNetworkWithAge.pdf", sep = "" ), plot = corrsPlot, width = 8, height = 6, units = 'in' )


qvalueMatrix <- matrix( p.adjust( pvalueMatrix, method = "bonferroni" ), nrow = nrow( pvalueMatrix ), ncol = ncol( pvalueMatrix ) )

qvalueData <- data.frame( qvalueMatrix )
colnames( qvalueData ) <- ages
qvalueData$CorticalLabels <- factor( corticalLabels, levels = rev( corticalLabels ) )

qvaluePlot <- ggplot( melt( qvalueData ), aes( x = variable, y = CorticalLabels, fill = value ) ) +
              geom_tile( colour = "darkred" ) +
              scale_fill_gradientn( name = "q-value", colours = heat.colors( 7 ) ) +
              scale_x_discrete( 'Age', labels = seq( from = ages[1], to = ages[length( ages )], by = 5 ), breaks = seq( from = ages[1], to = ages[length( ages )], by = 5 ) ) +
              scale_y_discrete( 'Cortical Labels' )
ggsave( filename = "qvalueHeatMap.pdf", plot = qvaluePlot, width = 10, height = 6, units = 'in' )

networkMaleData <- data.frame( networkMales )
colnames( networkMaleData ) <- ages
networkMaleData$CorticalLabels <- factor( corticalLabels, levels = rev( corticalLabels ) )

networkMalePlot <- ggplot( melt( networkMaleData ) ) +
               geom_tile( aes( x = variable, y = CorticalLabels, fill = value ), colour = "gray50", size = 0 ) +
               scale_fill_gradientn( name = "transitivity\nvalues", colours = heat.colors( 7 ) ) +
               scale_x_discrete( 'Age', labels = seq( from = ages[1], to = ages[length( ages )], by = 5 ), breaks = seq( from = ages[1], to = ages[length( ages )], by = 5 ) ) +
               scale_y_discrete( 'Cortical Labels' ) +
               ggtitle( "Male network" )
ggsave( filename = "maleNetwork.pdf", plot = networkMalePlot, width = 10, height = 6, units = 'in' )


networkFemaleData <- data.frame( networkFemales )
colnames( networkFemaleData ) <- ages
networkFemaleData$CorticalLabels <- factor( corticalLabels, levels = rev( corticalLabels ) )

networkFemalePlot <- ggplot( melt( networkFemaleData ) ) +
               geom_tile( aes( x = variable, y = CorticalLabels, fill = value ), colour = "gray50", size = 0 ) +
               scale_fill_gradientn( name = "transitivity\nvalues", colours = heat.colors( 7 ) ) +
               scale_x_discrete( 'Age', labels = seq( from = ages[1], to = ages[length( ages )], by = 5 ), breaks = seq( from = ages[1], to = ages[length( ages )], by = 5 ) ) +
               scale_y_discrete( 'Cortical Labels' ) +
               ggtitle( "Female network" )
ggsave( filename = "femaleNetwork.pdf", plot = networkFemalePlot, width = 10, height = 6, units = 'in' )

networkAllData <- data.frame( networkAll )
colnames( networkAllData ) <- ages
networkAllData$CorticalLabels <- factor( corticalLabels, levels = rev( corticalLabels ) )

networkAllPlot <- ggplot( melt( networkAllData ) ) +
               geom_tile( aes( x = variable, y = CorticalLabels, fill = value ), colour = "gray50", size = 0 ) +
               scale_fill_gradientn( name = "transitivity\nvalues", colours = heat.colors( 7 ) ) +
               scale_x_discrete( 'Age', labels = seq( from = ages[1], to = ages[length( ages )], by = 5 ), breaks = seq( from = ages[1], to = ages[length( ages )], by = 5 ) ) +
               scale_y_discrete( 'Cortical Labels' ) +
               ggtitle( "Both genders network" )
ggsave( filename = "allNetwork.pdf", plot = networkAllPlot, width = 10, height = 6, units = 'in' )

#
# pvals <- rep( NA, nrow( networkAll ) )
# for ( n in 1:32 )
#   {
#   dd<-summary( lm( networkFemales[n,] ~ I(ages) + I(ages)^2 ) )
#   dd<-summary( lm( networkMales[n,] ~  I(ages) + I(ages)^2 ) )
#   dd<-summary( lm( networkAll[n,] ~  I(ages) + I(ages^2) ) )
#   pvals[n]<-coefficients(dd)[3,4]
#   }
# print( "transitivity with age" )
# print( p.adjust( pvals, method = 'BH' ) )
#

#   |======================================================================| 100%Age: 10
# [1] "L superior frontal 0.410526315789474 0.014"
# [1] "R superior parietal 0.4 0.042"
#   |======================================================================| 100%Age: 15
# [1] "L superior frontal 0.284954811270601 0.046"
# [1] "R superior frontal 0.326450742240216 0.004"
# [1] "L inferior 0.330555555555555 0.027"
#   |======================================================================| 100%Age: 20
# [1] "R superior frontal 0.311038961038961 0.013"
#   |======================================================================| 100%Age: 25
#   |======================================================================| 100%Age: 30
# [1] "L superior temporal 1 0.022"
# [1] "L inferior 0.417582417582418 0.045"
# [1] "L orbital frontal 0.5 0.003"
# [1] "L superior parietal 0.6 0.003"
#   |======================================================================| 100%Age: 35
# [1] "R superior temporal 0.80952380952381 0.01"
# [1] "L inferior parietal 0.265079365079365 0.039"
#   |======================================================================| 100%Age: 40
# [1] "R infero temporal 0.214869281045752 0.044"
# [1] "R inferior parietal 0.678571428571429 0.001"
#   |======================================================================| 100%Age: 45
# [1] "R superior frontal 0.334848484848485 0.022"
# [1] "R inferior 0.423931623931624 0.021"
#   |======================================================================| 100%Age: 50
# [1] "R parahippocampal 0.566666666666667 0.026"
# [1] "R superior frontal 0.335930735930736 0.017"
#   |======================================================================| 100%Age: 55
# [1] "L occipital 0.555555555555556 0.026"
# [1] "R occipital 0.866666666666667 0.009"
#   |======================================================================| 100%Age: 60
# [1] "R occipital 0.722222222222222 0.001"
# [1] "R orbital frontal 0.435087719298246 0.016"
#   |======================================================================| 100%Age: 65
#   |======================================================================| 100%Age: 70
#   |======================================================================| 100%Age: 75
# [1] "R parahippocampal 0.833333333333333 0.021"
# [1] "L frontal pole 0.395238095238095 0.008"
# [1] "L inferior parietal 0.329365079365079 0.049"
#   |======================================================================| 100%Age: 80
# [1] "L middle frontal 0.39047619047619 0.017"
# [1] "R inferior 0.533333333333333 0.001"
# [1] "R inferior parietal 0.335064935064935 0.003"
# [1] "R postcentral 0.297727272727273 0.018"
#   |======================================================================| 100%Age: 85
# [1] "R inferior parietal 0.256818181818182 0.017"
#   |======================================================================| 100%Age: 90




