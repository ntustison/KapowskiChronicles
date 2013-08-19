library( kernlab )
library( randomForest )
library( ggplot2 )

trainingData <- read.csv( "../trainingBrainSegmentationPosteriors2Projections.csv" )
testingData <- read.csv( "../testingBrainSegmentationPosteriors2Projections.csv" )
# trainingData <- read.csv( "../trainingCorticalThicknessProjections.csv" )
# testingData <- read.csv( "../testingCorticalThicknessProjections.csv" )

# Remove ID, gender, and site

trainingData <- trainingData[, 3:ncol( testingData )]
testingData <- testingData[, 3:ncol( testingData )]

brainAgeVM <- rvm( AGE ~ ., data = trainingData, type = "regression",
                    kernel = "laplacedot",
                    verbosity = 1, tol = .Machine$double.eps,
                    minmaxdiff = 1e-3, cross = 0, fit = TRUE,
                    na.action = na.omit , iterations = 1000 )
predictedAge <- predict( brainAgeVM, testingData )

# brainAgeRF <- randomForest( AGE ~ ., data = trainingData,
#                     na.action = na.omit, replace = FALSE, ntree = 2000 )
# predictedAge <- predict( brainAgeRF, testingData )

correlation <- cor( testingData$AGE, predictedAge )
meanerr <- mean( abs(  testingData$AGE - predictedAge ) )

cat( "Correlation (true vs. predicted age): ", correlation, ", mean error ", meanerr , "\n", sep = '' )

plotData <- data.frame( TrueAge = testingData$AGE, PredictedAge = predictedAge )

brainAgeRegression <- lm( testingData$AGE ~ 1 + predictedAge )

annotation <- paste0( "r = ", signif( correlation, 2 ), ", mean error = ", signif( meanerr, 2 ), " years" )

brainAgePlot <- ggplot( plotData, aes( x = TrueAge, y = PredictedAge ) ) +
                stat_smooth( colour = "navyblue", formula = y ~ 1 + x, method = "lm",
                  size = 1, n = 1000, level = 0.95, se = TRUE, fullrange = TRUE, fill = 'black', alpha = 0.5 ) +
                geom_point( aes( colour = "darkred" ), size = 4, alpha = 0.5 ) +
                scale_x_continuous( "True age (years)", breaks = seq( 0, 100, by = 10 ), labels = seq( 0, 100, by = 10 ), limits = c( 0, 100 ) ) +
                scale_y_continuous( "Predicted age (years)", breaks = seq( 0, 100, by = 10 ), labels = seq( 0, 100, by = 10 ), limits = c( 0, 100 ) ) +
                scale_colour_manual( labels = annotation, values = "darkred" ) +
                theme( legend.justification = c( 1, 0 ), legend.position = c( 1, 0 ), legend.title = element_blank(), legend.background = element_blank(), legend.key = element_blank(), legend.text = element_text( colour = 'navyblue' ) ) +
                guides( colour = guide_legend( override.aes = list( shape = NA ) ), keywidth = 0 )
ggsave( filename = paste( "../brainAgeBrainSegmentationPosteriors2.pdf", sep = "" ), plot = brainAgePlot, width = 6, height = 6, units = 'in' )



# plot the results (Bland Altman)
#
# ageDifference <- ( testingData$AGE - predictedAge )
# meanAgeDifference <- mean( ageDifference )
# stdAgeDifference <- sd( ageDifference )
#
# upperLineIntercept <- meanAgeDifference + 0.95 * stdAgeDifference
# lowerLineIntercept <- meanAgeDifference - 0.95 * stdAgeDifference
#
# plotData <- data.frame( TrueAge = testingData$AGE, AgeDifference = ageDifference )
#
# brainAgePlot <- ggplot( plotData, aes( x = TrueAge, y = AgeDifference ) ) +
#                 geom_point( colour = "darkred", size = 4, alpha = 0.75 ) +
#                 geom_hline( aes( yintercept = upperLineIntercept ), colour = "navyblue", linetype = "dashed" ) +
#                 geom_hline( aes( yintercept = lowerLineIntercept ), colour = "navyblue", linetype = "dashed" ) +
#                 scale_x_continuous( "True age", breaks = seq( 10, 80, by = 10 ), labels = seq( 10, 80, by = 10 ), limits = c( 10, 80 ) ) +
#                 scale_y_continuous( "True - predicted age", breaks = seq( -20, 20, by = 5 ), labels = seq( -20, 20, by = 5 ), limits = c( -20, 20 ) ) +
#                 ggtitle( "True vs. Predicted Age" )
# ggsave( filename = paste( "../brainAge.pdf", sep = "" ), plot = brainAgePlot, width = 8, height = 6, units = 'in' )
#

if( FALSE )
  {
  img<-as.matrix( trainingData[ ,c(1,3:ncol(trainingData) ) ] )
  testimg<-as.matrix( testingData[ ,c(1,3:ncol(trainingData) ) ] )
  myregression <- sparseRegression( img , data.frame( AGE=trainingData$AGE , SEX=trainingData$SEX ), "AGE",  sparseness=-0.1, nvecs=10, its=10, cthresh=0)
  ff<-data.frame(AGE=trainingData$AGE, myregression$umatrix )
  tt<-glm( AGE ~ ., data = ff )
  cor.test( predict( tt ) , ff$AGE )

  newu<-  testimg %*% ( as.matrix( myregression$eigenanatomyimages ) )
  colnames( newu ) <- colnames( myregression$umatrix  )
  names( newu ) <- names( myregression$umatrix  )
  ff<-data.frame(  newu )
  predictedAge2 <- predict( tt, newdata=ff )
  correlation <- cor( testingData$AGE, predictedAge2 )
  meanerr<-(  mean( abs(  testingData$AGE - predictedAge2 ) ) )
  cat( "Correlation (true vs. predicted age): ", correlation, " mean err " , meanerr , "\n", sep = '' )
  }
