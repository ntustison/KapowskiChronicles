library( ANTsR )

#####################################################################################
mysurf <- renderSurfaceFunction(list(template,brain), alphasurf=c(0.3,0.3), surfval=0.5, smoothsval=1.5,  alphafunc=1, mycol="cadetblue1")
young<-15
mid<-35
mid2<-55
eld<-75
testagemat<-rbind( c( young , mid ) , rev( c( mid, mid2 ) ))
testagemat<-rbind( testagemat , rev( c( mid, eld ) ))
# testagemat<-rbind( ( c( young, mid ) ), ( c( mid, mid2 ) ))
# testagemat<-rbind( testagemat , ( c( mid, eld ) ))
pvageMat<-matrix( rep( NA, nrow( testagemat ) * length(corticalLabels) ), nrow=length(corticalLabels) )
for ( ntests in 1:nrow( testagemat ) )
# for ( ntests in 1:1 )
  {
maximumNumberOfPermutations<-1000
initialNetworkDifference <- c();
pb <- txtProgressBar( min = 0, max = maximumNumberOfPermutations, style = 3 )
testages<-testagemat[ntests,]
permutationCount <- rep( 0, numberOfLabels )
for( permutation in 0:maximumNumberOfPermutations )
  {
  samplesAge <- resultsCombined$AGE
  if( permutation > 0 & permutation < maximumNumberOfPermutations )
    {
      samplesAge <- sample( resultsCombined$AGE )
    }
    # localtransitivity betweeness pagerank
  gdensity <- 0.25
  ageDifference1 <- abs( samplesAge - testages[1] )
  resultsSubsetBasedOnAgeDifference1 <- subset( resultsCombined, ageDifference1 <= 10 )
  thicknessValues1 <- resultsSubsetBasedOnAgeDifference1[,thicknessColumns]
  ageDifference2 <- abs( samplesAge - testages[2] )
  resultsSubsetBasedOnAgeDifference2 <- subset( resultsCombined, ageDifference2 <= 10 )
  thicknessValues2 <- resultsSubsetBasedOnAgeDifference2[,thicknessColumns]
  thicknessValues1 <- residuals( lm( as.matrix( thicknessValues1 ) ~ resultsSubsetBasedOnAgeDifference1$SITE + resultsSubsetBasedOnAgeDifference1$SEX ) )
  thicknessValues2 <- residuals( lm( as.matrix( thicknessValues2 ) ~ resultsSubsetBasedOnAgeDifference2$SITE + resultsSubsetBasedOnAgeDifference2$SEX ) )
  correlationThicknessMatrix1 <- cor( thicknessValues1 )
  correlationThicknessMatrix2 <- cor( thicknessValues2 )
  myg1 <- makeGraph( correlationThicknessMatrix1 , gdensity )$localtransitivity # $walktrapcomm$modularity
  myg2 <- makeGraph( correlationThicknessMatrix2 , gdensity )$localtransitivity # $walktrapcomm$modularity
  networkDifference <- ( myg1 - myg2 )
  if( permutation == 0 )
    {
      initialNetworkDifference <- networkDifference
      mysign<-as.numeric( initialNetworkDifference > 0 )
      mysign[ mysign == 0 ] <- -1
      ( initialNetworkDifference )
    } else {
      permutationCount <- permutationCount + as.numeric( (networkDifference) > initialNetworkDifference )
    }
    if( permutation == 0 & FALSE )
      {
      locations<-list( vertices=centroids$vertices )
      myg1 <- makeGraph( correlationThicknessMatrix1 , gdensity )
      par3d(userMatrix=id, windowRect=c(0,0,512,512), zoom=0.7)
      par3d(userMatrix=lateralLeft, windowRect=c(0,0,512,512), zoom=0.7)
      renderNetwork( myg1$adjacencyMatrix , locations )
      fn<-paste('figs/temp_community_',testages[1],'.png',sep='')
      rgl.snapshot(fn)
      rgl.pop()
      fn<-paste('figs/temp_community_',testages[1],'.pdf',sep='')
      pdf(fn)
      plot( myg1$walktrapcomm, myg$mygraph )
      dev.off()
      myg2 <- makeGraph( correlationThicknessMatrix2 , gdensity )
      renderNetwork( myg2$adjacencyMatrix , locations )
      fn<-paste('figs/temp_community_',testages[2],'.png',sep='')
      rgl.snapshot(fn)
      rgl.pop()
      fn<-paste('figs/temp_community_',testages[2],'.pdf',sep='')
      pdf(fn)
      plot( myg2$walktrapcomm, myg$mygraph )
      dev.off()
      }
  setTxtProgressBar(pb, permutation )
} # end permutation
pvals <-  permutationCount / maximumNumberOfPermutations
pvageMat[,ntests]<-pvals
cat( "Ages: ", testages, "\n", sep = ' ' );
for ( ff in 1:length(mysign) ) if (  pvals[ff ] < 0.05 ) print(paste(corticalLabels[ff],initialNetworkDifference[ff],pvals[ff]))#
}

qq<-pvageMat
qq[ qq == 0 ]<-1/maximumNumberOfPermutations
qq<-matrix( as.numeric(p.adjust( qq , method='BH') <= 0.05 ) , ncol=nrow(testagemat) )
rownames( qq ) <- corticalLabels
colnames( qq ) <- testagemat[,1]
pdf('qvalues_network_ages.pdf')
pheatmap( qq , cluster_rows=F, cluster_cols=F)
dev.off()
pheatmap( qq , cluster_rows=F, cluster_cols=F)
