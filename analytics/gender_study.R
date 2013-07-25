library(boot)
library(ANTsR)
library(igraph)
n<-256 # for graph
mymaxperm<-5
sigma<-3
ageth<-5000*sigma # sigma*sigma
ages<-seq(10,80,by=2)# as.numeric(10:80)
sigval<-0.01
for ( controlvol in c( TRUE ) )
{

makegraph <- function( myrsfnetworkcorrs  ) {
    ewt<-c(0)
    neig<-nrow( myrsfnetworkcorrs )
    if ( neig == 0 ) return( 0 )
    corthresh<-0.000001
    for ( x in c( 1:neig ) )
      for ( y in c( ( x : neig ) ) )
        if ( myrsfnetworkcorrs[ x, y] > corthresh & myrsfnetworkcorrs[ x, y] < 1 )
          ewt<-c(ewt,myrsfnetworkcorrs[x,y])
    ewt<-ewt[2:length(ewt)]
    adjmat<-as.matrix( ( myrsfnetworkcorrs > corthresh & myrsfnetworkcorrs < 1 ) , nrow=neig,ncol=neig)
    g1<-graph.adjacency( adjmat,mode=c("undirected"))
#    gmetric<-betweenness(g1,normalized=T,weights=1/ewt) # betweenness( g1 )
#    gmetric<-closeness(g1,normalized=T,weights=1/ewt) # betweenness( g1 )
#    gmetric<-page.rank( g1 ,weights=1/ewt)$vector
    gmetric<-transitivity( g1 ,type="local",isolates=c("zero"))
#    gmetric<-degree( g1 )
    return (  gmetric )
    }
corw2 <- function( mat, weights , nuis )
  {
  if ( missing( mat ) | missing( weights ) )
  {
  print( args(corw) )
  return( 1 )
  }
  cormat<-matrix(rep(NA,ncol(mat)*ncol(mat)),ncol=ncol(mat))
  for ( x in 1:ncol(mat) )
    {
  for ( y in 1:ncol(mat) )
    {
      cormat[x,y]<-corr( cbind(mat[,x], mat[,y]), w = weights/max(weights) )
#    cormat[x,y]<-sqrt( summary( lm( mat[,x] ~ mat[,y] + nuis ), weights = weights/sum(weights)  )$r.squared )
    }
    }
  return( cormat )
  }
dati<-read.csv('labelresultsI.csv')
datk<-read.csv('labelresultsK.csv')
datn<-read.csv('labelresultsN.csv')
dato<-read.csv('labelresultsO.csv')
allth<-rbind( dati, datk )
allth<-rbind( allth, datn )
allth<-rbind( allth, dato )
allth$SITE<-as.factor( allth$SITE )
#
# for permutation!
# print("PERMUTING!") ; allth$AGE<-allth$AGE[ sample( 1:nrow(allth) ) ]
#
thcols<-6:ncol(allth)
wt_ages<-rep(NA,length(ages))
corrs<-rep(NA,length(ages))
ct<-1
ncols<-length(thcols)
pgenmat<-matrix(rep(NA,ncols*length(ages)),ncol=length(ages))
genmat<-matrix(rep(NA,ncols*length(ages)),ncol=length(ages))
thkmat<-matrix(rep(NA,ncols*length(ages)),ncol=length(ages))
ntw1<-thkmat
ntw2<-thkmat
for ( age in ages )
  {
  dage<- (allth$AGE - age)
  mycond<- abs( dage ) < ageth
  ddage<-subset(allth, mycond )
  thk<-ddage[,thcols]
  dage<- (ddage$AGE - age)
  cweights<-exp( -1.0 * dage * dage / sigma^2 )
  cweights<-cweights/sum(cweights)
  wage<-sum( ddage$AGE * cweights )
  wt_ages[ct]<-wage
  corthk<-corw2( thk , weights = cweights  )
  ############################################
  pct<-rep(0,ncol(thk))
  maxperm<-mymaxperm
  for ( perm in 0:maxperm ) {
  myth<-residuals( lm( as.matrix(thk) ~ ddage$SITE ) )
  testvals<-( ddage$SEX )
  if ( perm > 0 & perm < maxperm ) testvals<-sample( ddage$SEX )
  w1<-testvals == 1
  w2<-testvals == 2

  temp<-corw2( myth[ w1 , ] ,  cweights[ w1 ] )
  subnet0 <- reduceNetwork( temp, N=n )
  g0<-makegraph( subnet0$network ) # mean(subnet0$network[ subnet0$network >0 ])
  ntw1[,ct]<-g0

  temp<-corw2( myth[ w2, ] ,  cweights[  w2 ] )
  subnet0 <- reduceNetwork( temp, N=n )
  g0<-makegraph( subnet0$network )
  ntw2[,ct]<-g0

  dif<-abs( ntw1[,ct] -  ntw2[,ct] )
  if ( perm == 0 | perm == maxperm ) odif<-dif else  pct<-pct+as.numeric( dif >= odif )
  }
  print( age )
  print( pct / maxperm )
  ############################################

  for ( corlab in 1:32 ) {
    #
    myth<-thk[ ,corlab]
    myth<-residuals( lm( thk[ ,corlab] ~ ddage$SITE ) )
    meandiff<-sum( myth[ ddage$SEX == 2 ] * cweights[  ddage$SEX == 2 ] )
    meandiff<-meandiff - sum( myth[ ddage$SEX == 1 ] * cweights[  ddage$SEX == 1 ] )
    thkmat[ corlab, ct ]<-meandiff
    myage<-as.numeric( ddage$AGE )
#    if ( controlvol ) momo<-summary( lm( thk[,corlab] ~ SEX +  I(myage) + I(myage^2) + SITE + VOLUME , weights = cweights , data=ddage ) ) else momo<-summary( lm( thk[,corlab] ~ SEX +  I(myage) + I(myage^2) + SITE, weights = cweights , data=ddage ) )
    if ( controlvol ) momo<-summary( lm( thk[,corlab] ~ SEX + SITE + VOLUME , weights = cweights , data=ddage ) ) else momo<-summary( lm( thk[,corlab] ~ SEX + SITE, weights = cweights , data=ddage ) )
    val<-coef(  momo )[2,4]
    pgenmat[corlab,ct]<-val
#    if ( val < 0.000000001 ) print( (momo ) )
    val<-coef(  momo )[2,3]
    genmat[corlab,ct]<-val
  }
#  antsImageWrite(as.antsImage(corthk),paste('age',age,'network.nii.gz',sep=''))
#  image( corthk )
  corrs[ct]<-mean( cor( thk ) )
  n<-200
  subnet0 <- reduceNetwork( corthk, N=n )
  corrs[ct]<-mean(subnet0$network[ subnet0$network >0 ])
#  corrs[ct]<-sum( cor( thk )[ cor(thk) > 0.85 ]  )
  meanthk <- mean( apply( thk, FUN=mean, MARGIN=2 ) )
#  print(paste('nsub:',nrow(thk),'age:',age,'w-age:',wage,'mean-corr:',corrs[ct], 'mean-thickness:',meanthk  ) )
  ct<-ct+1
  }
mdl<-lm(  corrs ~  wt_ages + I((wt_ages)^2)  + I((wt_ages)^3)   + I((wt_ages)^4)   )
# print(summary(mdl))
plot(wt_ages,corrs,main='avg thickness network correlation against age')
points( wt_ages, predict(mdl), col='red', type='l')
pdf('thickness_corr_with_age.pdf',width=10,height=7)
plot(wt_ages,corrs,main='avg thickness network correlation against age')
points( wt_ages, predict(mdl), col='red', type='l')
dev.off()


library(randomForest)
library(e1071)
myform<-paste("AGE~SEX+SITE+VOLUME+",paste("I(",colnames(allth)[thcols],collapse='^2+',")"),"^2+",paste("I(",colnames(allth)[thcols],collapse='+',")"),"",sep='')
mydat<-data.frame( allth )
th<-0.08
mycond<-rnorm( nrow( allth ) )
train<-subset( mydat, mycond >= th  )
test<-subset( mydat, !mycond >= th  ); dim(test); dim(train)
mdl<-glm( myform , data = train )
predcols<-c(3,5:ncol(train))
mdl<-svm( y=train$AGE , x=train[,predcols] )
predage <- predict( mdl, test[,predcols]  )
trueage<-allth$AGE[ !mycond >= th ]
# plot( predage, trueage )
print( mean( abs( predage - trueage ) ) )

cortical_labels <- c( "L occipital", "R occipital",
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
                     "L postcentral",       "R postcentral" );

library(pheatmap)
qmat<-p.adjust( pgenmat, method="bonf" )
qmat<-matrix(as.numeric( qmat < sigval ) , nrow=nrow(genmat) )
qmat[ ( qmat == 1 ) & ( genmat < 0 ) ]<-( -1 )
qgenmat<-data.frame( qmat )
row.names( qgenmat ) <- cortical_labels;
# row.names(qgenmat)<-names(thk)
colnames(qgenmat)<-ages
genmat<-data.frame(genmat)
row.names( genmat ) <- cortical_labels;
# row.names(genmat)<-names(thk)
colnames(genmat)<-ages
pdf(paste('sex_v_age',controlvol,'.pdf',sep=''),width=10,height=7)
rr<-max(abs(genmat))
mybreaks<-c(0:100)/100*rr*2-rr
mybreaks[1]<-( -1.0 * rr )
mybreaks[101]<-(  rr )
pheatmap( genmat ,breaks=mybreaks, cluster_rows = F , cluster_cols = F , display_numbers = F)
dev.off()
pheatmap( qgenmat , cluster_rows = F , cluster_cols = F , show_rownames = T, show_colnames = T )
pdf(paste('sex_v_age',controlvol,'_p.pdf',sep=''),width=10,height=7)
pheatmap( qgenmat , cluster_rows = F , cluster_cols = F  )
dev.off()
}


rr<-max(abs(rbind(ntw1,ntw2)))
mybreaks<-c(0:100)/100*rr*2-rr
mybreaks[1]<-( -1.0 * rr )
mybreaks[101]<-(  rr )
ct<-1
for ( ntdif in list( ntw1 , ntw2 ) )
{
# row.names(ntdif)<-names(thk)
row.names( ntdif ) <- cortical_labels;
colnames(ntdif)<-ages
pdf(paste('network_sex_v_age_',ct,'.pdf',sep=''),width=10,height=7)
pheatmap( ntdif ,breaks=mybreaks, cluster_rows = F , cluster_cols = F , display_numbers = F)
dev.off()
ct<-ct+1
}
