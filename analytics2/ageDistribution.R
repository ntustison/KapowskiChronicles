library( ggplot2 )

resultsIXI <- read.csv( 'labelresultsANTsI.csv' )
resultsKirby <- read.csv( 'labelresultsANTsK.csv' )
resultsNKI <- read.csv( 'labelresultsANTsN.csv' )
resultsOasis <- read.csv( 'labelresultsANTsO.csv' )

plotData <- data.frame( cbind( SubjectTimePoint = results$TimePoint, Gender = gender, Site = site ) )
plotData <- transform( plotData, Gender = factor( Gender ) )
plotData$Site <- factor( site, levels = c( "IXI-Guys", "IXI-HH", "IXI-IOP", "MMRR", "NKI", "OASIS" ) )

myPlot <- ggplot( data = plotData, aes( x = Age, fill = Site ) ) +
             geom_histogram( binwidth = 2 ) +
             geom_bar() +
             facet_wrap( ~ Site ) +
             scale_x_continuous( "Age (years)", breaks = seq( 0, 100, by = 10 ), labels = seq( 0, 100, by = 10 ), limits = c( 0, 100 ) ) +
             scale_y_continuous( "Count" ) +
#              scale_colour_manual( values = c( male_color, female_color ), breaks = c( 1, 2 ), labels = c( "Male", "Female" ) ) +
              scale_fill_manual( values = c("#F8766D", "#C09B00", "#39B600", "#00C1A3", "#00B0F6", "#C77CFF"), breaks = 1:6, label = c( "IXI-Guys", "IXI-HH", "IXI-IOP", "MMRR", "NKI", "OASIS" ) ) +
#                theme( legend.justification = c( 0, 0 ), legend.position = c( 0, 0 ) ) +
             ggtitle( paste( "Age distribution by site", sep = "" ) )
ggsave( filename = paste0( "ageDistribution.pdf" ), plot = myPlot, width = 8, height = 6, units = 'in' )

# =============================================================
# this brings out the specific age - thickness relationship
# and also computes the "peak age" which is interesting ....
# =============================================================

if( FALSE )
  {
  thicknessValues <-as.matrix( resultsCombined[,6:37] )
  meanthickness<-apply(thicknessValues,FUN=mean,MARGIN=1)
  demog<-data.frame( resultsCombined , meanthickness )
  mysub<-meanthickness>1. & meanthickness<4 & demog$AGE > 1
  subdemog<-subset(demog,mysub)
  myform1<-paste("meanthickness ~ 1 + SEX + SITE + VOLUME ")
  myform2<-paste(myform1," + I(AGE) + I(AGE^2) + I(AGE^4)  ")
  mdl1<-lm( as.formula(myform1),data=subdemog)
  mdl2<-lm( as.formula(myform2),data=subdemog)
  pv<-anova(mdl1,mdl2)$P[[2]]
  visreg(mdl2, xvar=c("AGE","SEX"),main=as.character(pv))

  for ( x in 6:37 ) {
   myform1<-paste(colnames( resultsCombined )[x]," ~ 1 + SEX + SITE + VOLUME ")
   myform2<-paste(myform1," + I(AGE) + I(AGE^2)  ")
   print(myform2)
   mdl1<-lm( as.formula(myform1),data=subdemog)
   mdl2<-lm( as.formula(myform2),data=subdemog)
   thresid<-residuals( lm( myform1 ,data=subdemog) )
   mdl3<-lm( thresid ~ 1 + I(AGE) + I(AGE^2) , data=subdemog)
   mydf<-data.frame(AGE=c(10:160/2) )
   pa<-predict( mdl3 , newdata=mydf )
   peakage<-mydf$AGE[ which.max( pa ) ]
   pv<-anova(mdl1,mdl2)$P[[2]]
   visreg(mdl2, xvar=c("AGE","VOLUME"),main=paste(peakage,as.character(pv)))
   }
 }
