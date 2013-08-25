library( ggplot2 )
library( ANTsR )
library( visreg )

############################################
#
# Section 3.4:  Gender and Age Relationships with DiReCT Cortical Thickness
#   Look at basic regression results, i.e.
#     thickness ~ 1 + AGE + AGE^2 + GENDER + SITE + VOLUME
#
############################################

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

male_color = "darkred";
female_color = "navyblue";

gender <- cut( resultsCombined$SEX, breaks = c( 0.5, 1.5, 2.5 ), label = c( "male", "female" ) );
gender[which( gender == 1 )] <- 'male';
gender[which( gender == 2 )] <- 'female';

for( i in 1:32 )
  {
  myformula <- as.formula( paste0( "LABEL_", i, " ~ 1 + AGE + I(AGE^2) + VOLUME + SITE + SEX" ) )
  cat( paste0( corticalLabels[i], " ~ 1 + SITE + SEX + AGE + I(AGE^2) + VOLUME" ), "\n" )
  mylm <- lm( myformula, data = resultsCombined )
  mysummary <- summary( mylm )
  print( mysummary$coefficients )
  cat( "------------------------------------------\n" )

#   thickness <- mylm$residuals
  thickness <- resultsCombined[,i+5]

  plotData <- data.frame( cbind( Age = resultsCombined$AGE, Thickness = thickness, Gender = gender ) )
  plotData <- transform( plotData, Gender = factor( Gender ) );

  thickPlot <- ggplot( plotData, aes( x = Age, y = Thickness, group = Gender ) ) +
               stat_smooth( aes( group = Gender, colour = Gender ), formula = y ~ 1 + x + I(x^2) , method = "lm", size = 1, n = 1000, level = 0.95, se = TRUE, fullrange = TRUE, fill = 'black', alpha = 0.5 ) +
               geom_point( data = plotData, aes( colour = Gender, shape = Gender ), size = 3, alpha = 0.5 ) +
               scale_x_continuous( "Age (years)", breaks = seq( 10, 90, by = 10 ), labels = seq( 10, 90, by = 10 ), limits = c( 10, 90 ) ) +
#                scale_y_continuous( "Thickness (mm)", breaks = seq( -1, 1, by = 1 ), labels = seq( -1, 1, by = 1 ), limits = c( -1, 1 ) ) +
               scale_y_continuous( "Thickness (mm)", breaks = seq( 0, 6, by = 1 ), labels = seq( 0, 6, by = 1 ), limits = c( 0, 6 ) ) +
               scale_colour_manual( values = c( male_color, female_color ), breaks = c( 1, 2 ), labels = c( "Male", "Female" ) ) +
               scale_shape_manual( values = c( 18, 16 ), breaks = c( 1, 2 ), labels = c( "Male", "Female" ) ) +
               theme( legend.justification = c( 0, 0 ), legend.position = c( 0, 0 ) ) +
               ggtitle( paste( "Cortical thickness (", corticalLabels[i], ")", sep = "" ) )
  ggsave( filename = paste0( "../figs/label", i, "_results.pdf" ), plot = thickPlot, width = 8, height = 6, units = 'in' )
  }

# =============================================================
# this brings out the specific age - thickness relationship
# and also computes the "peak age" which is interesting ....
# =============================================================

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












