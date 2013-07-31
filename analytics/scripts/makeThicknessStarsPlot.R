# load required libraries
library( boot )
library( ggplot2 )
library( scales )
library( igraph )
library( reshape )
library( ANTsR )
library( ape )

resultsIXI <- read.csv( '../labelresultsI.csv' )
resultsKirby <- read.csv( '../labelresultsK.csv' )
resultsNKI <- read.csv( '../labelresultsN.csv' )
resultsOasis <- read.csv( '../labelresultsO.csv' )

resultsCombined <- rbind( resultsIXI, resultsKirby, resultsNKI, resultsOasis );
resultsCombined$SITE <- as.factor( resultsCombined$SITE )

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

thicknessStars <- resultsCombined[order(resultsCombined$AGE),6:37]
thicknessStarsLeft <- resultsCombined[order(resultsCombined$AGE),seq( 6, 37, by = 2 )]
thicknessStarsRight <- resultsCombined[order(resultsCombined$AGE),seq( 7, 37, by = 2 )]

thicknessStarsLeft <- ( thicknessStarsLeft - min( thicknessStars ) ) / ( max( thicknessStars ) - min( thicknessStars ) )
thicknessStarsRight <- ( thicknessStarsRight - min( thicknessStars ) ) / ( max( thicknessStars ) - min( thicknessStars ) )
thicknessStars <- ( thicknessStars - min( thicknessStars ) ) / ( max( thicknessStars ) - min( thicknessStars ) )

pdf( '../figs/thicknessStars.pdf' )
stars( t( thicknessStars ), scale = FALSE, full = TRUE, labels = corticalLabels,
       draw.segments = TRUE, main = 'Thickness Across Labels' )
dev.off()

pdf( '../figs/thicknessStarsLeft.pdf' )
stars( t( thicknessStarsLeft ), scale = FALSE, full = TRUE, labels = corticalLabels[seq( 1, 32, by = 2 )],
       draw.segments = TRUE, main = 'Thickness Across Labels' )
dev.off()

pdf( '../figs/thicknessStarsRight.pdf' )
stars( t( thicknessStarsRight ), scale = FALSE, full = TRUE, labels = corticalLabels[seq( 2, 32, by = 2 )],
       draw.segments = TRUE, main = 'Thickness Across Labels' )
dev.off()
