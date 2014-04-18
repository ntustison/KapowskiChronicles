#! /usr/bin/perl -w

use File::Basename;
use File::Path;
use File::Spec;

my @roiLabels = ( 1002, 1003, 1005, 1006, 1007, 1008, 1009, 1010, 1011,
                  1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020,
                  1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029,
                  1030, 1031, 1034, 1035,
                  2002, 2003, 2005, 2006, 2007, 2008, 2009, 2010, 2011,
                  2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020,
                  2021, 2022, 2023, 2024, 2025, 2026, 2027, 2028, 2029,
                  2030, 2031, 2034, 2035
                 );
my @roiNames = ( "left caudal anterior cingulate",
                 "left caudal middle frontal",
                 "left cuneus",
                 "left entorhinal",
                 "left fusiform",
                 "left inferior parietal",
                 "left inferior temporal",
                 "left isthmus cingulate",
                 "left lateral occipital",
                 "left lateral orbitofrontal",
                 "left lingual",
                 "left medial orbitofrontal",
                 "left middle temporal",
                 "left parahippocampal",
                 "left paracentral",
                 "left pars opercularis",
                 "left pars orbitalis",
                 "left pars triangularis",
                 "left pericalcarine",
                 "left postcentral",
                 "left posterior cingulate",
                 "left precentral",
                 "left precuneus",
                 "left rostral anterior cingulate",
                 "left rostral middle frontal",
                 "left superior frontal",
                 "left superior parietal",
                 "left superior temporal",
                 "left supramarginal",
                 "left transverse temporal",
                 "left insula",
                 "right caudal anterior cingulate",
                 "right caudal middle frontal",
                 "right cuneus",
                 "right entorhinal",
                 "right fusiform",
                 "right inferior parietal",
                 "right inferior temporal",
                 "right isthmus cingulate",
                 "right lateral occipital",
                 "right lateral orbitofrontal",
                 "right lingual",
                 "right medial orbitofrontal",
                 "right middle temporal",
                 "right parahippocampal",
                 "right paracentral",
                 "right pars opercularis",
                 "right pars orbitalis",
                 "right pars triangularis",
                 "right pericalcarine",
                 "right postcentral",
                 "right posterior cingulate",
                 "right precentral",
                 "right precuneus",
                 "right rostral anterior cingulate",
                 "right rostral middle frontal",
                 "right superior frontal",
                 "right superior parietal",
                 "right superior temporal",
                 "right supramarginal",
                 "right transverse temporal",
                 "right insula"
                 );

my $basedir = "/Users/ntustison/Desktop/ImportanceVisuals/";
my $csvfile =  "${basedir}/rfImportance.csv";
my $dktLabels = "${basedir}/antsMalfLabels.nii.gz";

open( FILE, "<${csvfile}" );
my @contents = <FILE>;
close( FILE );

my $freeSurferBrain = "${basedir}/rfImportanceFreeSurfer.nii.gz";
`cp $dktLabels $freeSurferBrain`;

my $antsBrain = "${basedir}/rfImportanceANTs.nii.gz";
`cp $dktLabels $antsBrain`;

my $maxImportance = 0;

for( my $i = 0; $i < @roiNames; $i++ )
  {
  my $roiName = $roiNames[$i];
  my @tmp = split( ' ', $roiName );
  $roiName = join( '.', @tmp );

  my $roiLabel = $roiLabels[$i];


  print "$roiName\n";

  my $idx = -1;
  for( my $j = 3; $j < @contents; $j++ )
    {
    my @line = split( ',', $contents[$j] );
    if( $line[0] =~ m/^$roiName$/ )
      {
      $idx = $j;
      }
    }
  if( $idx < 0 )
    {
    next;
    }

  my @line = split( ',', $contents[$idx] );
  my $antsImportance = $line[1];
  my $freeSurferImportance = $line[2];

  if( $antsImportance > $maxImportance )
    {
    $maxImportance = $antsImportance;
    }
  if( $freeSurferImportance > $maxImportance )
    {
    $maxImportance = $antsImportance;
    }

  `UnaryOperateImage 3 $freeSurferBrain r 0 $freeSurferBrain $roiLabel $freeSurferImportance`;
  `UnaryOperateImage 3 $antsBrain r 0 $antsBrain $roiLabel $antsImportance`;
  }

my $colorMask = "${basedir}/colorMmask.nii.gz";
my $totalMask = "${basedir}/totalMmask.nii.gz";

`ThresholdImage 3 $dktLabels $totalMask 0 0 0 1`;
`ThresholdImage 3 $dktLabels $colorMask 1000 1000000 1 0`;
`ImageMath 3 $colorMask MD $colorMask 1`;

`ImageMath 3 $antsBrain m $colorMask $antsBrain`;
`ImageMath 3 $antsBrain GD $antsBrain 2`;
`ImageMath 3 $freeSurferBrain m $colorMask $freeSurferBrain`;
`ImageMath 3 $freeSurferBrain GD $freeSurferBrain 2`;

my $freeSurferBrainRGB = "${basedir}/rfImportanceFreeSurferRGB.nii.gz";
my $antsBrainRGB = "${basedir}/rfImportanceANTsRGB.nii.gz";

`ConvertScalarImageToRGB 3 $antsBrain $antsBrainRGB $colorMask hot none 0 $maxImportance 0 255`;
`ConvertScalarImageToRGB 3 $freeSurferBrain $freeSurferBrainRGB $colorMask hot none 0 $maxImportance 0 255`;

print "antsSurf -s $totalMask -f [$antsBrainRGB,$colorMask,1] -d 1 \n";
print "antsSurf -s $totalMask -f [$freeSurferBrainRGB,$colorMask,1] -d 1 \n";