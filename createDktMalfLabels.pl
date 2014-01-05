#/usr/bin/perl -w

use strict;

use Cwd 'realpath';
use Switch;
use File::Find;
use File::Basename;
use File::Path;
use File::Spec;
use FindBin qw($Bin);

my $usage = qq{
  Usage: createIXILabels.pl <outputDir>
 };

my ( $outputDirectory ) = @ARGV;

my $dktDirectory = '/home/njt4n/share/Data/Public/OASIS-TRT-20/';

#  ======== need to change per cohort =======================

my $thicknessDirectory = '/home/njt4n/share/Data/Public/IXI/ThicknessAnts/';
my $dktXfrmPrefix = "${dktDirectory}/SyNToIXITemplate/T_template2_BrainCerebellumx";
my $template = "${thicknessDirectory}/../TotalTemplate/T_template2_BrainCerebellum.nii.gz";

#  ==========================================================

my @dktInverseWarps = <${dktXfrmPrefix}*InverseWarp.nii.gz>;


my $ANTSPATH = "/home/njt4n/share/Pkg/ANTs/bin/bin/";
my $UTILSPATH = "/home/njt4n/share/Pkg/Utilities/bin/";

if( ! -d $outputDirectory )
  {
  mkpath( $outputDirectory, {verbose => 0, mode => 0755} ) or
    die "Can't create output directory $outputDirectory\n\t";
  }

my @suffixList = ( ".nii.gz" );

my @commandFiles = <${thicknessDirectory}/*command.sh>;

my $count = 0;
for( my $i = 0; $i < @commandFiles; $i++ )
  {
  print "$commandFiles[$i]\n";

  my ( $filename, $directories, $suffix ) = fileparse( $commandFiles[$i], ".sh" );

  my $prefix = $filename;
  $prefix =~ s/command//;

  my $mask = "${directories}/${prefix}BrainExtractionMask.nii.gz";
  my $n4 = "${directories}/${prefix}BrainSegmentation0N4.nii.gz";
  my $skullStripped = "${outputDirectory}/${prefix}SkullStripped.nii.gz";
  my $malfImage = "${outputDirectory}/${prefix}-malf.nii.gz";

  if( -e $malfImage )
    {
    next;
    }

  my $commandFile = "${outputDirectory}/${prefix}command.sh";

  open( FILE, ">${commandFile}" );
  print FILE "#!/bin/sh\n";
  print FILE "export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1\n";
  print FILE "\n";

  print FILE "${UTILSPATH}/BinaryOperateImages 3 $n4 x $mask $skullStripped\n";

  if( ! -e "${outputDirectory}/T_templatex${prefix}Warped.nii.gz" )
    {
    print FILE "${ANTSPATH}/antsRegistrationSyN.sh -d 3 -f $template -m $skullStripped -o ${outputDirectory}/T_templatex${prefix}\n";
    }

  # register the images

  my @labelImages = ();
  my @warpedImages = ();
  for( my $j = 1; $j <= @dktInverseWarps; $j++ )
    {
    my $fileId = "OASIS-TRT-20-${j}";

    my $dktVolume = "${dktDirectory}/Volumes/${fileId}_brain.nii.gz";
    my $dktLabels = "${dktDirectory}/Labels/${fileId}_DKT31_CMA_labels.nii.gz";

    my $outputPrefix = "${outputDirectory}/${prefix}x${fileId}";

    my @args = ( "${ANTSPATH}/antsApplyTransforms",
                   '-d', 3,
                   '-i', $dktLabels,
                   '-r', $n4,
                   '-o', "${outputPrefix}LabelsWarped.nii.gz",
                   '-n', 'NearestNeighbor',
                   '-t', "[${outputDirectory}/T_templatex${prefix}0GenericAffine.mat,1]",
                   '-t', "${outputDirectory}/T_templatex${prefix}1InverseWarp.nii.gz",
                   '-t', "${dktXfrmPrefix}${fileId}1Warp.nii.gz",
                   '-t', "${dktXfrmPrefix}${fileId}0GenericAffine.mat"
                 );
    print FILE "@{args}\n";

    my @args = ( "${ANTSPATH}/antsApplyTransforms",
                   '-d', 3,
                   '-i', $dktVolume,
                   '-r', $n4,
                   '-o', "${outputPrefix}Warped.nii.gz",
                   '-n', 'Linear',
                   '-t', "[${outputDirectory}/T_templatex${prefix}0GenericAffine.mat,1]",
                   '-t', "${outputDirectory}/T_templatex${prefix}1InverseWarp.nii.gz",
                   '-t', "${dktXfrmPrefix}${fileId}1Warp.nii.gz",
                   '-t', "${dktXfrmPrefix}${fileId}0GenericAffine.mat"
                 );
    print FILE "@{args}\n";

    push( @labelImages, "${outputPrefix}LabelsWarped.nii.gz" );
    push( @warpedImages, "${outputPrefix}Warped.nii.gz" );
    }

  print FILE "${ANTSPATH}/jointfusion 3 -g @warpedImages -l @labelImages -m Joint[0.1,2] $n4 $malfImage\n";

  print FILE "rm ${outputDirectory}/T_templatex${prefix}1Warp.nii.gz ${outputDirectory}/T_templatex${prefix}1InverseWarp.nii.gz ${outputDirectory}/T_templatex${prefix}0GenericAffine.mat";
  print FILE "rm @labelImages @warpedImages\n";

  print FILE "\n";
  close( FILE );

  print "** labels ${filename}\n";
  $count++;
  system( "qsub -N ${count}ixidkt -v ANTSPATH=$ANTSPATH -q standard -l nodes=1:ppn=1 -l mem=10gb -l walltime=50:00:00 $commandFile" );
  }
