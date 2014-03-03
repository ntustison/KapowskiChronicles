#! /usr/bin/perl -w

use strict;
use Cwd 'realpath';

use File::Find;
use File::Basename;
use File::Path;
use File::Spec;
use FindBin qw($Bin);

my $usage = qq{
Usage: processLogFiles.pl <inputDir> <outputFile> <isFreesurfer>
 };

my ( $inputDirectory, $outputFile, $isFreesurfer ) = @ARGV;

open( FILE, ">${outputFile}" );
print FILE "Id,ElapsedTime\n";

my @logFiles = <${inputDirectory}/*.o*>;

for( my $i = 0; $i < @logFiles; $i++ )
  {
  open( FILE2, "<${logFiles[$i]}" );
  my @contents = <FILE2>;
  close( FILE2 );

  my $elapsedTime = 0;

  if( $isFreesurfer == 0 )
    {
    my $timeString = ${contents[-3]};
    chomp( $timeString );
    my @tokens = split( ' ', $timeString );
    if( @tokens < 3 )
      {
      print "${logFiles[$i]}\n";
      next;
      }
    $elapsedTime = ( $tokens[3] / 3600 );

    my $index = -1;
    for( my $j = 0; $j < @contents; $j++ )
      {
      if( $contents[$j] =~ m/N4Truncated0/ )
        {
        $index = $j;
        last;
        }
      }
    if( $index == -1 )
      {
      print "${logFiles[$i]}\n";
      next;
      }

    @tokens = split( ' ', $contents[$index] );
    for( my $j = 0; $j < @tokens; $j++ )
      {
      if( $tokens[$j] =~ m/N4Truncated0/ )
        {
        $index = $j;
        last;
        }
      }
    ( my $id = $tokens[$index] ) =~ s/N4Truncated0\.nii\.gz//;
    my @comps = split( '/', $id );
    $id = ${comps[-1]};
    print FILE "$id,$elapsedTime\n";
    }
  else
    {
    my $idString = ${contents[-2]};
    chomp( $idString );
    my @tokens = split( ' ', $idString );
    my $id = $tokens[2];

    my $timeString = ${contents[-3]};
    chomp( $timeString );
    @tokens = split( ' ', $timeString );
    $elapsedTime = $tokens[2];
    print FILE "$id,$elapsedTime\n";
    }
 }

close( FILE );

