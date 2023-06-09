#!/usr/bin/perl -W
# ancestry_calculator.pl <stdin.txt

# calculate average depth from DP

# input format (tab delimited)
# ind \t chr \t start \t stop \t size \t ancestry
# 2831	1	15059	1056083	1041024	1
# 2831	1	1058192	1934433	876241	2
# 2835	1	1934442	2421900	487458	3

# where ancestry codes are:
# 1 - hom RF
# 2 - het RF/YF
# 3 - hom YF

# the output is formatted as:
# ind \t PropRF \t PropYF

use strict;
use warnings;

my $totalSize = 0 ;
my $sizeRF = 0 ;
my $sizeYF = 0 ;
my $tally = 0 ;

my ( $ind, $chr, $start, $stop, $size, $ancestry, $prior ) ;

print STDOUT "Ind \t PropRF \t PropYF\n" ;

while ( my $line = <STDIN> ) {
	chomp $line;

	( $ind, $chr, $start, $stop, $size, $ancestry ) = split ( "\t", $line ) ;

	if ( not defined $prior ) {
		$prior = $ind ;
		$tally += $size ;
		if ( $ancestry == 1 ) {
			$sizeRF += (2 * $size) ;
		} elsif ( $ancestry == 2 ) {
			$sizeRF += $size ;
			$sizeYF += $size ;
		} elsif ( $ancestry == 3 ) {
			$sizeYF += (2 * $size) ;
		}
	} elsif ( defined $prior ) {
		if ( $ind eq $prior ) {
			$tally += $size ;
			if ( $ancestry == 1 ) {
				$sizeRF += (2 * $size) ;
			} elsif ( $ancestry == 2 ) {
				$sizeRF += $size ;
				$sizeYF += $size ;
			} elsif ( $ancestry == 3 ) {
				$sizeYF += (2 * $size) ;
			}
		} elsif ( $ind ne $prior ) {
			$totalSize = 2 * $tally ;
			print STDOUT "$prior\t", sprintf("%.4f", $sizeRF/$totalSize), "\t", sprintf("%.4f", $sizeYF/$totalSize), "\n";
			$prior = $ind ;
			$tally = $size ;
			$totalSize = 0 ;
			$sizeYF = 0 ;
			$sizeRF = 0 ;
			if ( $ancestry == 1 ) {
				$sizeRF += (2 * $size) ;
			} elsif ( $ancestry == 2 ) {
				$sizeRF += $size ;
				$sizeYF += $size ;
			} elsif ( $ancestry == 3 ) {
				$sizeYF += (2 * $size) ;
			} 
		}
	}
}

$totalSize = 2 * $tally ;
print STDOUT "$prior\t",sprintf("%.4f", $sizeRF/$totalSize), "\t", sprintf("%.4f", $sizeYF/$totalSize), "\n";
