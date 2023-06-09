#!perl

use warnings;
use File::Basename; 

if (@ARGV < 1) { print "\nUSE: perl program.pl infile > outfile\n\n"; exit; }

# Chr_position order    ID  Ancestry
# 1_122579053 7312 11535 NA
# 2_158780 7313 11535 HET

my $infile = shift(@ARGV);
$filename = basename($infile);

open IN, $infile;
my @lines_raw = <IN>; chomp @lines_raw; close IN;
my @lines_filtered = grep !/^NA/, @lines_raw;

my $first_line = shift(@lines_filtered);
$first_line =~ s/\s+/\t/g;
my @first_parts = split(/\t/, $first_line); chomp @first_parts;
my @first_identifiers = split(/_/, $first_parts[0]); chomp @first_identifiers;

my $current_chromosome = $first_identifiers[0];
my $previous_chromosome = $first_identifiers[0];
my $start = $first_identifiers[1];
my $stop = $first_identifiers[1];
my $ancestry = $first_parts[3];

my $snps = 1;

$filename = basename($infile);
($outname) = fileparse($filename,'\..*'); 

open($out, ">", $outname."_blocks.txt");

#print $out "\#chromosome\tblock_start\tblock_stop\tancestry_of_block\tblock_size\tnum_snps_in_block\n";

foreach my $line (@lines_filtered)
{
    chomp $line;
    $line =~ s/\s+/\t/g;
    my @parts = split(/\t/, $line); chomp @parts;
    my @identifiers = split(/_/, $parts[0]); chomp @identifiers;
    $current_chromosome = $identifiers[0];

    if ($current_chromosome eq $previous_chromosome)
    {
        if ($ancestry eq $parts[3])
        {
            $stop = $identifiers[1];
            $snps++;
        }

        else
        {
            my $block_size = $stop - $start;
            print $out "$previous_chromosome\t$start\t$stop\t$ancestry\t$block_size\t$snps\n";
            $previous_chromosome = $identifiers[0];
            $start = $identifiers[1];
            $stop = $identifiers[1];
            $ancestry = $parts[3];
            $snps = 1;
        }
    }

    else
    {
        my $block_size = $stop - $start;
        print $out "$previous_chromosome\t$start\t$stop\t$ancestry\t$block_size\t$snps\n";
        $previous_chromosome = $identifiers[0];
        $start = $identifiers[1];
        $stop = $identifiers[1];
        $ancestry = $parts[3];
        $snps = 1;
    }
}

my $block_size = $stop - $start;

print $out "$previous_chromosome\t$start\t$stop\t$ancestry\t$block_size\t$snps\n";

close $out;

exit;
