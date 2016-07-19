=head1 Description

    This script is used to generate normalized kmer
    feature vectors for each CRM and neg seq given
    both of CRM and neg seq fasta files as input.
    Also, the output file format is consistent to libsvm
    format.

=head1 Usage

    perl kmer2lib.pl PosFastaFile NegFastaFile OutDir

=cut


use strict;
use warnings;
use File::Basename;
use FindBin qw($Bin);

die `pod2text $0` if (@ARGV != 3);

# die "Usage:\n\tperl $0 PosFastaFile NegFastaFile OutDir\n" unless @ARGV==3;

# script to generate kmer features
my $kmer = "$Bin/kmers.pl";

my $posFile = $ARGV[0];
my $negFile = $ARGV[1];
my $outdir = $ARGV[2];

my $posFileName = basename($posFile);

# generate 6mer on including reverse complement
my @posKmerOut = `perl $kmer --rc $posFile 6`;
my @negKmerOut = `perl $kmer --rc $negFile 6`;

# output normalized CRM kmer features to file
open OUT1,">$outdir/$posFileName.andNeg.rc.6mer";
for my $line (@posKmerOut)
{
    print OUT1 "1 ";
    my $count = 1;
    chomp($line);
    my @array = split /\s+/,$line;
    my $vectorLen = vecLen(@array[3..$#array]); # vector size
    $vectorLen = 1e-50 if $vectorLen == 0; # handle zero vector, same as to print out all zeros for the norm kmer vector
    for my $elem (@array[3..$#array])
    {
        my $normValue = $elem / $vectorLen; # normalize feature value
        if ($normValue > 0)
        {
            print OUT1 "$count\:";
            printf OUT1 "%.5f ",$normValue;
        }
        $count = $count + 1;
    }
    print OUT1 "\n";
}

# output normalized neg seq kmer features to file
for my $line (@negKmerOut)
{
    print OUT1 "-1 ";
    my $count = 1;
    chomp($line);
    my @array = split /\s+/,$line;
    my $vectorLen = vecLen(@array[3..$#array]);
    $vectorLen = 1e-50 if $vectorLen == 0;
    for my $elem (@array[3..$#array])
    {
        my $normValue = $elem / $vectorLen;
        if ($normValue > 0)
        {
            print OUT1 "$count\:";
            printf OUT1 "%0.5f ",$normValue;
        }
        $count = $count + 1;
    }
    print OUT1 "\n";
}

close OUT1;

##===== 
# get kmer feature vector size
##=====
sub vecLen {
    my @vector = @_;
    my $sum=0;
    foreach (@vector){
        $sum += $_**2;
    }
    $sum = sqrt($sum);
}
