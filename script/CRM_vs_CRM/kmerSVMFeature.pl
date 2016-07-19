=head1 Description

    This script is used to generate normalized kmer 
    features in libsvm format for both of test and 
    training data.

=head1 Usage

    perl kmerSVMFeature.pl CRMname Outdir

=cut

use strict;
use warnings;
use File::Basename;
use FindBin qw($Bin);

# die "Usage:\n\tperl $0 CRMname Outdir\n" unless @ARGV==2;
die `pod2text $0` if (@ARGV!=2);

# create dir at slave nodes
# my $tmpDir = $ENV{TMPDIR};
# my $tmpModelDir = "$tmpDir/output";
# `mkdir $tmpModelDir` unless (-e "$tmpModelDir");

# script to generate normalized kmer feature in libsvm format given CRM and neg fasta files
my $kmer2lib = "$Bin/kmer2lib.pl";

my $crmName = $ARGV[0]; # CRMname
my $outdir = $ARGV[1];

# 10 trials
for (my $k=1;$k<=10;$k++)
{
    # `mkdir $tmpModelDir/time$k` unless (-e "$tmpModelDir/time$k");
    warn "trial $k\n";
    # 5 fold
    for (my $i=1;$i<=5;$i++)
    {
        # my $curDir = "$tmpModelDir/time$k/fold$i";
        my $homeDir = "$outdir/$crmName/time$k/fold$i";
        # `mkdir $curDir` unless (-e "$curDir");
        # copy data from master node to slave node
        # `cp $homeDir/test.crm.fasta $curDir/test.crm.fasta`;
        # `cp $homeDir/test.neg.fasta $curDir/test.neg.fasta`;
        # `cp $homeDir/train.crm.fasta $curDir/train.crm.fasta`;
        # `cp $homeDir/train.neg.fasta $curDir/train.neg.fasta`;
        
        # generate kmer features
        `perl $kmer2lib $homeDir/test.crm.fasta $homeDir/test.neg.fasta $homeDir`;
        `perl $kmer2lib $homeDir/train.crm.fasta $homeDir/train.neg.fasta $homeDir`;
        # `mv $curDir/* $homeDir/`;
    }
}



