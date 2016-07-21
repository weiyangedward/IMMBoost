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
die `pod2text $0` if (@ARGV!=4);

# script to generate normalized kmer feature in libsvm format given CRM and neg fasta files
my $kmer2lib = "$Bin/kmer2lib.pl";

my $crmName = $ARGV[0]; # CRMname
my $outdir = $ARGV[1];
my $times = $ARGV[2];
my $nfolds = $ARGV[3];

# 10 trials
for (my $k=1;$k<=$times;$k++)
{

    warn "time $k...\n";
    # 5 fold
    for (my $i=1;$i<=$nfolds;$i++)
    {
        my $homeDir = "$outdir/$crmName/time$k/fold$i";
        # generate kmer features
        `perl $kmer2lib $homeDir/test.crm.fasta $homeDir/test.neg.fasta $homeDir`;
        `perl $kmer2lib $homeDir/train.crm.fasta $homeDir/train.neg.fasta $homeDir`;

    }
}



