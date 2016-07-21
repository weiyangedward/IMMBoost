=head1 Description

    This script trains svm with kmer features and 
    predict on test data. AUC will be computed for
    the predict results. Note that libsvm library 
    is needed.

=head1 Usage

    perl kmerSVM.pl CRMname Outdir

=cut

use strict;
use warnings;
use File::Basename;
use FindBin qw($Bin);

die `pod2text $0` if (@ARGV!=3);


# my $kmer2lib = "$Bin/kmer2lib.pl"; # generate kmer features
my $sub = "$Bin/../../src/libsvm-3.21/tools/subset.py"; # sample a subset of training data
my $scale = "$Bin/../../src/libsvm-3.21/svm-scale"; # scale test and training data features
my $libsvm = "$Bin/libsvmPval.py";

my $crm = $ARGV[0];
my $outdir = $ARGV[1];
my $times = $ARGV[2];


# warn "$crm\n";
# 10 trials
for (my $k=1;$k<=$times;$k++)
{
   # 5fold
    warn "time $k...\n";
    for (my $i=1;$i<=5;$i++)
    {
        my $homeDir = "$outdir/$crm/time$k/fold$i";

        # sample a subset (50) of training data to find best parameter
        `python $sub $homeDir/train.crm.fasta.andNeg.rc.6mer 10 $homeDir/train.crm.fasta.andNeg.rc.6mer.sub1 $homeDir/train.crm.fasta.andNeg.rc.6mer.sub2`;
        # train svm and pred on test data
        `python $libsvm $homeDir/train.crm.fasta.andNeg.rc.6mer $homeDir/test.crm.fasta.andNeg.rc.6mer $homeDir/train.crm.fasta.andNeg.rc.6mer.model $homeDir/test.crm.fasta.andNeg.rc.6mer.pred $homeDir/train.crm.fasta.andNeg.rc.6mer.sub1 2>$homeDir/libsvm.log`;

        ##=== store label for test data =====
        my @label = ();
        # `cp $curDir/test.label $homeDir/test.label`;
        open LAB,"$homeDir/test.label" or die "cannot open $homeDir/test.label";
        while (<LAB>)
        {
            chomp(my $line = $_);
            my @array = split /\s+/,$line;
            push @label,$array[1];
        }
        close LAB;
        ##=== store predicted score ====
        my @predValue = ();
        open PRED,"$homeDir/test.crm.fasta.andNeg.rc.6mer.pred" or die "cannot open $homeDir/test.crm.fasta.andNeg.rc.6mer.pred";
        while (<PRED>)
        {
            chomp(my $line = $_);
            my $conf = $1 if $line =~ /\[(\S+?)\]/;
            push @predValue,$conf;
        }
        close PRED;
        ##== output test label and predicted score ====
        open OUT1,">$homeDir/test.crm.fasta.andNeg.rc.6mer.pred.lab" or die "cannot open $homeDir/test.crm.fasta.andNeg.rc.6mer.pred.lab";
        for (my $j=0;$j<=$#label;$j++)
        {
            print OUT1 "$label[$j] $predValue[$j]\n";
        }
        close OUT1;
        ##==== compute auc score ====
        `Rscript $Bin/auc.R $homeDir/test.crm.fasta.andNeg.rc.6mer.pred.lab $homeDir/test.crm.fasta.andNeg.rc.6mer.pred.lab.auc`;
    }
}
##===== averaging auc score =====
open OUT,">$outdir/$crm/kmerSVM.average.auc";
my $sumAUC = 0;
# 10trials
for (my $k=1;$k<=$times;$k++)
{
    # 5folds
    for (my $i=1;$i<=5;$i++)
    {
        my $homeDir = "$outdir/$crm/time$k/fold$i";
        open AUC,"$homeDir/test.crm.fasta.andNeg.rc.6mer.pred.lab.auc" or die "cannot open $homeDir/test.crm.fasta.andNeg.rc.6mer.pred.lab.auc\n";
        while (<AUC>)
        {
            chomp(my $second = <AUC>);
            my $aucValue = (split /\s+/,$second)[-1];
            $sumAUC = $sumAUC + $aucValue;
        }
        close AUC;
    }
}
# average AUC over 10trials x 5folds
my $averageAUC = $sumAUC / (5*$times);
print OUT "$averageAUC\n";
close OUT;


