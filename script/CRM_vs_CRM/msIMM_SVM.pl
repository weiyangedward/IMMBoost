=head1 Description
  
  This script trains a linear SVM using liblinear 
  library on msIMM features (msIMM-SVM), and predict 
  on test data. AUC is computed for predict results.

=head1 Usage

  perl msIMM_SVM.pl CRMname outdir

=cut

use strict;
use warnings;
use File::Basename;
use FindBin qw($Bin);

die `pod2text $0` if (@ARGV!=4);

my $libsvm = "$Bin/libsvmPval.py";
my $liblinear = "$Bin/libPval.py";
my $scale = "$Bin/../../src/libsvm-3.21/svm-scale";
my $sub = "$Bin/../../src/libsvm-3.21/tools/subset.py";

my $crm = $ARGV[0];
my $outdir = $ARGV[1];
my $times = $ARGV[2];
my $nfolds = $ARGV[3];


for (my $k=1;$k<=$times;$k++)
{
    warn "time $k...\n";
    for (my $i=1;$i<=$nfolds;$i++)
    {
        # warn "$k, $i\n";
        my $homeDir = "$outdir/$crm/time$k/fold$i";
        ##======================
        # prepare inputs for svm
        #=======================
        `$scale -s $homeDir/train.ensembFeat.lib.range $homeDir/train.ensembFeat.lib > $homeDir/train.ensembFeat.lib.scaled`;
        `$scale -r $homeDir/train.ensembFeat.lib.range $homeDir/test.ensembFeat.lib > $homeDir/test.ensembFeat.lib.scaled`;
        `$scale -r $homeDir/train.ensembFeat.lib.range $homeDir/train.ensembFeat.Dmel.lib > $homeDir/train.ensembFeat.Dmel.lib.scaled`;
        ##======================================================
        # train IMM-SVM and predict on test data using liblinear
        #=======================================================
        `python $liblinear $homeDir/train.ensembFeat.lib.scaled $homeDir/test.ensembFeat.lib.scaled $homeDir/train.ensembFeat.lib.scaled.model $homeDir/test.ensembFeat.lib.scaled.pred.confidentScore $homeDir/train.ensembFeat.Dmel.lib.scaled 2>$homeDir/liblinear.log`;

        #====================================================
        # train IMM-SVM and predict on test data using libsvm
        #====================================================
        # warn "train IMM-SVM and predict on test data...\n";
        # `python $libsvm $homeDir/train.ensembFeat.lib.scaled $homeDir/test.ensembFeat.lib.scaled $homeDir/train.ensembFeat.lib.scaled.model $homeDir/test.ensembFeat.lib.scaled.pred.confidentScore $homeDir/train.ensembFeat.Dmel.lib.scaled 2>$homeDir/liblinear.log`;

        my @label = ();
        open TEST,"$homeDir/test.ensembFeat.lib" or die "cannot open $homeDir/test.ensembFeat.lib";
        while (<TEST>){
            chomp(my $line = $_);
            my @array = split /\s+/,$line;
            push @label,$array[0];
        }
        close TEST;

        my @predValue = ();
        open PRED,"$homeDir/test.ensembFeat.lib.scaled.pred.confidentScore" or die "cannot open $homeDir/test.ensembFeat.lib.scaled.pred.confidentScore";
        while (<PRED>){
            chomp(my $line = $_);
            my $conf = $1 if $line =~ /\[(\S+?)\]/;
            push @predValue,$conf;
        }
        close PRED;

        open OUT1,">$homeDir/test.ensembFeat.lib.scaled.svm.pred.label" or die "cannot open $homeDir/test.ensembFeat.lib.scaled.svm.pred.label";
        for (my $j=0;$j<=$#label;$j++){
            print OUT1 "$label[$j] $predValue[$j]\n";
        }
        close OUT1;

        #=============
        # compute AUC
        #=============
        `Rscript $Bin/auc.R $homeDir/test.ensembFeat.lib.scaled.svm.pred.label $homeDir/test.ensembFeat.lib.scaled.svm.pred.label.auc`;
    }
}

open OUT,">$outdir/$crm/IMM_SVM.average.auc";
my $sumAUC = 0;
for (my $k=1;$k<=$times;$k++)
{
    for (my $i=1;$i<=$nfolds;$i++)
    {
        my $homeDir = "$outdir/$crm/time$k/fold$i";
        open AUC,"$homeDir/test.ensembFeat.lib.scaled.svm.pred.label.auc" or die "cannot open $homeDir/test.ensembFeat.lib.scaled.svm.pred.label.auc\n";
        while (<AUC>)
        {
            chomp(my $second = <AUC>);
            my $aucValue = 0;
            eval{$aucValue = (split /\s+/,$second)[-1]}; die $@ if $@;
            $sumAUC = $sumAUC + $aucValue;
        }
        close AUC;
    }
}
my $averageAUC = $sumAUC / ($nfolds*$times);
print OUT "$averageAUC\n";
close OUT;
