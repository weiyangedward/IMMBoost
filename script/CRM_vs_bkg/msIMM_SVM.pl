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

die `pod2text $0` if (@ARGV!=3);

my $libsvm = "$Bin/libsvmPval.py";
my $liblinear = "$Bin/libPval.py";
my $scale = "$Bin/../../src/libsvm-3.21/svm-scale";
my $sub = "$Bin/../../src/libsvm-3.21/tools/subset.py";

my $crm = $ARGV[0];
my $outdir = $ARGV[1];
my $times = $ARGV[2];

# 10trials
for (my $k=1;$k<=$times;$k++)
{

    warn "time $k ...\n";
    # 5folds
    for (my $i=1;$i<=5;$i++)
    {
        my $homeDir = "$outdir/$crm/time$k/fold$i";
        
        # create scale range from training data
        `$scale -s $homeDir/train.ensembFeat.lib.range $homeDir/train.ensembFeat.filGroup2.lib > $homeDir/train.ensembFeat.lib.scaled`;
        `$sub $homeDir/train.ensembFeat.lib.scaled 10 $homeDir/train.ensembFeat.lib.scaled.sub1 $homeDir/train.ensembFeat.lib.scaled.sub2`;
        # scale test data using scale range
        `$scale -r $homeDir/train.ensembFeat.lib.range $homeDir/test.ensembFeat.filGroup2.lib > $homeDir/test.ensembFeat.lib.scaled`;

        `python $liblinear $homeDir/train.ensembFeat.lib.scaled $homeDir/test.ensembFeat.lib.scaled $homeDir/train.ensembFeat.lib.scaled.model $homeDir/test.ensembFeat.lib.scaled.pred.confidentScore $homeDir/train.ensembFeat.lib.scaled.sub1 > $homeDir/liblinear.log`;
        # `python $libsvm $homeDir/train.ensembFeat.lib.scaled $homeDir/test.ensembFeat.lib.scaled $homeDir/train.ensembFeat.lib.scaled.model $homeDir/test.ensembFeat.lib.scaled.pred.confidentScore $homeDir/train.ensembFeat.lib.scaled.sub1 > $homeDir/liblinear.log`;

        # store label and pred score
        my @label = ();
        open TEST,"$homeDir/test.ensembFeat.filGroup2.lib" or die "cannot open $homeDir/test.ensembFeat.filGroup2.lib";
        while (<TEST>)
        {
            chomp(my $line = $_);
            my @array = split /\s+/,$line;
            push @label,$array[0];
        }
        close TEST;
        my @predValue = ();
        open PRED,"$homeDir/test.ensembFeat.lib.scaled.pred.confidentScore" or die "cannot open $homeDir/test.ensembFeat.lib.scaled.pred.confidentScore";
        while (<PRED>)
        {
            chomp(my $line = $_);
            my $conf = $1 if $line =~ /\[(\S+?)\]/;
            push @predValue,$conf;
        }
        close PRED;
        # output label and pred score to file
        open OUT1,">$homeDir/test.ensembFeat.lib.scaled.svm.pred.label" or die "$homeDir/test.ensembFeat.lib.scaled.svm.pred.label";
        for (my $j=0;$j<=$#label;$j++){
            print OUT1 "$label[$j] $predValue[$j]\n";
        }
        close OUT1;
        # compute AUC
        `Rscript $Bin/auc.R $homeDir/test.ensembFeat.lib.scaled.svm.pred.label $homeDir/test.ensembFeat.lib.scaled.svm.pred.label.auc`;
    }
}

# average over AUC on 10trials x 5folds
open OUT,">$outdir/$crm/IMM_SVM.average.auc";
my $sumAUC = 0;
for (my $k=1;$k<=$times;$k++)
{
    for (my $i=1;$i<=5;$i++)
    {
        my $homeDir = "$outdir/$crm/time$k/fold$i";
        open AUC,"$homeDir/test.ensembFeat.lib.scaled.svm.pred.label.auc" or die "cannot open $homeDir/test.ensembFeat.lib.scaled.svm.pred.label.auc\n";
        while (<AUC>)
        {
            chomp(my $second = <AUC>);
            my $aucValue = (split /\s+/,$second)[-1];
            $sumAUC = $sumAUC + $aucValue;
        }
        close AUC;
    }
}
my $averageAUC = $sumAUC / (5*$times);
print OUT "$averageAUC\n";
close OUT;
