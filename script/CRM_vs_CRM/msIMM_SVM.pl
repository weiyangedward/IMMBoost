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

# die "Usage: perl $0 mapping IMMmodelDir\n" unless @ARGV==2;
die `pod2text $0` if (@ARGV!=2);

#my $train = "/home/n-z/weiyang4/bin/liblinear-1.95/train";
#my $pred = "/home/n-z/weiyang4/bin/liblinear-1.95/predict";

# my $tmpDir = $ENV{TMPDIR};
# my $tmpModelDir = "$tmpDir/output";
# `mkdir $tmpModelDir` unless (-e "$tmpModelDir");

my $libsvm = "$Bin/libsvmPval.py";
my $liblinear = "$Bin/libPval.py";
my $scale = "$Bin/../../src/libsvm-3.21/svm-scale";
my $sub = "$Bin/../../src/libsvm-3.21/tools/subset.py";

my $crm = $ARGV[0];
my $outdir = $ARGV[1];


for (my $k=1;$k<=10;$k++)
{
    # `mkdir $tmpModelDir/time$k` unless (-e "$tmpModelDir/time$k");
    for (my $i=1;$i<=5;$i++)
    {
        # my $homeDir = "$tmpModelDir/time$k/fold$i";
        # `mkdir $homeDir` unless (-e "$homeDir");
        my $homeDir = "$outdir/$crm/time$k/fold$i";
        ##====== Copy files from home dir to tmp dir =======##
        # `cp $homeDir/train.ensembFeat.filGroup2.lib $homeDir/train.ensembFeat.filGroup2.lib`;
        # `cp $homeDir/test.ensembFeat.filGroup2.lib $homeDir/test.ensembFeat.filGroup2.lib`;
        # `cp $homeDir/train.ensembFeat.filGroup2.Dmel.lib $homeDir/train.ensembFeat.filGroup2.Dmel.lib`;
        ##====== prepare inputs for liblinear =========##
        `$scale -s $homeDir/train.ensembFeat.lib.range $homeDir/train.ensembFeat.filGroup2.lib > $homeDir/train.ensembFeat.lib.scaled`;
        # `$sub $homeDir/train.ensembFeat.lib.scaled 50 $homeDir/train.ensembFeat.lib.scaled.sub1 $homeDir/train.ensembFeat.lib.scaled.sub2`;
        `$scale -r $homeDir/train.ensembFeat.lib.range $homeDir/test.ensembFeat.filGroup2.lib > $homeDir/test.ensembFeat.lib.scaled`;
        `$scale -r $homeDir/train.ensembFeat.lib.range $homeDir/train.ensembFeat.filGroup2.Dmel.lib > $homeDir/train.ensembFeat.filGroup2.Dmel.lib.scaled`;
        ##======== liblinear prediction ========##
        #`python $liblinear $homeDir/train.ensembFeat.lib.scaled $homeDir/test.ensembFeat.lib.scaled $homeDir/train.ensembFeat.lib.scaled.model $homeDir/test.ensembFeat.lib.scaled.pred.confidentScore $homeDir/train.ensembFeat.filGroup2.Dmel.lib.scaled > $homeDir/liblinear.log`;
        `python $libsvm $homeDir/train.ensembFeat.lib.scaled $homeDir/test.ensembFeat.lib.scaled $homeDir/train.ensembFeat.lib.scaled.model $homeDir/test.ensembFeat.lib.scaled.pred.confidentScore $homeDir/train.ensembFeat.filGroup2.Dmel.lib.scaled > $homeDir/liblinear.log`;
        my @label = ();
        open TEST,"$homeDir/test.ensembFeat.filGroup2.lib";
        while (<TEST>){
            chomp(my $line = $_);
            my @array = split /\s+/,$line;
            push @label,$array[0];
        }
        close TEST;
        my @predValue = ();
        open PRED,"$homeDir/test.ensembFeat.lib.scaled.pred.confidentScore";
        while (<PRED>){
            chomp(my $line = $_);
            my $conf = $1 if $line =~ /\[(\S+?)\]/;
            push @predValue,$conf;
        }
        close PRED;
        open OUT1,">$homeDir/test.ensembFeat.lib.scaled.svm.pred.label";
        for (my $j=0;$j<=$#label;$j++){
            print OUT1 "$label[$j] $predValue[$j]\n";
        }
        close OUT1;
        `Rscript $Bin/auc.R $homeDir/test.ensembFeat.lib.scaled.svm.pred.label $homeDir/test.ensembFeat.lib.scaled.svm.pred.label.auc`;
        # `cp $homeDir/* $homeDir/`;
    }
}

open OUT,">$outdir/$crm/IMM_SVM.average.auc";
my $sumAUC = 0;
for (my $k=1;$k<=10;$k++){
    for (my $i=1;$i<=5;$i++){
        # my $homeDir = "$tmpModelDir/time$k/fold$i";
        my $homeDir = "$outdir/$crm/time$k/fold$i";
        open AUC,"$homeDir/test.ensembFeat.lib.scaled.svm.pred.label.auc" or die "cannot open $homeDir/test.ensembFeat.lib.scaled.svm.pred.label.auc\n";
        while (<AUC>){
            chomp(my $second = <AUC>);
            my $aucValue = (split /\s+/,$second)[-1];
            $sumAUC = $sumAUC + $aucValue;
        }
        close AUC;
    }
}
my $averageAUC = $sumAUC / 50;
print OUT "$averageAUC\n";
close OUT;
# `mv $tmpModelDir/SVM.average.auc $modelDir/$crm/SVM.average.auc`;
