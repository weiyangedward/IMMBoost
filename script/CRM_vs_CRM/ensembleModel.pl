=head1 Description
  
  This script combines the prediction from msIMM, 
  msIMM-RF and msIMM-SVM to generate an ensemble 
  model, which normalizes and averages predictions 
  from all models and make a final predicton.

=head1 Usage

  perl ensembleModel.pl CRMList Outdir

=cut

use strict;
use warnings;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use FindBin qw($Bin);

# die "Usage: perl $0 mapping outdir\n" unless @ARGV==2;

die `pod2text $0` if (@ARGV!=4);

my $crm = $ARGV[0];
my $outdir = $ARGV[1];
my $times = $ARGV[2];
my $nfolds = $ARGV[3];


# 10trials
for (my $k=1;$k<=$times;$k++)
{
    warn "time $k...\n";
    # 5folds
    for (my $i=1;$i<=$nfolds;$i++)
    {
        my $homeDir = "$outdir/$crm/time$k/fold$i";

        ##===== use equal weight for predictions from each model =====
        my $immTrainAUC = 1;
        # my $kmerSVMTrainAUC = 1;
        my $rfTrainAUC = 1;
        my $svmTrainAUC = 1;

        ##=== read the order of RF's test data ===
        my %id2pos = ();
        open IN,"$homeDir/test.ensembFeat.filGroup2" or die "cannot open $homeDir/test.ensembFeat.filGroup2";
        my $pos = 0;
        while (my $line = <IN>)
        {
            chomp($line);
            next if $line =~ /^crm/;
            my ($id,$score) = split(/\s+/,$line,2);
            $id2pos{$id} = $pos;
            $pos ++;
        }
        ##=== sync the order of RF's pred to input test data ====
        open RF,"$homeDir/test.ensembFeat.pred" or die "cannot open $homeDir/test.ensembFeat.pred\n";
        my @rf2pred = ();
        while (my $line = <RF>)
        {
            my ($lab, $score) = split /\s+/,$line;
            push @rf2pred,$score;
        }
        close RF;
        ##======= sync the order of SVM's pred ========##
        open SVM,"$homeDir/test.ensembFeat.lib.scaled.svm.pred.label" or die "cannot open $homeDir/test.ensembFeat.lib.scaled.svm.pred.label\n";
        my @svm2pred = ();
        while (my $line = <SVM>)
        {
            my ($lab, $score) = split /\s+/,$line;
            push @svm2pred,$score;
        }
        close SVM;
        ##======== sync RF's and SVM's order to test label ========
        my @normRFpred = norm(\@rf2pred);
        my @normSVMpred = norm(\@svm2pred);
        my @rerangeSVMpred = ();
        my @rerangeRFpred = ();
        my @testLab = ();
        open IN,"$homeDir/test.label" or die "cannot open $homeDir/test.label";
        while (my $line = <IN>)
        {
            chomp($line);
            my ($id, $lab) = split /\s+/,$line;
            push @testLab,$lab;
            if (exists $id2pos{$id}){
                push @rerangeSVMpred,$normSVMpred[$id2pos{$id}];
                push @rerangeRFpred,$normRFpred[$id2pos{$id}];
            }
        }
        close IN;
        # normalize msIMM pred
        my @imm2pred = ();
        open IMM,"$homeDir/test.crm.and.neg.IMM.pred.lab" or die "$homeDir/test.crm.and.neg.IMM.pred.lab";
        while (my $line = <IMM>)
        {
            my ($id, $score, $lab) = split /\s+/,$line;
            push @imm2pred,$score;
        }
        close IMM;
        my @normIMMpred = norm(\@imm2pred);
        # normalize kmerSVM pred
        # my @kmerSVM2pred = ();
        # open SVM,"$homeDir/kmerSVM.test.pred";
        # while (my $line = <SVM>)
        # {
        #    my ($lab, $score) = split /\s+/,$line;
        #    push @kmerSVM2pred,$score;
        # }
        # close SVM;
        # my @normkmerSVMpred = norm(\@kmerSVM2pred);

        open OUT,">$homeDir/blend.test.pred";
        for (my $i=0; $i<=$#testLab; $i++)
        {
            print OUT "$testLab[$i] ";
            my $newPredWithoutKmerSVM = $immTrainAUC * $normIMMpred[$i] + $rfTrainAUC * $rerangeRFpred[$i] + $svmTrainAUC * $rerangeSVMpred[$i];
            print OUT "$newPredWithoutKmerSVM\n";
        }
        close OUT;
        `Rscript $Bin/auc.R $homeDir/blend.test.pred $homeDir/blend.test.pred.auc`;
    }
}

# average AUC over 10trirals x 5folds and output to file
open OUT,">$outdir/$crm/ensembleModel.average.auc";
my $sumAUC = 0;
for (my $k=1;$k<=$times;$k++)
{

    for (my $i=1;$i<=$nfolds;$i++)
    {
        my $homeDir ="$outdir/$crm/time$k/fold$i";
        open AUC,"$homeDir/blend.test.pred.auc" or die "cannot open $homeDir/blend.test.pred.auc\n";
        while (<AUC>)
        {
            chomp(my $second = <AUC>);
            my $aucValue = (split /\s+/,$second)[-1];
            $sumAUC = $sumAUC + $aucValue;
        }
        close AUC;
    }
}
my $averageAUC = $sumAUC / ($nfolds * $times);
print OUT "$averageAUC\n";
close OUT;


sub norm {
    my $array = shift;
    my @new = ();
    my $max = max @{$array};
    my $min = min @{$array};
    my $len = $max - $min;
    for my $val (@{$array}){
        my $newVal = ($val - $min) / $len;
        push @new,$newVal;
    }
    return @new;
}
