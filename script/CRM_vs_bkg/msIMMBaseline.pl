=head1 Description

    This script is used to train msIMMs and score
    test seq in a given CRM set. A total of 10trials 
    x 5fold cross validation are performed, and an 
    average AUC over these 50 test is computed. This 
    result can be used as a baseline to evaluate 
    ensemble models. Note that BioPerl is needed is 
    this script.

=head1 Usage

    perl msIMMBaseline.pl CRMname outdir

=cut

use strict;
use warnings;
use FindBin qw($Bin);
use Bio::SeqIO;
use File::Basename;

my $train = "$Bin/../../src/imm/bin/imm_build";
my $pred = "$Bin/../../src/imm/bin/imm_score";


die `pod2text $0` if (@ARGV!=4);


my $crmName = $ARGV[0]; # CRMname
my $outdir = $ARGV[1]; # outdir
my $times = $ARGV[2];
my $nfolds = $ARGV[3];

my $crmDir = "$outdir/$crmName";

##========== 
# 10trials x 5fold cross validation
# for each cross validation, a msIMM model is trained,
# then score test seq, and a AUC is computed
## ===========
# 10 trials 
for (my $k=1;$k<=$times;$k++)
{
    # `mkdir $tmpModelDir/time$k` unless (-e "$tmpModelDir/time$k");
    warn "time $k ...\n";
    # 5 fold
    for (my $i=1;$i<=$nfolds;$i++)
    {
        # my $curDir = "$tmpModelDir/time$k/fold$i";
        # `mkdir $curDir` unless (-e "$curDir");
        my $homeDir = "$crmDir/time$k/fold$i";

        #===== hash test data label =====
        my %testID = ();
        my $testLabel = "$homeDir/test.label";
        open IN,$testLabel or die "cannot open $testLabel";
        while (my $line = <IN>)
        {
            chomp $line;
            my ($id, $lab) = split /\s+/,$line;
            $testID{$id} = $lab; # hash test data label
        }
        close IN;
        #====== Train msIMM and score test seq ======
        `$train -r -k 6 < $homeDir/train.neg.fasta > $homeDir/train.msNeg.model`;    ## filGroup
        `$train -r -k 6 < $homeDir/train.crm.fasta > $homeDir/train.msCRM.model`;
        `$pred -f -n $homeDir/train.msNeg.model $homeDir/train.msCRM.model $homeDir/test.crm.and.neg.fasta > $homeDir/test.crm.and.neg.IMM.pred`;
        #====== output file with CRM id, score, label
        open OUT,">$homeDir/test.crm.and.neg.IMM.pred.lab";
        open PRED,"$homeDir/test.crm.and.neg.IMM.pred" or die "cannot open $homeDir/test.crm.and.neg.IMM.pred";
        while (<PRED>)
        {
            chomp(my $line = $_);
            my ($id, $score) = split /\s+/,$line;
            print OUT "$id\t$score\t$testID{$id}\n";
        }
        close PRED;
        close OUT;
        # compute AUC given the pred score
        `Rscript $Bin/IMMauc.R $homeDir/test.crm.and.neg.IMM.pred.lab $homeDir/test.crm.and.neg.IMM.pred.lab.auc`;
    }
}
##===== 
# Compute average AUC
# sum of AUC among 10 trials x 5 folds,
# and divided by 50
##=======
my $sumAUC = 0;
# 10 trials
for (my $k=1;$k<=$times;$k++)
{
    # 5 fold
    for (my $i=1;$i<=$nfolds;$i++)
    {
        # my $homeDir = "$tmpModelDir/time$k/fold$i";
        my $homeDir = "$crmDir/time$k/fold$i";
        open AUC,"$homeDir/test.crm.and.neg.IMM.pred.lab.auc" or die "cannot open $homeDir/test.crm.and.neg.IMM.pred.lab.auc";
        while (<AUC>)
        {
            chomp(my $second = <AUC>);
            my $aucValue = (split /\s+/,$second)[-1];
            $sumAUC = $sumAUC + $aucValue;
        }
        close AUC;
    }
}
my $aveAUC = $sumAUC / ($nfolds*$times);  ## a total of 10x5 AUCs
# output average AUC to file
open OUT,">$crmDir/IMM.average.auc";
print OUT "$aveAUC\n";
close OUT;

