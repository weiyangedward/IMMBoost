=head1

This script is used to train msIMMs and score
test seq in a given CRM set. A total of 10 trials x 5 fold
cross validation are performed, and an average
AUC over these 50 test is computed. This result
can be used as a baseline to evaluate ensemble models.
Note that BioPerl is needed is this script.

=cut

use strict;
use warnings;
use lib '/shared-mounts/sinhas/lib/';
use FindBin qw($Bin);
use Bio::SeqIO;
use File::Basename;

my $pred = "/home/weiyang4/software/SCRMshaw/bin/imm_score";
my $train = "/home/weiyang4/software/SCRMshaw/bin/imm_build";

die "Usage:\n\tperl $0 CRMname fileDir\n" unless @ARGV==2;

# create dir at slave nodes
my $tmpDir = $ENV{TMPDIR};
my $tmpModelDir = "$tmpDir/output";
`mkdir $tmpModelDir` unless (-e "$tmpModelDir");


my $crmName = $ARGV[0]; # CRMname
my $fileDir = $ARGV[1]; # outdir

my $crmDir = "$fileDir/$crmName";

##========== 
# 10trials x 5fold cross validation
# for each cross validation, a msIMM model is trained,
# then score test seq, and a AUC is computed
## ===========
# 10 trials 
for (my $k=1;$k<=10;$k++)
{
    `mkdir $tmpModelDir/time$k` unless (-e "$tmpModelDir/time$k");
    warn "$k ....\n";
    # 5 fold
    for (my $i=1;$i<=5;$i++)
    {
        my $curDir = "$tmpModelDir/time$k/fold$i";
        `mkdir $curDir` unless (-e "$curDir");
        my $homeDir = "$crmDir/time$k/fold$i";
        
        # copy data from master node to slave nodes
        `cp $homeDir/test.label $curDir/test.label`; # test data label
        `cp $homeDir/train.neg.fasta $curDir/train.neg.fasta`; # neg training data
        `cp $homeDir/train.crm.fasta $curDir/train.crm.fasta`; # pos training data
        `cp $homeDir/test.crm.and.neg.fasta $curDir/test.crm.and.neg.fasta`; # test data

        #===== hash test data label =====
        my %testID = ();
        my $testLabel = "$curDir/test.label";
        open IN,$testLabel;
        while (my $line = <IN>)
        {
            chomp $line;
            my ($id, $lab) = split /\s+/,$line;
            $testID{$id} = $lab; # hash test data label
        }
        close IN;
        #====== Train msIMM and score test seq ======
        `$train -r -k 6 < $curDir/train.neg.fasta > $curDir/train.msNeg.model`;    ## filGroup
        `$train -r -k 6 < $curDir/train.crm.fasta > $curDir/train.msCRM.model`;
        `$pred -f -n $curDir/train.msNeg.model $curDir/train.msCRM.model $curDir/test.crm.and.neg.fasta > $curDir/test.crm.and.neg.IMM.pred`;
        #====== output file with CRM id, score, label
        open OUT,">$curDir/test.crm.and.neg.IMM.pred.lab";
        open PRED,"$curDir/test.crm.and.neg.IMM.pred";
        while (<PRED>)
        {
            chomp(my $line = $_);
            my ($id, $score) = split /\s+/,$line;
            print OUT "$id\t$score\t$testID{$id}\n";
        }
        close PRED;
        close OUT;
        # compute AUC given the pred score
        `Rscript $Bin/auc.R $curDir/test.crm.and.neg.IMM.pred.lab $curDir/test.crm.and.neg.IMM.pred.lab.auc`;
        # copy AUC and pred files back to master node
        `cp $curDir/test.crm.and.neg.IMM.pred.lab.auc $homeDir/test.crm.and.neg.IMM.pred.lab.auc`;
        `cp $curDir/test.crm.and.neg.IMM.pred.lab $homeDir/test.crm.and.neg.IMM.pred.lab`;
    }
}
##===== 
# Compute average AUC
# sum of AUC among 10 trials x 5 folds,
# and divided by 50
##=======
my $sumAUC = 0;
# 10 trials
for (my $k=1;$k<=10;$k++)
{
    # 5 fold
    for (my $i=1;$i<=5;$i++)
    {
        my $curDir = "$tmpModelDir/time$k/fold$i";
        open AUC,"$curDir/test.crm.and.neg.IMM.pred.lab.auc";
        while (<AUC>)
        {
            chomp(my $second = <AUC>);
            my $aucValue = (split /\s+/,$second)[-1];
            $sumAUC = $sumAUC + $aucValue;
        }
        close AUC;
    }
}
my $aveAUC = $sumAUC / 50;  ## a total of 10x5 AUCs
# output average AUC to file
open OUT,">$tmpModelDir/IMM.average.auc";
print OUT "$aveAUC\n";
close OUT;
`cp $tmpModelDir/IMM.average.auc $crmDir/IMM.average.auc`;

