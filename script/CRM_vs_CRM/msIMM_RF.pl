=head1 Description
  
  This script train a randomforest model on training 
  data with msIMM features and pred on test data.
  AUC scores will also be computed for the predict results.

=head1 Usage

  perl msIMM_RF.pl CRMList Outdir

=cut

use strict;
use warnings;
use File::Basename;
use FindBin qw($Bin);


die `pod2text $0` if (@ARGV!=3);


my $crm = $ARGV[0]; # CRMname
my $outdir = $ARGV[1]; # outdir
my $times = $ARGV[2];


sub pred {
    # 10trials
    for (my $k=1;$k<=$times;$k++) 
    {
        warn "time $k...\n";
        # 5folds
        for (my $i=1;$i<=5;$i++) 
        {    

            my $homeDir = "$outdir/$crm/time$k/fold$i";


            # generate label files
            `awk \'{print \$1}\' $homeDir/train.ensembFeat.filGroup2.lib > $homeDir/train.ensembFeat.lib.lab`;
            `awk \'{print \$1}\' $homeDir/test.ensembFeat.filGroup2.lib > $homeDir/test.ensembFeat.lib.lab`;
            # train randomforest and pred on test data
            `Rscript $Bin/randomForest.R $homeDir/train.ensembFeat.filGroup2 $homeDir/train.ensembFeat.lib.lab $homeDir/test.ensembFeat.filGroup2 $homeDir/test.ensembFeat.lib.lab $homeDir/test.ensembFeat.pred $homeDir/train.ensembFeat.modelImportance`;
            # compute AUC scores
            `Rscript $Bin/auc.R $homeDir/test.ensembFeat.pred $homeDir/test.ensembFeat.pred.auc`;
        }
    }

    ##===== average AUC over 10trials x 5folds =======##
    open OUT,">$outdir/$crm/IMM_RF.average.auc";
    my $sumAUC = 0;
    # 10trials
    for (my $k=1;$k<=$times;$k++)
    {
        # 5folds
        for (my $i=1;$i<=5;$i++)
        { 
            my $homeDir = "$outdir/$crm/time$k/fold$i";
            open AUC,"$homeDir/test.ensembFeat.pred.auc" or die "cannot open $homeDir/test.ensembFeat.pred.auc\n";
            while (<AUC>)
            {
                chomp(my $second = <AUC>);
                my $aucValue = (split /\s+/,$second)[-1];
                $sumAUC = $sumAUC + $aucValue;
            }
            close AUC;
        }
    }
    my $averageAUC = $sumAUC / (5*$times); # average over 10 trials x 5-folds
    print OUT "$averageAUC\n";
    close OUT;
    # `mv $tmpModelDir/RF.average.auc $outdir/$crm/RF.average.auc`; 
}

# call pred() function
pred();

