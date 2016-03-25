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

# if (@ARGV!=2){
#     die "Usage: perl $0 CRMname outdir\n";
# }

die `pod2text $0` if (@ARGV!=2);

# create an output dir in slave node
# my $tmpDir = $ENV{TMPDIR};
# my $tmpModelDir = "$tmpDir/output";
# `mkdir $tmpModelDir` unless (-e "$tmpModelDir");

my $crm = $ARGV[0]; # CRMname
my $outdir = $ARGV[1]; # outdir


sub pred {
    # 10trials
    for (my $k=1;$k<=10;$k++) 
    {
        # warn "trial $k\n";
        # `mkdir $tmpModelDir/time$k` unless (-e "$tmpModelDir/time$k");
        # 5folds
        for (my $i=1;$i<=5;$i++) 
        {    
            # my $curDir = "$tmpModelDir/time$k/fold$i";
            # `mkdir $curDir` unless (-e "$curDir");
            my $homeDir = "$outdir/$crm/time$k/fold$i";
            # copy data from master node to slave node
            # `cp $homeDir/train.ensembFeat.filGroup2.lib $curDir/train.ensembFeat.filGroup2.lib`;
            # `cp $homeDir/test.ensembFeat.filGroup2.lib $curDir/test.ensembFeat.filGroup2.lib`;
            # `cp $homeDir/train.ensembFeat.filGroup2 $curDir/train.ensembFeat.filGroup2`;
            # `cp $homeDir/test.ensembFeat.filGroup2 $curDir/test.ensembFeat.filGroup2`;

            # generate label files
            `awk \'{print \$1}\' $homeDir/train.ensembFeat.filGroup2.lib > $homeDir/train.ensembFeat.lib.lab`;
            `awk \'{print \$1}\' $homeDir/test.ensembFeat.filGroup2.lib > $homeDir/test.ensembFeat.lib.lab`;
            # train randomforest and pred on test data
            `Rscript $Bin/randomForest.R $homeDir/train.ensembFeat.filGroup2 $homeDir/train.ensembFeat.lib.lab $homeDir/test.ensembFeat.filGroup2 $homeDir/test.ensembFeat.lib.lab $homeDir/test.ensembFeat.pred $homeDir/train.ensembFeat.modelImportance`;
            # compute AUC scores
            `Rscript $Bin/auc.R $homeDir/test.ensembFeat.pred $homeDir/test.ensembFeat.pred.auc`;
            # `cp $homeDir/test.ensembFeat.pred.auc $homeDir/test.ensembFeat.pred.auc`;
            # `mv $homeDir/test.ensembFeat.pred $homeDir/test.ensembFeat.pred`;
            # `mv $homeDir/train.ensembFeat.modelImportance $homeDir/train.ensembFeat.modelImportance`;
        }
    }

    ##===== average AUC over 10trials x 5folds =======##
    open OUT,">$outdir/$crm/IMM_RF.average.auc";
    my $sumAUC = 0;
    # 10trials
    for (my $k=1;$k<=10;$k++)
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
    my $averageAUC = $sumAUC / 50; # average over 10 trials x 5-folds
    print OUT "$averageAUC\n";
    close OUT;
    # `mv $tmpModelDir/RF.average.auc $outdir/$crm/RF.average.auc`; 
}

# call pred() function
pred();

