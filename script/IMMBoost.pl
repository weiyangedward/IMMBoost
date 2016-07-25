=head1 Description:
  
  This script wrap the pipeline for IMMBoost (CRM_vs_bkg and CRM_vs_CRM)

=head1 Usage:

  perl IMMBoost.pl [options] CRMList Outdir Datadir crmGroupTable

  --task <str>      What task to perform? default=crm_vs_crm. There are two modes:
                          1) "--task crm_vs_bkg": classify CRM from background 
                                  genomic seq; 
                          2) "--task crm_vs_crm": classify CRM from other CRM seq
  
  --step <str>      What steps to run? default=12345678.
                        1. prepare data
                        2. IMM score feature
                        3. IMM prediction
                        4. IMM-SVM prediction
                        5. IMM-RF prediction 
                        6. generate kmer-SVM feature
                        7. kmer-SVM prediction
                        8. IMM-Ensemble
  --nfolds <int>    To perform n-fold cross validation. default=5.
  --ktimes <int>    To repeat n-fold cross validation for k times. default=2.

=head1 Example:

  perl IMMBoost.pl -task crm_vs_bkg CRMsetsList.txt ../sampleOutput/CRM_vs_bkg/ ../sampleData/CRMsets/ CRM.group.V3.txt

=cut 

use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;

my ($task, $step, $ktimes, $nfolds);

$task = "crm_vs_crm";
$nfolds = 5;
$ktimes = 2;
$step = "12345678";

GetOptions(
  "task:s"=>\$task,
  "step:s"=>\$step,
  "nfolds:i"=>\$nfolds,
  "ktimes:i"=>\$ktimes
);

die `pod2text $0` if (@ARGV<4);

my $crmLst = $ARGV[0]; # a list file contains path to the dir of CRM datasets
my $outdir = $ARGV[1]; # outdir
my $datadir = $ARGV[2];
my $crmGroupTable = $ARGV[3];

#==============================
# check if BioPerl is installed
#==============================
my $require_package = eval{
    require Bio::SeqIO;
    # require Bio::DB::Fasta;
    Bio::SeqIO->import();
    #  Bio::DB::Fasta->import();
    1;
};

if ($require_package){
    warn "Bio::SeqIO loaded and imported successfully\n";
    # warn "Bio::DB::Fasta loaded and imported successfully\n";
}
else{
    die("Bio::SeqIO cannot be imported. Please check to see if the modules have already been installed.\n");
}

#================================================================================
#
# Task: CRM_vs_CRM 
#
#================================================================================
sub CRM_vs_CRM {
  my @crmNames = ();
  open IN,$crmLst or die "cannot open $crmLst";
  while (<IN>)
  {
    chomp(my $indir = $_);
    my $crm = (split /\//,$indir)[-1];
    push @crmNames,$crm;
  }

  if ($step =~/1/)
  {
    if (-d "$outdir")
    {
      `rm -rf $outdir/*`;
    }
    #============================================
    # Step1:
    # for each CRM set, train a msIMM model and 
    # generate 10trials x 5folds cross validation 
    # test and training data
    #============================================
    warn "Step1 ----------------------------------\n";
    warn "prepare Model And Data ...\n";
    open IN,$crmLst or die "cannot open $crmLst";
    while (<IN>)
    {
      chomp(my $indir = $_);
      my $crm = (split /\//,$indir)[-1];
      `mkdir -p $outdir/$crm`;
      warn "train IMM model on $crm\n";
      `perl $Bin/CRM_vs_CRM/prepareModelAndData.pl $Bin/../$indir $outdir $ktimes $nfolds`;
    }
    close IN;
    warn "Step1 is done!\n";
  }

  if ($step =~/2/)
  {
    ##=============================================
    # Step2: 
    # Generate msIMM score as features for each seq
    ##=============================================
    warn "Step2 ----------------------------------\n";
    warn "Generate msIMM score features ...\n";
    for my $crm (@crmNames)
    {
      warn "$crm\n";
      for (my $k=1;$k<=$ktimes;$k++)
      {
        for (my $i=1;$i<=$nfolds;$i++)
        {
          warn "time$k fold$i\n";
          `perl $Bin/CRM_vs_CRM/generateIMMScoreFeature.pl $crm $outdir $k $i`;
        }
      }
    }
    warn "filter msCRMs in the same group...\n";
    `perl $Bin/CRM_vs_CRM/filterGroupCRM.pl $datadir $outdir $crmGroupTable $ktimes $nfolds`;

    warn "filter Dmel CRMs in the same group...\n";
    `perl $Bin/CRM_vs_CRM/filterGroupCRM.DmelTrainData.pl $datadir $outdir $crmGroupTable $ktimes $nfolds`;
    warn "Step2 is done!\n";
  }

  if ($step =~/3/)
  {
    warn "Step3 ----------------------------------\n";
    ##==============
    # Step3.1:
    # msIMM baseline
    ##==============
    warn "Step3.1 msIMM prediction ...\n";
    for my $crm (@crmNames){
      warn "$crm\n";
      `perl $Bin/CRM_vs_CRM/msIMMBaseline.pl $crm $outdir $ktimes $nfolds`;
    }
    warn "Step3.1 is done!\n";

    ##====================
    # Step3.2:
    # Summarize msIMM AUCs
    ##====================
    warn "Step3.2 Summarize msIMM AUCs...\n";
    combine_auc($outdir, \@crmNames, "summaryAUC_msIMMBaseline.txt", "IMM.average.auc");
    warn "Step3.2 is done!\n";
  }

  if ($step=~/4/)
  {
    warn "Step4 ----------------------------------\n";
    ##==================
    # Step4.1:
    # IMM-SVM prediction 
    ##==================
    warn "Step4.1 IMM-SVM prediction...\n";
    for my $crm (@crmNames)
    {
      warn "$crm\n";
      `perl $Bin/CRM_vs_CRM/msIMM_SVM.pl $crm $outdir $ktimes $nfolds`;
    }
    warn "Step4.1 is done!\n";


    ##============================
    # Step4.2:
    # Summarize IMM-SVM prediction 
    ##============================
    warn "Step4.2 Summarize IMM-SVM AUCs...\n";
    combine_auc($outdir, \@crmNames, "summaryAUC_IMM_SVM.txt", "IMM_SVM.average.auc");
    warn "Step4.2 is done!\n";
  }

  if ($step =~/5/)
  {
    warn "Step5 ----------------------------------\n";
    ##=================
    # Step5.1:
    # IMM-RF prediction 
    ##=================
    warn "Step5.1 IMM-RF prediction...\n";
    for my $crm (@crmNames){
      warn "$crm\n";
      `perl $Bin/CRM_vs_CRM/msIMM_RF.pl $crm $outdir $ktimes $nfolds`;
    }
    warn "Step5.1 is done!\n";

    ##=====================
    # Step5.2:
    # Summarize IMM-RF AUCs 
    ##=====================
    warn "Step5.2 Summarize IMM-RF AUCs...\n";
    combine_auc($outdir, \@crmNames, "summaryAUC_IMM_RF.txt", "IMM_RF.average.auc");
    warn "Step5.2 is done!\n";
  }

  if ($step =~/6/)
  {
    warn "Step6 ----------------------------------\n";
    ##==================================
    # Step6.1:
    # generate kmer features for kmerSVM 
    ##==================================
    warn "Step6.1 generate kmer features for kmerSVM...\n";
    for my $crm (@crmNames)
    {
      warn "$crm\n"; 
      `perl $Bin/CRM_vs_CRM/kmerSVMFeature.pl $crm $outdir $ktimes $nfolds`;
    }
    warn "Step6.1 is done!\n";
  }


  if ($step =~/7/)
  {
    warn "Step7 ----------------------------------\n";
    ##==================
    # Step7.1:
    # kmerSVM prediction 
    ##==================
    warn "Step7.1 kmerSVM prediction...\n";
    for my $crm (@crmNames)
    {
      warn "$crm\n";
      `perl $Bin/CRM_vs_CRM/kmerSVM.pl $crm $outdir $ktimes $nfolds`;
    }
    warn "Step7.1 is done!\n";

    ##=====================
    # Step7.2:
    # Summarize kmerSVM AUC 
    ##=====================
    warn "Step7.2 Summarize kmerSVM AUC...\n";
    combine_auc($outdir, \@crmNames, "summaryAUC_kmerSVM.txt", "kmerSVM.average.auc");
    warn "Step7.2 is done!\n";
  }


  if ($step =~/8/)
  {
    warn "Step8 ----------------------------------\n";
    ##===================
    # Step8.1:
    # Ensemble prediction 
    ##===================
    warn "Step8.1 Ensemble prediction...\n";
    for my $crm (@crmNames)
    {
      warn "$crm\n";
      `perl $Bin/CRM_vs_CRM/ensembleModel.pl $crm $outdir $ktimes $nfolds`;
    }
    warn "Step8.1 is done!\n";


    ##===========================
    # Step8.2:
    # summarize ensembleModel AUC 
    ##===========================
    warn "Step8.2 summarize ensembleModel AUC...\n";
    combine_auc($outdir, \@crmNames, "summaryAUC_ensembleModel.txt", "ensembleModel.average.auc");
    warn "Step8.2 is done!\n";
  }

  warn "All steps are doen!\n";
}


#================================================================================
#
# Task: CRM_vs_bkg 
#
#================================================================================
sub CRM_vs_bkg {

  my @crmNames = (); # an array to store CRMnames
  open IN,$crmLst or die "cannot open $crmLst";
  while (<IN>)
  {
    chomp(my $indir = $_);
    my $crm = (split /\//,$indir)[-1];
    push @crmNames,$crm;
  }

  if ($step =~/1/)
  {
    #============================================
    # Step1:
    # for each CRM set, train a msIMM model and 
    # generate 10trials x 5folds cross validation 
    # test and training data
    #============================================
    warn "Step1 ----------------------------------\n";
    warn "prepare Model And Data...\n";

    open IN,$crmLst or die "cannot open $crmLst";
    while (<IN>)
    {
      chomp(my $indir = $_);
      my $crm = (split /\//,$indir)[-1];
      `mkdir -p $outdir/$crm`;
      warn "train IMM model on $crm\n";
      my $prepareModelAndData = "$Bin/CRM_vs_bkg/prepareModelAndData.pl";
      `perl $prepareModelAndData $Bin/../$indir $outdir $ktimes $nfolds`;
    }
    close IN;
    warn "Step1 is done!\n";
  }

  if ($step =~/2/)
  {
    ##=============================================
    # Step2: 
    # Generate msIMM score as features for each seq
    ##=============================================
    warn "Step2 ----------------------------------\n";
    warn "Generate msIMM score features...\n";
    for my $crm (@crmNames)
    {
     warn "$crm\n";
     for (my $k=1;$k<=$ktimes;$k++)
     {
       for (my $i=1;$i<=$nfolds;$i++)
       {
          warn "time $k fold $i\n";
          `perl $Bin/CRM_vs_bkg/generateIMMScoreFeature.pl $crm $outdir $k $i`;
       }
     }
    }
    # filter CRM training seq that are in the same group as test seq
    `perl $Bin/CRM_vs_bkg/filterGroupCRM.pl $datadir $outdir $crmGroupTable $ktimes $nfolds`;
    warn "Step2 is done!\n";
  }


  if ($step =~/3/)
  {
    ##=============== 
    # Step3.1:
    # msIMM baseline
    ##===============
    warn "Step3 ----------------------------------\n";
    warn "Step3.1 msIMM baseline...\n";
    for my $crm (@crmNames){
      warn "$crm\n";
      `perl $Bin/CRM_vs_bkg/msIMMBaseline.pl $crm $outdir $ktimes $nfolds`;
    }
    warn "Step3.1 is done!\n";

    ##===================== 
    # Step3.2:
    # Summarize msIMM AUCs
    ##=====================
    warn "Step3.2 Summarize msIMM AUCs...\n";
    combine_auc($outdir, \@crmNames, "summaryAUC_msIMMBaseline.txt", "IMM.average.auc");
    warn "Step3.2 is done!\n";
  }

  if ($step =~/4/)
  {
    ##==================
    # Step4.1:
    # IMM-SVM prediction 
    ##==================
    warn "Step4 ----------------------------------\n";
    warn "Step4.1 IMM-SVM prediction...\n";
    for my $crm (@crmNames)
    {
      warn "$crm\n";
      `perl $Bin/CRM_vs_bkg/msIMM_SVM.pl $crm $outdir $ktimes $nfolds`;
    }
    warn "Step4.1 is done!\n";


    ##============================= 
    # Step4.2:
    # Summarize IMM-SVM prediction 
    ##=============================
    warn "Step4.2 Summarize IMM-SVM AUCs...\n";
    combine_auc($outdir, \@crmNames, "summaryAUC_IMM_SVM.txt", "IMM_SVM.average.auc");
    warn "Step4.2 is done!\n";
  }

  if ($step =~/5/)
  {
    ##=================== 
    # Step5.1:
    # IMM-RF prediction 
    ##===================
    warn "Step5 ----------------------------------\n";
    warn "Step5.1 IMM-RF prediction...\n";
    for my $crm (@crmNames){
      warn "$crm\n";
      `perl $Bin/CRM_vs_bkg/msIMM_RF.pl $crm $outdir $ktimes $nfolds`;
    }
    warn "Step5.1 is done!\n";

    ##====================== 
    # Step5.2:
    # Summarize IMM-RF AUCs 
    ##======================
    warn "Step5.2 Summarize IMM-RF AUCs...\n";
    combine_auc($outdir, \@crmNames, "summaryAUC_IMM_RF.txt", "IMM_RF.average.auc");
    warn "Step5.2 is done!\n";
  }


  if ($step =~/6/)
  {
    ##===================================
    # Step6.1:
    # generate kmer features for kmerSVM 
    ##===================================
    warn "Step6 ----------------------------------\n";
    warn "Step6.1 generate kmer features for kmerSVM...\n";
    for my $crm (@crmNames)
    {
      warn "$crm\n";
      `perl $Bin/CRM_vs_bkg/kmerSVMFeature.pl $crm $outdir $ktimes $nfolds`;
    }
    warn "Step6.1 is done!\n";
  }


  if ($step =~/7/)
  {
    ##====================
    # Step7.1:
    # kmerSVM prediction 
    ##====================
    warn "Step7 ----------------------------------\n";
    warn "Step7.1 kmerSVM prediction...\n";
    for my $crm (@crmNames)
    {
      `perl $Bin/CRM_vs_bkg/kmerSVM.pl $crm $outdir $ktimes $nfolds`;
    }
    warn "Step7.1 is done!\n";


    ##===================== 
    # Step7.2:
    # Summarize kmerSVM AUC 
    ##=====================
    warn "Step7.2 Summarize kmerSVM AUC...\n";
    combine_auc($outdir, \@crmNames, "summaryAUC_kmerSVM.txt", "kmerSVM.average.auc");
    warn "Step7.2 is done!\n";
  }

  if ($step =~/8/)
  {
    ##===================== 
    # Step8.1:
    # Ensemble prediction 
    ##=====================
    warn "Step8 ----------------------------------\n";
    warn "Step8.1 Ensemble prediction...\n";
    for my $crm (@crmNames)
    {
      warn "$crm\n";
      `perl $Bin/CRM_vs_bkg/ensembleModel.pl $crm $outdir $ktimes $nfolds`;
    }
    warn "Step8.1 is done!\n";


    ##===========================
    # Step8.2:
    # summarize ensembleModel AUC 
    ##===========================
    warn "Step8.2 summarize ensembleModel AUC...\n";
    combine_auc($outdir, \@crmNames, "summaryAUC_ensembleModel.txt", "ensembleModel.average.auc");
    warn "Step8.2 is done!\n";
  }

  warn "All steps are doen!\n";
}

#================================
# helper function to combine AUCs 
# and output to a file
#================================
sub combine_auc {
  my $outdir = shift(@_);
  my $crmNames = shift(@_);
  my $out_file = shift(@_);
  my $auc_file = shift(@_);

  open OUT,">$outdir/$out_file" or die "cannot open $outdir/$out_file";
  for my $crm (@{$crmNames})
  {
     my $auc_path = "$outdir/$crm/$auc_file";
     if (-e $auc_path)
     {
         print OUT "$crm\t";
         open IN,$auc_path or die "cannot open $auc_path\n";
         while (<IN>)
         {
             chomp(my $second = $_);
             my $aucValue = (split /\s+/,$second)[-1];
             print OUT "$aucValue\n";
         }
         close IN;
     }
     else{
         print OUT "$crm\t-\n";
     }
  }
  close OUT;
}

#=============
# main routine
#=============
sub IMMBoost {
  if ($task eq "crm_vs_bkg"){
    CRM_vs_bkg();
  }
  elsif ($task eq "crm_vs_crm"){
    CRM_vs_CRM();
  }
  else {
    die `pod2text $0`;
  }
}

# call IMMBoost()
IMMBoost();