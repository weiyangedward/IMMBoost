=head1 Description:
  
  This script wrap the pipeline for IMMBoost (
  CRM_vs_bkg)

=head1 Usage:

  perl IMMBoost.pl [options] CRMList Outdir Datadir crmGroupTable

  --task <str>      Two modes: 1) "--task crm_vs_bkg": classify CRM from background genomic seq; 2) "--task crm_vs_crm": classify CRM from other CRM seq

=head1 Example:

  perl IMMBoost.pl -task crm_vs_bkg CRMsetsList.txt ../sampleOutput/CRM_vs_bkg/ ../sampleData/CRMsets/ CRM.group.V3.txt

=cut

use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;

my ($task);

GetOptions(
  "task:s"=>\$task,
);

die `pod2text $0` if (@ARGV!=4);



my $crmLst = $ARGV[0]; # a list file contains path to the dir of CRM datasets
my $outdir = $ARGV[1]; # outdir
my $datadir = $ARGV[2];
my $crmGroupTable = $ARGV[3];

#===== CRM_vs_CRM =======
sub CRM_vs_CRM {
  #=========
  # Step1:
  # for each CRM set, train a msIMM model and 
  # generate 10trials x 5folds cross validation 
  # test and training data
  #========
  warn "Step1 prepare Model And Data...\n";
  my @crmNames = ();
  open IN,$crmLst;
  while (<IN>)
  {
      chomp(my $indir = $_);
      my $crm = (split /\//,$indir)[-1];
      push @crmNames,$crm;
     `mkdir $outdir/$crm` unless (-d "$outdir/$crm");
     warn "train IMM model on $crm\n";
     `perl $Bin/CRM_vs_CRM/prepareModelAndData.pl $indir $outdir`;
  }
  close IN;
  warn "Step1 is done!\n";

  ##======== 
  # Step2: 
  # Generate msIMM score as features for each seq
  ##========
  warn "Step2 Generate msIMM score features...\n";
  for my $crm (@crmNames)
  {
     warn "$crm\n";
     for (my $k=1;$k<=10;$k++)
     {
         for (my $i=1;$i<=5;$i++)
         {
             warn "trial$k fold$i\n";
             `perl $Bin/CRM_vs_CRM/generateIMMScoreFeature.pl $crm $outdir $k $i`;
         }
     }
  }
  warn "fil group...\n";
  `perl $Bin/CRM_vs_CRM/filterGroupCRM.pl $datadir $outdir $crmGroupTable`;
  warn "fil Dmel...\n";
  `perl $Bin/CRM_vs_CRM/filterGroupCRM.DmelTrainData.pl $datadir $outdir $crmGroupTable`;
  warn "Step2 is done!\n";

  ##======== 
  # Step3.1:
  # msIMM baseline
  ##========
  warn "Step3.1 msIMM baseline...\n";
  for my $crm (@crmNames){
    warn "$crm\n";
    `perl $Bin/CRM_vs_CRM/msIMMBaseline.pl $crm $outdir`;
  }
  warn "Step3.1 is done!\n";

  ##====== 
  # Step3.2:
  # Summarize msIMM AUCs
  ##======
  warn "Step3.2 Summarize msIMM AUCs...\n";
  open OUT,">$outdir/summaryAUC_msIMMBaseline.txt";
  for my $crm (@crmNames)
  {
         if (-e "$outdir/$crm/IMM.average.auc")
         {
             print OUT "$crm\t";
             open IN,"$outdir/$crm/IMM.average.auc";
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
  warn "Step3.2 is done!\n";

  ##======== 
  # Step4.1:
  # IMM-SVM prediction 
  ##========
  warn "Step4.1 IMM-SVM prediction...\n";
  for my $crm (@crmNames)
  {
    warn "$crm\n";
    `perl $Bin/CRM_vs_CRM/msIMM_SVM.pl $crm $outdir`;
  }
  warn "Step4.1 is done!\n";



  ##======== 
  # Step4.2:
  # Summarize IMM-SVM prediction 
  ##========
  warn "Step4.2 Summarize IMM-SVM AUCs...\n";
  open OUT,">$outdir/summaryAUC_IMM_SVM.txt";
  for my $crm (@crmNames)
  {
         if (-e "$outdir/$crm/IMM_SVM.average.auc")
         {
             print OUT "$crm\t";
             open IN,"$outdir/$crm/IMM_SVM.average.auc";
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
  warn "Step4.2 is done!\n";


  ##======== 
  # Step5.1:
  # IMM-RF prediction 
  ##========
  warn "Step5.1 IMM-RF prediction...\n";
  for my $crm (@crmNames){
    warn "$crm\n";
    `perl $Bin/CRM_vs_CRM/msIMM_RF.pl $crm $outdir`;
  }
  warn "Step5.1 is done!\n";

  ##======== 
  # Step5.2:
  # Summarize IMM-RF AUCs 
  ##========
  warn "Step5.2 Summarize IMM-RF AUCs...\n";
  open OUT,">$outdir/summaryAUC_IMM_RF.txt";
  for my $crm (@crmNames)
  {
         if (-e "$outdir/$crm/IMM_RF.average.auc")
         {
             print OUT "$crm\t";
             open IN,"$outdir/$crm/IMM_RF.average.auc";
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
  warn "Step5.2 is done!\n";



  ##======= 
  # Step6.1:
  # generate kmer features for kmerSVM 
  ##=======
  warn "Step6.1 generate kmer features for kmerSVM...\n";
  for my $crm (@crmNames)
  {
    warn "$crm\n";
    `perl $Bin/CRM_vs_CRM/kmerSVMFeature.pl $crm $outdir`;
  }
  warn "Step6.1 is done!\n";


  ##======== 
  # Step6.2:
  # kmerSVM prediction 
  ##========
  warn "Step6.2 kmerSVM prediction...\n";
  for my $crm (@crmNames)
  {
    warn "$crm\n";
    `perl $Bin/CRM_vs_CRM/kmerSVM.pl $crm $outdir`;
  }
  warn "Step6.2 is done!\n";

  ##======= 
  # Step6.3:
  # Summarize kmerSVM AUC 
  ##=======
  warn "Step6.3 Summarize kmerSVM AUC...\n";
  open OUT,">$outdir/summaryAUC_kmerSVM.txt";
  for my $crm (@crmNames)
  {
     if (-e "$outdir/$crm/kmerSVM.average.auc")
     {
         print OUT "$crm\t";
         open IN,"$outdir/$crm/kmerSVM.average.auc";
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
  warn "Step6.3 is done!\n";


  ##=========== 
  # Step7.1:
  # Ensemble prediction 
  ##===========
  warn "Step7.1 Ensemble prediction...\n";
  for my $crm (@crmNames)
  {
    warn "$crm\n";
    `perl $Bin/CRM_vs_CRM/ensembleModel.pl $crm $outdir`;
  }
  warn "Step7.1 is done!\n";


  ##==========
  # Step7.2:
  # summarize ensembleModel AUC 
  ##==========
  warn "Step7.1 summarize ensembleModel AUC...\n";
  open OUT,">$outdir/summaryAUC_ensembleModel.txt";
  for my $crm (@crmNames)
  {
     if (-e "$outdir/$crm/ensembleModel.average.auc")
     {
         print OUT "$crm\t";
         open IN,"$outdir/$crm/ensembleModel.average.auc" or die "cannot open $outdir/$crm/ensembleModel.average.auc\n";
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
  warn "Step7.2 is done!\n";

}

#======= CRM_vs_bkg =========
sub CRM_vs_bkg {

  #=========
  # Step1:
  # for each CRM set, train a msIMM model and 
  # generate 10trials x 5folds cross validation 
  # test and training data
  #========
  warn "Step1 prepare Model And Data...\n";
  my @crmNames = (); # an array to store CRMnames
  open IN,$crmLst;
  while (<IN>)
  {
    chomp(my $indir = $_);
    my $crm = (split /\//,$indir)[-1];
    push @crmNames,$crm;
    `mkdir $outdir/$crm` unless (-d "$outdir/$crm");
    warn "train IMM model on $crm\n";
    my $prepareModelAndData = "$Bin/CRM_vs_bkg/prepareModelAndData.pl";
    `perl $prepareModelAndData $indir $outdir`;
  }
  close IN;
  warn "Step1 is done!\n";

  ##======== 
  # Step2: 
  # Generate msIMM score as features for each seq
  ##========
  warn "Step2 Generate msIMM score features...\n";
  for my $crm (@crmNames)
  {
   warn "$crm\n";
   for (my $k=1;$k<=10;$k++)
   {
     for (my $i=1;$i<=5;$i++)
     {
        warn "trial $k fold $i\n";
        `perl $Bin/CRM_vs_bkg/generateIMMScoreFeature.pl $crm $outdir $k $i`;
     }
   }
  }
  # filter CRM training seq that are in the same group as test seq
  `perl $Bin/CRM_vs_bkg/filterGroupCRM.pl $datadir $outdir $crmGroupTable`;
  warn "Step2 is done!\n";

  ##======== 
  # Step3.1:
  # msIMM baseline
  ##========
  warn "Step3.1 msIMM baseline...\n";
  for my $crm (@crmNames){
    warn "$crm\n";
    `perl $Bin/CRM_vs_bkg/msIMMBaseline.pl $crm $outdir`;
  }
  warn "Step3.1 is done!\n";

  ##====== 
  # Step3.2:
  # Summarize msIMM AUCs
  ##======
  warn "Step3.2 Summarize msIMM AUCs...\n";
  open OUT,">$outdir/summaryAUC_msIMMBaseline.txt";
  for my $crm (@crmNames)
  {
    if (-e "$outdir/$crm/IMM.average.auc")
    {
        print OUT "$crm\t";
        open IN,"$outdir/$crm/IMM.average.auc";
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
  warn "Step3.2 is done!\n";


  ##======== 
  # Step4.1:
  # IMM-SVM prediction 
  ##========
  warn "Step4.1 IMM-SVM prediction...\n";
  for my $crm (@crmNames)
  {
    warn "$crm\n";
    `perl $Bin/CRM_vs_bkg/msIMM_SVM.pl $crm $outdir`;
  }
  warn "Step4.1 is done!\n";


  ##======== 
  # Step4.2:
  # Summarize IMM-SVM prediction 
  ##========
  warn "Step4.2 Summarize IMM-SVM AUCs...\n";
  open OUT,">$outdir/summaryAUC_IMM_SVM.txt";
  for my $crm (@crmNames)
  {
         if (-e "$outdir/$crm/IMM_SVM.average.auc")
         {
             print OUT "$crm\t";
             open IN,"$outdir/$crm/IMM_SVM.average.auc";
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
  warn "Step4.2 is done!\n";

  ##======== 
  # Step5.1:
  # IMM-RF prediction 
  ##========
  warn "Step5.1 IMM-RF prediction...\n";
  for my $crm (@crmNames){
    warn "$crm\n";
    `perl $Bin/CRM_vs_bkg/msIMM_RF.pl $crm $outdir`;
  }
  warn "Step5.1 is done!\n";

  ##======== 
  # Step5.2:
  # Summarize IMM-RF AUCs 
  ##========
  warn "Step5.2 Summarize IMM-RF AUCs...\n";
  open OUT,">$outdir/summaryAUC_IMM_RF.txt";
  for my $crm (@crmNames)
  {
         if (-e "$outdir/$crm/IMM_RF.average.auc")
         {
             print OUT "$crm\t";
             open IN,"$outdir/$crm/IMM_RF.average.auc";
             while (<IN>){
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
  warn "Step5.2 is done!\n";


  ##======= 
  # Step6.1:
  # generate kmer features for kmerSVM 
  ##=======
  warn "Step6.1 generate kmer features for kmerSVM...\n";
  for my $crm (@crmNames)
  {
    warn "$crm\n";
    `perl $Bin/CRM_vs_bkg/kmerSVMFeature.pl $crm $outdir`;
  }
  warn "Step6.1 is done!\n";


  ##======== 
  # Step6.2:
  # kmerSVM prediction 
  ##========
  warn "Step6.2 kmerSVM prediction...\n";
  for my $crm (@crmNames)
  {
    `perl $Bin/CRM_vs_bkg/kmerSVM.pl $crm $outdir`;
  }
  warn "Step6.2 is done!\n";


  ##======= 
  # Step6.3:
  # Summarize kmerSVM AUC 
  ##=======
  warn "Step6.3 Summarize kmerSVM AUC...\n";
  open OUT,">$outdir/summaryAUC_kmerSVM.txt";
  for my $crm (@crmNames)
  {
     if (-e "$outdir/$crm/kmerSVM.average.auc")
     {
         print OUT "$crm\t";
         open IN,"$outdir/$crm/kmerSVM.average.auc";
         while (<IN>){
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
  warn "Step6.3 is done!\n";

  ##=========== 
  # Step7.1:
  # Ensemble prediction 
  ##===========
  warn "Step7.1 Ensemble prediction...\n";
  for my $crm (@crmNames)
  {
    warn "$crm\n";
    `perl $Bin/CRM_vs_bkg/ensembleModel.pl $crm $outdir`;
  }
  warn "Step7.1 is done!\n";


  ##==========
  # Step7.2:
  # summarize ensembleModel AUC 
  ##==========
  warn "Step7.1 summarize ensembleModel AUC...\n";
  open OUT,">$outdir/summaryAUC_ensembleModel.txt";
  for my $crm (@crmNames)
  {
     if (-e "$outdir/$crm/ensembleModel.average.auc")
     {
         print OUT "$crm\t";
         open IN,"$outdir/$crm/ensembleModel.average.auc" or die "cannot open $outdir/$crm/ensembleModel.average.auc\n";
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
  warn "Step7.2 is done!\n";
}

sub IMMBoost {
  if ($task == "crm_vs_bkg"){
    CRM_vs_bkg();
  }
  elsif ($task == "crm_vs_crm"){
    CRM_vs_CRM();
  }
  else {
    die `pod2text $0`;
  }
}

# call IMMBoost()
IMMBoost();