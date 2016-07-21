=head1 Description
  
  This script creates cross validation sets for
    classification test with CRMs of one specific expression
    domain as positive data and other CRMs as negative data.
    Note that BioPerl is needed.

=head1 Usage

  perl prepareModelAndData.pl CRMdir Outdir

=cut

use strict;
use warnings;
use FindBin qw($Bin);
use List::Util qw(shuffle);
use Bio::SeqIO;
use POSIX;

die `pod2text $0` if (@ARGV!=3);

my $crmDir = $ARGV[0];
my $outdir = $ARGV[1];
my $times = $ARGV[2];

my $train = "$Bin/../../src/imm/bin/imm_build";
my $pred = "$Bin/../../src/imm/bin/imm_score";
my $auc = "$Bin/IMMauc.R";

my $crmFa = "$crmDir/fasta/CRM.fasta";
my $negFa = "$crmDir/fasta/negCRM.fasta";
my $msCRM = "$crmDir/fasta/msCRM.fasta";
my $msNeg = "$crmDir/fasta/negmsCRM.fasta";
my $crm = (split /\//,$crmDir)[-1];
`mkdir -p $outdir/$crm`;
my (%crmSeq,@crmSeqID,%negSeq,@negSeqID,
    $crmNum,$negNum,@crmSeqIDrand,@negSeqIDrand);

##==============
# read in CRMs seq
##==============
my %crmGeneralID = ();  ## msCRM general ID
my %msCRMSeq = ();      ## msCRM ID
sub parseCRM {
    my $fa = Bio::SeqIO->new(-file=>$crmFa,-format=>'Fasta');
    while (my $nextSeq = $fa->next_seq()){
        my $id = $nextSeq->id();
        my $seq = $nextSeq->seq();
        $crmSeq{$id} = $seq;
        push @crmSeqID,$id;
        $crmNum ++;
    }

    my $fa2 = Bio::SeqIO->new(-file=>$msCRM,-format=>'Fasta');
    while (my $curSeq = $fa2->next_seq()){
        my $id = $curSeq->id();
        my $seq = $curSeq->seq();
        $msCRMSeq{$id} = $seq;
        my $generalID = $2 if $id =~ /(\w+?)_(\S+)/;
        $crmGeneralID{$generalID} = 1;
    }
}
##=====================
# read in negative data 
##=====================
my %negGeneralID = ();
my %msNegSeq = ();
sub parseNeg {
    my $fa = Bio::SeqIO->new(-file=>$negFa,-format=>'Fasta');
    while (my $nextSeq = $fa->next_seq()){
        my $id = $nextSeq->id();
        my $seq = $nextSeq->seq();
        $negSeq{$id} = $seq;
        push @negSeqID,$id;
        $negNum ++;
    }

    my $fa1 = Bio::SeqIO->new(-file=>$msNeg,-format=>'Fasta');
    while (my $nextSeq = $fa1->next_seq()){
        my $id = $nextSeq->id();
        my $seq = $nextSeq->seq();
        my $generalID = $2 if $id =~ /(\w+?)_(\S+)/;
        $msNegSeq{$id} = $seq;
        $negGeneralID{$generalID} = 1;
    }
}
##==================
# Create directories 
##==================
sub createDir{
    for (my $k=1; $k<=$times; $k++)
    {
        `mkdir -p $outdir/$crm/time$k/`;
        for (my $i = 1; $i<=5; $i++)
        {
            `mkdir -p $outdir/$crm/time$k/fold$i`;
        }
    }
}
##===============================================
# Create 10 trials x 5-fold cross validation sets
# by randomly partition the CRM and neg sets
##===============================================
sub sepData {
    
    for (my $k=1;$k<=$times;$k++)
    {
        my ($first, $second) = shufData();
        my @crmSeqIDrand = @$first;
        my @negSeqIDrand = @$second;
        my $crmStart = 0;
        my $crmEnd = 0;
        my $negStart = 0;
        my $negEnd = 0;

        for (my $i = 1; $i <= 5; $i ++)
        {
            $crmEnd = floor(($i/5) * $crmNum) - 1;
            $negEnd = floor(($i/5) * $negNum) - 1;

            open OUT1,">$outdir/$crm/time$k/fold$i/test.crm.fasta";
            open OUT2,">$outdir/$crm/time$k/fold$i/train.crm.fasta";
            open OUT3,">$outdir/$crm/time$k/fold$i/test.label";
            open OUT4,">$outdir/$crm/time$k/fold$i/train.label";

            my %trainCRM = ();
            for (my $j=0; $j<=$crmNum-1;$j++){
                my $newID = $crmSeqIDrand[$j];
                ##===== output test CRM ======
                if ($j <= $crmEnd && $j >= $crmStart)
                {
                    print OUT1 ">$crmSeqIDrand[$j]\n$crmSeq{$crmSeqIDrand[$j]}\n";
                    print OUT3 "$newID\t1\n";
                }
                # mark other CRM seq as in training set
                else
                {   
                    $trainCRM{$crmSeqIDrand[$j]} = 1;
                }
            }
            close OUT1;
            ##=== output msCRM train data ===
            for my $id (keys %msCRMSeq)
            {
                my $generalID = $2 if $id =~ /(\w+?)_(\S+)/;
                if ($trainCRM{$generalID})
                {
                    print OUT2 ">$id\n$msCRMSeq{$id}\n";
                    print OUT4 "$id\t1\n";
                }
            }
            close OUT2;

            open OUT1,">$outdir/$crm/time$k/fold$i/test.neg.fasta";
            open OUT2,">$outdir/$crm/time$k/fold$i/train.neg.fasta";

            my %trainNeg = ();
            for (my $j=0;$j<=$negNum-1;$j++)
            {
                my $newID = $negSeqIDrand[$j];
                ##===== output neg test ======##
                if ($j <= $negEnd && $j >= $negStart){
                    print OUT1 ">$negSeqIDrand[$j]\n$negSeq{$negSeqIDrand[$j]}\n";
                    print OUT3 "$newID\t-1\n";
                }
                # mark other CRM seq as in training set
                else
                {  
                    $trainNeg{$negSeqIDrand[$j]} = 1;
                }
            }
            # output neg seq to file
            for my $id (keys %msNegSeq)
            {
                my $generalID = $2 if $id =~ /(\w+?)_(\S+)/;
                if ($trainNeg{$generalID})
                {
                    print OUT2 ">$id\n$msNegSeq{$id}\n";
                    print OUT4 "$id\t-1\n";
                }
            }
            close OUT1;
            close OUT2;
            close OUT3;
            close OUT4;

            $crmStart = $crmEnd + 1;
            $negStart = $negEnd + 1;
            `cat $outdir/$crm/time$k/fold$i/test.crm.fasta $outdir/$crm/time$k/fold$i/test.neg.fasta > $outdir/$crm/time$k/fold$i/test.crm.and.neg.fasta`;
            `cat $outdir/$crm/time$k/fold$i/train.crm.fasta $outdir/$crm/time$k/fold$i/train.neg.fasta > $outdir/$crm/time$k/fold$i/train.crm.and.neg.fasta`;
        }
    }
}

##==================================================
# train msIMM for each CRM set, these msIMMs will be 
# used to score test and training data, after which
# the pred score will be used as feature to train a 
# classification model
##==================================================
sub trainAllData {
    `mkdir -p $outdir/$crm/allData`;
    #==================
    # train msCRM model 
    #==================
    `$train -r -k 6 < $msNeg > $outdir/$crm/allData/neg.model`; 
    `$train -r -k 6 < $msCRM > $outdir/$crm/allData/crm.model`;
    `cp $msCRM $outdir/$crm/allData/msCRM.fasta`;
    `cp $msNeg $outdir/$crm/allData/msNeg.fasta`;   
}

##======= 
# shuffle CRM and neg set id 
##=======
sub shufData {
    my @crmSeqIDrand = shuffle @crmSeqID;
    my @negSeqIDrand = shuffle @negSeqID;
    return (\@crmSeqIDrand, \@negSeqIDrand);
}

sub run {
    parseCRM();
    parseNeg();
    createDir();
    sepData();
    trainAllData();
}

run();
