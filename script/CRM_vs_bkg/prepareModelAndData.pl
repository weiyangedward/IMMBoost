=head1 Description

    This script creates cross validation sets for
    classification test with CRM seq as positive data
    and background genomic seq as negative data
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

# die "Usage:\n\tperl $0 CRMdir Outdir\n" unless @ARGV==2;

die `pod2text $0` if (@ARGV!=4);

my $crmDir = $ARGV[0]; # CRMdir
my $outdir = $ARGV[1];
my $times = $ARGV[2];
my $nfolds = $ARGV[3];

# relative paths to IMM
my $train = "$Bin/../../src/imm/bin/imm_build";
my $pred = "$Bin/../../src/imm/bin/imm_score";
# my $auc = "$Bin/auc.R";

# data
my $crmFa = "$crmDir/fasta/CRM.fasta"; # CRM seq
my $negFa = "$crmDir/fasta/randomGenomicSeq.fasta"; # background seq
my $msCRM = "$crmDir/fasta/msCRM.fasta"; # msCRM seq
# my $msNeg = "$crmDir/fasta/neg.acc.lenGC.100.fasta";
my $msNeg = $negFa;
my $crm = (split /\//,$crmDir)[-1];
if (-d "$outdir/$crm")
{
    `rm -rf $outdir/$crm`;
}
`mkdir -p $outdir/$crm`;
# global variables
my (%crmSeq, @crmSeqID, %negSeq, @negSeqID, 
    $crmNum, $negNum, @crmSeqIDrand, @negSeqIDrand);

##==============
# read in CRMs seq
##==============
my %crmGeneralID = ();  ## msCRM general ID, e.g., bl_blMir7_OE
my %msCRMSeq = ();      ## msCRM ID
sub parseCRM {
    # store CRM seq
    my $fa1 = Bio::SeqIO->new(-file=>$crmFa,-format=>'Fasta');
    while (my $curSeq = $fa1->next_seq())
    {
        my $id = $curSeq->id();
        my $seq = $curSeq->seq();
        $crmSeq{$id} = $seq; # hash CRM seq
        push @crmSeqID,$id; # store CRM id
        $crmNum ++;
    }
    #  store msCRM seq
    my $fa2 = Bio::SeqIO->new(-file=>$msCRM,-format=>'Fasta');
    while (my $curSeq = $fa2->next_seq())
    {
        my $id = $curSeq->id(); # e.g., Dana_D_5
        my $seq = $curSeq->seq();
        $msCRMSeq{$id} = $seq; # hash msCRM seq
        my $generalID = $2 if $id =~ /(\w+?)_(\S+)/; # generalID as in CRM seq, e.g., D_5
        $crmGeneralID{$generalID} = 1; # hash generalID, this will be used to get msCRM seq as training data given CRM id
    }
}
##=====================
# read in negative data 
##=====================
sub parseNeg 
{
    my $fa1 = Bio::SeqIO->new(-file=>$negFa,-format=>'Fasta');
    while (my $curSeq = $fa1->next_seq())
    {
        my $id = $curSeq->id();
        my $seq = $curSeq->seq();
        $negSeq{$id} = $seq; # hash neg seq
        push @negSeqID,$id;
        $negNum ++;
    }
}
##======= 
# Create dir 
##========
sub createDir{
    # 10 trials
    for (my $k=1; $k<=$times; $k++)
    {
        `mkdir -p $outdir/$crm/time$k/`;
        # 5 fold
        for (my $i = 1; $i<=5; $i++)
        {
           `mkdir -p $outdir/$crm/time$k/fold$i`;
        }
    }
}
##==========
# Create 10 trials x 5-fold cross validation sets
# by randomly partition the CRM and neg sets
##==========
sub sepData {
    # 10 trials
    for (my $k=1;$k<=$times;$k++)
    {
        my ($first, $second) = shufData(); # get shufl CRM and neg id
        my @crmSeqIDrand = @$first;
        my @negSeqIDrand = @$second;

        my $crmStart = 0;
        my $crmEnd = 0;
        my $negStart = 0;
        my $negEnd = 0;

        # 5 fold
        for (my $i = 1; $i <= $nfolds; $i ++)
        {
            # use start and end index to get testing and training set from CRM and neg id array
            $crmEnd = floor(($i/$nfolds) * $crmNum) - 1;
            $negEnd = floor(($i/$nfolds) * $negNum) - 1;
            # output positive data
            open OUT1,">$outdir/$crm/time$k/fold$i/test.crm.fasta";
            open OUT2,">$outdir/$crm/time$k/fold$i/train.crm.fasta";
            open OUT3,">$outdir/$crm/time$k/fold$i/test.label"; # label of both of pos and neg test data
            open OUT4,">$outdir/$crm/time$k/fold$i/train.label"; # label of both of pos and neg training data

            # output testing CRM seq
            my %trainCRM = ();
            for (my $j=0; $j<=$crmNum-1;$j++)
            {
                my $newID = $crmSeqIDrand[$j];
                # put CRM seq that within cur range into testing set
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
            # output msCRM training seq
            for my $id (keys %msCRMSeq)
            {
                my $generalID = $2 if $id =~ /(\w+?)_(\S+)/;
                if ($trainCRM{$generalID})
                {
                    print OUT2 ">$id\n$msCRMSeq{$id}\n";
                    print OUT4 "$id\t1\n";
                }
            }
            close OUT1;
            close OUT2;
            # output neg data
            open OUT1,">$outdir/$crm/time$k/fold$i/test.neg.fasta";
            open OUT2,">$outdir/$crm/time$k/fold$i/train.neg.fasta";
            
            for (my $j=0;$j<=$negNum-1;$j++)
            {
                my $newID = $negSeqIDrand[$j];
                # output neg testing seq if within cur range
                if ($j <= $negEnd && $j >= $negStart)
                {
                    print OUT1 ">$negSeqIDrand[$j]\n$negSeq{$negSeqIDrand[$j]}\n";
                    print OUT3 "$newID\t-1\n";
                }
                # mark other as training seq
                else
                {
                    print OUT2 ">$newID\n$negSeq{$negSeqIDrand[$j]}\n";
                    print OUT4 "$newID\t-1\n";
                }
            }
            close OUT1;
            close OUT2;
            close OUT3;
            close OUT4;
            # update start index of CRM and neg id array
            $crmStart = $crmEnd + 1;
            $negStart = $negEnd + 1;
            # combine pos and neg test data
            `cat $outdir/$crm/time$k/fold$i/test.crm.fasta $outdir/$crm/time$k/fold$i/test.neg.fasta > $outdir/$crm/time$k/fold$i/test.crm.and.neg.fasta`;
            # combine pos and neg training data
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
    # use IMM to train both of pos and neg models
    `$train -r -k 6 < $msNeg > $outdir/$crm/allData/neg.model`;
    `$train -r -k 6 < $msCRM > $outdir/$crm/allData/crm.model`;
    # copy msCRM and neg seq files
    `cp $msCRM $outdir/$crm/allData/msCRM.fa`;
    `cp $msNeg $outdir/$crm/allData/msNeg.fa`;
}

##=========================== 
# shuffle CRM and neg set id 
##===========================
sub shufData {
    my @crmSeqIDrand = shuffle @crmSeqID; # shufl CRM id
    my @negSeqIDrand = shuffle @negSeqID; # shufl neg id
    return (\@crmSeqIDrand, \@negSeqIDrand);
}

# call functions
sub run {
    parseCRM(); # store CRM seq
    parseNeg();
    createDir();
    sepData();
    trainAllData();
}

run();
