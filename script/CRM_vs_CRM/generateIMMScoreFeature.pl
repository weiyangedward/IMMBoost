=head1 Description

    This script is used to generate msIMM score as 
    features for each seq in test and training sets. 
    These new features will be used to train IMMBoost 
    models to classify CRM seq. Note that BioPerl is 
    needed. Negative seq will be from other CRMs

=head1 Usage

  perl generateIMMScoreFeature.pl CRMname OutDir timeNum foldNum

=cut

use strict;
use warnings;
# use lib '/shared-mounts/sinhas/lib/';
use File::Basename;
use Bio::SeqIO;
use FindBin qw($Bin);

# die "Usage: perl $0 CRMname OutDir timeNum foldNum\n" unless @ARGV==4;
die `pod2text $0` if (@ARGV!=4);

# my $tmpDir = $ENV{TMPDIR};
# my $tmpModelDir = "$tmpDir/model";
# `mkdir $tmpModelDir` unless (-e "$tmpModelDir");

my $train = "$Bin/../../src/imm/bin/imm_build";
my $pred = "$Bin/../../src/imm/bin/imm_score";

my $crm = $ARGV[0];
my $modelDir = $ARGV[1];
my $timeNum = $ARGV[2];
my $foldNum = $ARGV[3];

# get CRMnames
my @crmNames;
sub getNames 
{
    while (<$modelDir/*>)
    {
        chomp(my $crmDir = $_);
        my $crmName = basename($crmDir);
        push @crmNames,$crmName;
    }
}

# generate msIMM features for each seq in test and training sets
sub pred {
        my $crmDir = "$modelDir/$crm";
        my $k = $timeNum;
        my $i = $foldNum;
        my $homeDir = "$crmDir/time$k/fold$i";  
        my $trainFa = "$crmDir/time$k/fold$i/train.crm.and.neg.fasta";
        my $testFa = "$crmDir/time$k/fold$i/test.crm.and.neg.fasta";
            
        ##==== hash train data =======##
        my %trainSeq = ();
        my $fa1 = Bio::SeqIO->new(-file=>$trainFa,-format=>'Fasta');
        while (my $curSeq = $fa1->next_seq()){
            my $id = $curSeq->id();
            my $seq = $curSeq->seq();
            $trainSeq{$id} = $seq;
        }      
        ##====== hash test data ======
        my %testSeq = ();
        my $fa2 = Bio::SeqIO->new(-file=>$testFa,-format=>'Fasta');
        while (my $curSeq = $fa2->next_seq()){
            my $id = $curSeq->id();
            my $seq = $curSeq->seq();
            $testSeq{$id} = $seq;
        }
        ##==== Build IMMs to generate IMM score features for training and test sequences 
        ##====

        # store all of the IMM output
        my %output = ();
        for my $otherCRM (@crmNames)
        {             
            my $crmFa = "$modelDir/$otherCRM/allData/msCRM.fasta";
            my $negFa = "$modelDir/$otherCRM/allData/msNeg.fasta";
            ##===== hash msCRM data ======
            my %crmSeq = ();
            my %crmGeneralID = ();
            my $fa3 = Bio::SeqIO->new(-file=>$crmFa,-format=>'Fasta');
            while (my $curSeq = $fa3->next_seq()){
                my $id = $curSeq->id();
                my $seq = $curSeq->seq();
                $crmSeq{$id} = $seq;
            }

            ##===== read msNeg data =====
            my %negSeq = ();
            my %negGeneralID = ();
            my $fa4 = Bio::SeqIO->new(-file=>$negFa,-format=>'Fasta');
            while (my $curSeq = $fa4->next_seq()){
                my $id = $curSeq->id();
                my $seq = $curSeq->seq();
                $negSeq{$id} = $seq;
            }

            ##==== for each training seq, generate pos and neg training seq to train a msIMM, and then score this seq
            ##====
            for my $id (sort keys %trainSeq)
            {
                ##==== Create test fasta files  ====
                open OUT,">$homeDir/$id.fa";                   
                print OUT ">$id\n$trainSeq{$id}\n";
                close OUT;
                ##==== read general ID ====
                my $thisSeqGeneralID = $id;
                if ($id =~ /(\w+?)_(\S+)/)
                {
                    $thisSeqGeneralID = $2;
                }
                #==== create CRM training seq for msIMM, where seq that are the same as target seq or in test data are filtered since target seq will be used to train a model to pred on test data
                #=====
                open OUT,">$homeDir/$otherCRM.crms.fa";
                for my $crmID (sort keys %crmSeq)
                {
                    my $generalID = $2 if $crmID =~ /(\w+?)_(\S+)/;
                    # add sequence to CRM training set if it is not the same as target sequences or in test set
                    if ($generalID ne $thisSeqGeneralID && !$testSeq{$generalID})
                    {
                        print OUT ">$crmID\n$crmSeq{$crmID}\n";
                    }
                }
                close OUT;

                #=== Create neg training seq for msIMM, where seq that are the same as target seq or in test data are filtered
                #====
                open OUT,">$homeDir/$otherCRM.neg.fa";
                for my $negID (sort keys %negSeq)
                {
                    my $generalID = $2 if $negID =~ /(\w+?)_(\S+)/;
                    ## add sequence to CRM training set if it is not the same as target sequences or in test set
                    if ($generalID ne $thisSeqGeneralID && !$testSeq{$generalID})
                    {
                        print OUT ">$negID\n$negSeq{$negID}\n";
                    }
                }
                close OUT;

                ##==== Train msIMM and score target seq ====
                `$train -r -k 6 < $homeDir/$otherCRM.neg.fa > $homeDir/$otherCRM.neg.model`;
                `$train -r -k 6 < $homeDir/$otherCRM.crms.fa > $homeDir/$otherCRM.crms.model`;
                my $crmModel = "$homeDir/$otherCRM.crms.model";
                my $negModel = "$homeDir/$otherCRM.neg.model";

                ## hash pred score
                $output{"train"}{$id}{$otherCRM} = `$pred -f -n $negModel $crmModel $homeDir/$id.fa`;

                # delete intermediate files
                `rm $homeDir/$otherCRM.neg.fa`;
                `rm $homeDir/$otherCRM.crms.fa`;
                `rm $homeDir/$id.fa`;
                `rm $homeDir/$otherCRM.crms.model`;
                `rm $homeDir/$otherCRM.neg.model`; 
                ##==== delete everything in tmp directory
                # `rm -rf $homeDir/*`;
            }
            
            ##==== for each test seq, generate pos and neg training seq to train a msIMM, and then score this seq
            ##====   
            for my $id (sort keys %testSeq)
            {
                ##= Create test fasta files  
                open OUT,">$homeDir/$id.fa";
                print OUT ">$id\n$testSeq{$id}\n";
                close OUT;
                ##== Create training CRMs for IMM ==
                open OUT,">$homeDir/$otherCRM.crms.fa";
                for my $crmID (sort keys %crmSeq)
                {
                    my $generalID = $2 if $crmID =~ /(\w+?)_(\S+)/;
                    if ($generalID ne $id && !$testSeq{$generalID})
                    {
                        print OUT ">$crmID\n$crmSeq{$crmID}\n";
                    }
                }
                close OUT;
                ##====== Create training neg data for IMM ==========##
                open OUT,">$homeDir/$otherCRM.neg.fa";
                for my $negID (sort keys %negSeq)
                {
                    my $generalID = $2 if $negID =~ /(\w+?)_(\S+)/;
                    if ($generalID ne $id && !$testSeq{$generalID})
                    {
                        print OUT ">$negID\n$negSeq{$negID}\n";
                    }
                }
                close OUT;
                ##==== Train msIMM and score target seq ====
                `$train -r -k 6 < $homeDir/$otherCRM.neg.fa > $homeDir/$otherCRM.neg.model`;
                `$train -r -k 6 < $homeDir/$otherCRM.crms.fa > $homeDir/$otherCRM.crms.model`;
                my $crmModel = "$homeDir/$otherCRM.crms.model";
                my $negModel = "$homeDir/$otherCRM.neg.model";
                ##==== Store prediction into a hash table ====##
                $output{"test"}{$id}{$otherCRM} = `$pred -f -n $negModel $crmModel $homeDir/$id.fa`;

                # delete intermediate files
                `rm $homeDir/$otherCRM.neg.fa`;
                `rm $homeDir/$otherCRM.crms.fa`;
                `rm $homeDir/$id.fa`;
                `rm $homeDir/$otherCRM.crms.model`;
                `rm $homeDir/$otherCRM.neg.model`;
                ##==== delete everything in tmp directory
                # `rm -rf $homeDir/*`;
            }
        }
        ##====== Output IMM features for training sequences ======##
        open OUT1,">$homeDir/train.ensembFeat";
        print OUT1 "crm\t";
        for my $crm (sort @crmNames){
            print OUT1 "$crm\t";
        }
        print OUT1 "\n";

        for my $id (sort keys %{$output{"train"}})
        {
            print OUT1 "$id\t";
            for my $crm (sort keys %{$output{"train"}{$id}})
            {
                my @score = split /\s+/,$output{"train"}{$id}{$crm};
                print OUT1 "$score[1]\t";
            }
            print OUT1 "\n";
        }
        close OUT1;

        ##======= output test seq msIMM features to file ====
        open OUT2,">$homeDir/test.ensembFeat";
        print OUT2 "crm\t";
        for my $crm (sort @crmNames){
            print OUT2 "$crm\t";
        }
        print OUT2 "\n";
        for my $id (sort keys %{$output{"test"}})
        {
            print OUT2 "$id\t";
            for my $crm (sort keys %{$output{"test"}{$id}})
            {
                my @score = split /\s+/,$output{"test"}{$id}{$crm};
                print OUT2 "$score[1]\t";
            }
            print OUT2 "\n";
        }
        close OUT2;
        ##====== Move feature files to local directories ======##
        # `mv $homeDir/test.ensembFeat $crmDir/time$k/fold$i/test.ensembFeat`;
        # `mv $homeDir/train.ensembFeat $crmDir/time$k/fold$i/train.ensembFeat`;
}

# call functions
getNames();
pred();

