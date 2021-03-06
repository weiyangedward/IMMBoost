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
use File::Basename;
use Bio::SeqIO;
use FindBin qw($Bin);

die `pod2text $0` if (@ARGV!=4);

my $train = "$Bin/../../src/imm/bin/imm_build";
my $pred = "$Bin/../../src/imm/bin/imm_score";

my $crm = $ARGV[0];
my $modelDir = $ARGV[1];
my $timeNum = $ARGV[2];
my $foldNum = $ARGV[3];

my %train_seq_label = ();
my %test_seq_label = ();

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

#======================================
# generate msIMM features for each seq 
# in test and training sets
#======================================
sub pred {
        my $crmDir = "$modelDir/$crm";
        my $k = $timeNum;
        my $i = $foldNum;
        my $homeDir = "$crmDir/time$k/fold$i";  
        my $trainFa = "$crmDir/time$k/fold$i/train.crm.and.neg.fasta";
        my $testFa = "$crmDir/time$k/fold$i/test.crm.and.neg.fasta";
            
        ##=================================
        # read train data into hashTable
        # key=id, value=seq_fa 
        #==================================
        my %trainSeq = ();
        my $fa1 = Bio::SeqIO->new(-file=>$trainFa,-format=>'Fasta');
        while (my $curSeq = $fa1->next_seq()){
            my $id = $curSeq->id();
            my $seq = $curSeq->seq();
            $trainSeq{$id} = $seq;
        }      
        ##==============
        # read test data into hashTable
        # key=id, value=seq_fa
        #===============
        my %testSeq = ();
        my $fa2 = Bio::SeqIO->new(-file=>$testFa,-format=>'Fasta');
        while (my $curSeq = $fa2->next_seq()){
            my $id = $curSeq->id();
            my $seq = $curSeq->seq();
            $testSeq{$id} = $seq;
        }
        ##==========================================
        # Build IMMs and predict IMM scores as features
        # for training and test sequences 
        ##==========================================

        # store all of the IMM output
        my %output = ();
        for my $otherCRM (@crmNames)
        {             
            # warn "$otherCRM ...\n";
            my $crmFa = "$modelDir/$otherCRM/allData/msCRM.fasta";
            my $negFa = "$modelDir/$otherCRM/allData/msNeg.fasta";
            ##===============================
            # read msCRM data into hashTable
            # key=id, value=seq_fa
            #================================
            my %crmSeq = ();
            my %crmGeneralID = ();
            my $fa3 = Bio::SeqIO->new(-file=>$crmFa,-format=>'Fasta');
            while (my $curSeq = $fa3->next_seq()){
                my $id = $curSeq->id();
                my $seq = $curSeq->seq();
                $crmSeq{$id} = $seq;
            }

            ##============================
            # read Neg data into hashTable
            # key=id, value=seq_fa
            #=============================
            my %negSeq = ();
            my %negGeneralID = ();
            my $fa4 = Bio::SeqIO->new(-file=>$negFa,-format=>'Fasta');
            while (my $curSeq = $fa4->next_seq()){
                my $id = $curSeq->id();
                my $seq = $curSeq->seq();
                $negSeq{$id} = $seq;
            }

            ##===========================================
            # for each training seq S:
            # 1. generate pos and neg training seq
            # 2. then use training to train IMM and
            # score on S
            ##===========================================
            for my $id (sort keys %trainSeq)
            {
                ##==== Create test fasta files  ====
                # warn "train seq: $id\n";
                open OUT,">$homeDir/$id.fa";
                print OUT ">$id\n$trainSeq{$id}\n";
                close OUT;
                ##==== read in general seq ID ====
                my $thisSeqGeneralID = $id;
                if ($id =~ /(\w+?)_(\S+)/)
                {
                    $thisSeqGeneralID = $2;
                }
                #=========================================
                # create CRM training seq for msIMM, 
                # where seq that are the same as 
                # target seq or in test data are filtered 
                # so that no cheating when using target seq 
                # to train a model and pred on test data
                #=========================================
                open OUT,">$homeDir/$otherCRM.crms.fa";
                for my $crmID (sort keys %crmSeq)
                {
                    my $generalID = $2 if $crmID =~ /(\w+?)_(\S+)/;
                    #==================================
                    # add sequence to CRM training set 
                    # if it is not the same as target 
                    # sequences or in test set
                    #==================================
                    if ($generalID ne $thisSeqGeneralID && !$testSeq{$generalID})
                    {
                        print OUT ">$crmID\n$crmSeq{$crmID}\n";
                    }
                }
                close OUT;

                #========================================
                # Create neg training seq for msIMM, 
                # where seq that are the same as 
                # target seq or in test data are filtered
                #========================================
                open OUT,">$homeDir/$otherCRM.neg.fa";
                for my $negID (sort keys %negSeq)
                {
                    my $generalID = $2 if $negID =~ /(\w+?)_(\S+)/;
                    ##================================= 
                    # add sequence to CRM training set 
                    # if it is not the same as target 
                    # sequences or in test set
                    #==================================
                    if ($generalID ne $thisSeqGeneralID && !$testSeq{$generalID})
                    {
                        print OUT ">$negID\n$negSeq{$negID}\n";
                    }
                }
                close OUT;

                ##================================
                # Train msIMM and score target seq 
                #=================================
                `$train -r -k 6 < $homeDir/$otherCRM.neg.fa > $homeDir/$otherCRM.neg.model`;
                `$train -r -k 6 < $homeDir/$otherCRM.crms.fa > $homeDir/$otherCRM.crms.model`;
                my $crmModel = "$homeDir/$otherCRM.crms.model";
                my $negModel = "$homeDir/$otherCRM.neg.model";

                ## add pred score to hashTable
                $output{"train"}{$id}{$otherCRM} = `$pred -f -n $negModel $crmModel $homeDir/$id.fa`;

                # delete intermediate files
                `rm $homeDir/$otherCRM.neg.fa`;
                `rm $homeDir/$otherCRM.crms.fa`;
                `rm $homeDir/$id.fa`;
                `rm $homeDir/$otherCRM.crms.model`;
                `rm $homeDir/$otherCRM.neg.model`; 
            }
            
            ##====================================
            # for each test seq S: 
            # 1. generate pos and neg training seq
            # 2. then use training to train IMM and
            # score on S
            ##=====================================
            for my $id (sort keys %testSeq)
            {
                ## Create test fasta files  
                open OUT,">$homeDir/$id.fa";
                # warn "test seq: $id\n";
                print OUT ">$id\n$testSeq{$id}\n";
                close OUT;
                ## Create training CRMs for IMM ==
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
                ##================================
                # Create training neg data for IMM 
                #=================================
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
                ##==== add pred score to hash table ====
                $output{"test"}{$id}{$otherCRM} = `$pred -f -n $negModel $crmModel $homeDir/$id.fa`;

                # delete intermediate files
                `rm $homeDir/$otherCRM.neg.fa`;
                `rm $homeDir/$otherCRM.crms.fa`;
                `rm $homeDir/$id.fa`;
                `rm $homeDir/$otherCRM.crms.model`;
                `rm $homeDir/$otherCRM.neg.model`;
            }
        }

        # read in train and test seq labels
        read_labels();
        ##==========================================
        # Output IMM features for training sequences 
        #===========================================
        open OUT1,">$homeDir/train.ensembFeat";
        open TRAIN_lib,">$homeDir/train.ensembFeat.lib";
        open TRAIN_lib_DMEL,">$homeDir/train.ensembFeat.Dmel.lib";

        print OUT1 "crm\t";
        for my $crm (sort @crmNames){
            print OUT1 "$crm\t";
        }
        print OUT1 "\n";

        for my $id (sort keys %{$output{"train"}})
        {
            print OUT1 "$id\t";
            if ($train_seq_label{$id} eq "1")
            {
                print TRAIN_lib "1 ";
                if ($id =~ /^Dmel/)
                {
                    print TRAIN_lib_DMEL "1 ";
                }
            }
            else
            {
                print TRAIN_lib "-1 ";
                if ($id =~ /^Dmel/)
                {
                    print TRAIN_lib_DMEL "-1 ";
                }
            }

            my $feature_count = 1;
            for my $crm (sort keys %{$output{"train"}{$id}})
            {
                my @score = split /\s+/,$output{"train"}{$id}{$crm};
                print OUT1 "$score[1]\t";
                print TRAIN_lib "$feature_count\:$score[1] ";
                if ($id =~ /^Dmel/)
                {
                    print TRAIN_lib_DMEL "$feature_count\:$score[1] ";
                }
                $feature_count++;
            }
            print OUT1 "\n";
            print TRAIN_lib "\n";
            
            if ($id =~ /^Dmel/)
            {
                print TRAIN_lib_DMEL "\n";
            }
        }
        close OUT1;
        close TRAIN_lib;
        close TRAIN_lib_DMEL;

        ##======================================
        # output test seq msIMM features to file 
        #=======================================
        open OUT2,">$homeDir/test.ensembFeat";
        open TEST_lib,">$homeDir/test.ensembFeat.lib";
        open TEST_lib_DMEL,">$homeDir/test.ensembFeat.Dmel.lib";

        print OUT2 "crm\t";
        for my $crm (sort @crmNames){
            print OUT2 "$crm\t";
        }
        print OUT2 "\n";

        for my $id (sort keys %{$output{"test"}})
        {
            print OUT2 "$id\t";
            if ($test_seq_label{$id} eq "1")
            {
                print TEST_lib "1 ";
                print TEST_lib_DMEL "1 ";
            }
            else
            {
                print TEST_lib "-1 ";
                print TEST_lib_DMEL "-1 ";
            }

            my $feature_count = 1;
            for my $crm (sort keys %{$output{"test"}{$id}})
            {
                my @score = split /\s+/,$output{"test"}{$id}{$crm};
                print OUT2 "$score[1]\t";
                print TEST_lib "$feature_count\:$score[1] ";
                print TEST_lib_DMEL "$feature_count\:$score[1] ";
                $feature_count++;
            }
            print OUT2 "\n";
            print TEST_lib "\n";
            print TEST_lib_DMEL "\n";
        }
        close OUT2;
        close TEST_lib;
        close TEST_lib_DMEL;
}

sub read_labels {
    my $train_label_file = "$modelDir/$crm/time$timeNum/fold$foldNum/train.label";
    my $test_label_file = "$modelDir/$crm/time$timeNum/fold$foldNum/test.label";

    open LAB,$train_label_file or die "cannot open $train_label_file";
    while (<LAB>){
        chomp(my $line = $_);
        my ($id, $lab) = split /\s+/,$line;
        $train_seq_label{$id} = $lab;
    }
    close LAB;

    open LAB,$test_label_file or die "cannot open $test_label_file";
    while (<LAB>){
        chomp(my $line = $_);
        my ($id, $lab) = split /\s+/,$line;
        $test_seq_label{$id} = $lab;
    }
    close LAB;
}

# call functions
getNames();
pred();

