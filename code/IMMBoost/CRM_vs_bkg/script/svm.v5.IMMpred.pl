=head1

This script is used to generate msIMM score as features
for each seq in test and training sets. These new features
will be used to train a new model to classify CRM seq.
Note that BioPerl is needed.

=cut

use strict;
use warnings;
use lib '/shared-mounts/sinhas/lib/';
use File::Basename;
use Bio::SeqIO;

die "Usage: perl\n\t$0 CRMname OutDir timeNum foldNum\n" unless @ARGV==4;

my $tmpDir = $ENV{TMPDIR};
my $tmpModelDir = "$tmpDir/model";
`mkdir $tmpModelDir` unless (-e "$tmpModelDir");

my $pred = "/home/weiyang4/software/SCRMshaw/bin/imm_score";
my $train = "/home/weiyang4/software/SCRMshaw/bin/imm_build";

my $crm = $ARGV[0]; # CRMname
my $modelDir = $ARGV[1]; # outdir
my $timeNum = $ARGV[2]; # trial
my $foldNum = $ARGV[3]; # fold

# get CRMnames
my @crmNames;
sub getNames {
    while (<$modelDir/*>)
    {
        chomp(my $crmDir = $_);
        my $crmName = basename($crmDir);
        push @crmNames,$crmName; # store CRMname in array
    }
}

# generate msIMM features for each seq in test and training sets
sub pred {
    my $crmDir = "$modelDir/$crm";
    my $k = $timeNum;
    my $i = $foldNum;    
    my $trainFa = "$crmDir/time$k/fold$i/train.crm.and.neg.fasta"; # training seq
    my $testFa = "$crmDir/time$k/fold$i/test.crm.and.neg.fasta"; # test seq
        
    ##==== hash training seq =====
    my %trainSeq = ();
    my $fa1 = Bio::SeqIO->new(-file=>$trainFa,-format=>'Fasta');
    while (my $curSeq = $fa1->next_seq()){
        my $id = $curSeq->id();
        my $seq = $curSeq->seq();
        $trainSeq{$id} = $seq; # hash seq e.g., msCRM id = Dyak_ftz_CE8024
    }
    
    ##====== hash test data ========
    my %testSeq = ();
    my $fa2 = Bio::SeqIO->new(-file=>$testFa,-format=>'Fasta');
    while (my $curSeq = $fa2->next_seq()){
        my $id = $curSeq->id();
        my $seq = $curSeq->seq();
        $testSeq{$id} = $seq; ## dmel ID, e.g., cad_14_construct
    }

    ##===== train msIMMs and score training and test seq =====
    my %output = ();   ## store all of the IMM output
    # loop through all CRMname
    for my $otherCRM (@crmNames)
    {
        my $crmFa = "$modelDir/$otherCRM/allData/msCRM.fa";
        my $negFa = "$modelDir/$otherCRM/allData/msNeg.fa";
            
        ##====== hash msCRM seq ==========
        my %crmSeq = ();
        my %crmGeneralID = ();
        my $fa3 = Bio::SeqIO->new(-file=>$crmFa,-format=>'Fasta');
        while (my $curSeq = $fa3->next_seq())
        {
            my $id = $curSeq->id();
            my $seq = $curSeq->seq();
            $crmSeq{$id} = $seq; ## msCRM id
        }

        ##===== hash msNeg seq =======
        my %negSeq = ();
        my $fa4 = Bio::SeqIO->new(-file=>$negFa,-format=>'Fasta');
        while (my $curSeq = $fa4->next_seq())
        {
            my $id = $curSeq->id();
            my $seq = $curSeq->seq();
            $negSeq{$id} = $seq;
        }
        ##==== for each training seq, generate pos and neg training seq to train a msIMM, and then score this seq
        ##====
        # loop through each training seq
        for my $id (sort keys %trainSeq)
        {
            # output target seq as test data
            open OUT,">$tmpModelDir/$id.fa";
            print OUT ">$id\n$trainSeq{$id}\n";
            close OUT;
            # convert msCRM id to CRM id
            my $thisSeqGeneralID = $id;  ## in case sequence is from negative data, which does not have msID
            if ($id =~ /(\w+?)_(\S+)/)
            {
                $thisSeqGeneralID = $2;
            }

            #==== create CRM training seq for msIMM, where seq that are the same as target seq or in test data are filtered since target seq will be used to train a model to pred on test data
            #=====
            open OUT,">$tmpModelDir/$otherCRM.crms.fa";
            # loop through each CRM seq in a CRM set
            for my $crmID (sort keys %crmSeq)
            {
                my $generalID = $2 if $crmID =~ /(\w+?)_(\S+)/;
                # add a CRM seq to training set if it is not the same as target seq nor in test data since target seq will be used to train a model to pred on test data
                if (($thisSeqGeneralID ne $generalID) && !$testSeq{$generalID})
                {
                    print OUT ">$crmID\n$crmSeq{$crmID}\n";
                }
            }
            close OUT;

            #=== Create neg training seq for msIMM, where seq that are the same as target seq or in test data are filtered
            #====
            open OUT,">$tmpModelDir/$otherCRM.neg.fa";
            # loop through each neg seq
            for my $negID (sort keys %negSeq)
            {
                if (($negID ne $id) && !$testSeq{$negID})
                {
                    print OUT ">$negID\n$negSeq{$negID}\n";
                }
            }
            close OUT;

            ##==== Train msIMM and score target seq ====
            `$train -r -k 6 < $tmpModelDir/$otherCRM.neg.fa > $tmpModelDir/$otherCRM.neg.model`;
            `$train -r -k 6 < $tmpModelDir/$otherCRM.crms.fa > $tmpModelDir/$otherCRM.crms.model`;
            my $crmModel = "$tmpModelDir/$otherCRM.crms.model";
            my $negModel = "$tmpModelDir/$otherCRM.neg.model";
            ## hash pred score
            $output{"train"}{$id}{$otherCRM} = `$pred -f -n $negModel $crmModel $tmpModelDir/$id.fa`;
            ## delete everything in the dir at slave node
            `rm -rf $tmpModelDir/*`;
        }

        ##==== for each test seq, generate pos and neg training seq to train a msIMM, and then score this seq
        ##====
        # loop through each test seq
        for my $id (sort keys %testSeq)
        {
            # target seq as output test data
            open OUT,">$tmpModelDir/$id.fa";
            print OUT ">$id\n$testSeq{$id}\n";
            close OUT;
            #==== create CRM training seq for msIMM, where seq that are the same as target seq or in test data are filtered since target seq will be used to train a model to pred on test data
            #=====
            open OUT,">$tmpModelDir/$otherCRM.crms.fa";
            for my $crmID (sort keys %crmSeq)
            {
                my $generalID = $2 if $crmID =~ /(\w+?)_(\S+)/;
                if ($id ne $generalID && !$testSeq{$generalID})
                {
                    print OUT ">$crmID\n$crmSeq{$crmID}\n";
                }
            }
            close OUT;

             #=== Create neg training seq for msIMM, where seq that are the same as target seq or in test data are filtered
            #====
            open OUT,">$tmpModelDir/$otherCRM.neg.fa";
            for my $negID (sort keys %negSeq)
            {
                if ($negID ne $id && !$testSeq{$negID})
                {
                    print OUT ">$negID\n$negSeq{$negID}\n";
                }
            }
            close OUT;

            ##==== Train msIMM and score target seq ====
            `$train -r -k 6 < $tmpModelDir/$otherCRM.neg.fa > $tmpModelDir/$otherCRM.neg.model`;
            `$train -r -k 6 < $tmpModelDir/$otherCRM.crms.fa > $tmpModelDir/$otherCRM.crms.model`;
            my $crmModel = "$tmpModelDir/$otherCRM.crms.model";
            my $negModel = "$tmpModelDir/$otherCRM.neg.model";
            # Store prediction into a hash table
            $output{"test"}{$id}{$otherCRM} = `$pred -f -n $negModel $crmModel $tmpModelDir/$id.fa`;
            # delete everything in tmp directory
            `rm -rf $tmpModelDir/*`;
        }
    }
    ##====== output training seq msIMM features to file ====
    open OUT1,">$tmpModelDir/train.ensembFeat";
    print OUT1 "crm\t";
    for my $crm (sort @crmNames)
    {
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
    open OUT2,">$tmpModelDir/test.ensembFeat";
    print OUT2 "crm\t";
    for my $crm (sort @crmNames)
    {
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

    ##====== Move feature files to master dir ======
    `mv $tmpModelDir/test.ensembFeat $crmDir/time$k/fold$i/test.ensembFeat`;
    `mv $tmpModelDir/train.ensembFeat $crmDir/time$k/fold$i/train.ensembFeat`;
    
}

# call functions
getNames();
pred();

