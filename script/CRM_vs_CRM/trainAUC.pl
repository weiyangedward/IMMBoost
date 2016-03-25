use strict;
use warnings;
use lib '/shared-mounts/sinhas/lib/';
use Bio::SeqIO;
use POSIX;
use FindBin qw($Bin);

die "Usage: perl mapping outdir\n" unless @ARGV==2;

my $tmpDir = $ENV{TMPDIR};
my $tmpModelDir = "$tmpDir/trainAUC";
`mkdir $tmpModelDir` unless (-e "$tmpModelDir");

my $crm = $ARGV[0];
my $outdir = $ARGV[1];

for (my $k=1;$k<=10;$k++){
for (my $i=1;$i<=5;$i++){
    my $homeDir = "$outdir/$crm/time$k/fold$i";
    `cp $homeDir/train.crm.fasta $tmpModelDir/train.crm.fasta`;
    `cp $homeDir/train.neg.filGroup.fasta $tmpModelDir/train.neg.filGroup.fasta`;
    `cp $homeDir/train.crm.fasta.andNeg.rc.6mer $tmpModelDir/train.crm.fasta.andNeg.rc.6mer`;
    `cp $homeDir/train.ensembFeat.filGroup2 $tmpModelDir/train.ensembFeat.filGroup2`;
    `cp $homeDir/train.ensembFeat.filGroup2.lib $tmpModelDir/train.ensembFeat.filGroup2.lib`;

    my $crmTrainFa = "$tmpModelDir/train.crm.fasta";
    my $negTrainFa = "$tmpModelDir/train.neg.filGroup.fasta";
    my $trainKmer = "$tmpModelDir/train.crm.fasta.andNeg.rc.6mer";
    my $trainIMMpred = "$tmpModelDir/train.ensembFeat.filGroup2";
    my $trainIMMpredLib = "$tmpModelDir/train.ensembFeat.filGroup2.lib";
    
    warn "========== read pos train fasta ============\n";
    my $crmPos = 0;
    my %crmSeq = ();
    my @crmSeqID = ();
    my %crmNum = ();
    my %crmGeneralID2id = ();
    my %crmID2pos = ();
    my $fa1 = Bio::SeqIO->new(-file=>$crmTrainFa,-format=>'Fasta');
    while (my $nextSeq = $fa1->next_seq()){
        my $id = $nextSeq->id();
        my $seq = $nextSeq->seq();
        $crmSeq{$id} = $seq;
        my $generalID = $2 if $id =~ /(\w+?)_(\S+)/;
#        push @crmSeqID,$generalID;         ## order should match kmer
#        $crmNum ++;
#        $crmNum{$generalID} = 1;
        if (!$crmNum{$generalID}){
            push @crmSeqID,$generalID;
            $crmNum{$generalID} = 1;
        }
        $crmID2pos{$id} = $crmPos;
        $crmPos ++;
        push @{$crmGeneralID2id{$generalID}},$id;
    }
    my $crmNumLen = keys %crmNum;
    warn "======== read neg train fasta =======\n";
    my $negPos = 0;
    my %negSeq = ();
    my @negSeqID = ();
    my %negNum = ();
    my %negGeneralID2id = ();
    my %negID2pos = ();
    my $fa2 = Bio::SeqIO->new(-file=>$negTrainFa,-format=>'Fasta');
    while (my $nextSeq = $fa2->next_seq()){
        my $id = $nextSeq->id();
        my $seq = $nextSeq->seq();
        $negSeq{$id} = $seq;
        my $generalID = $2 if $id =~ /(\w+?)_(\S+)/;
#        push @negSeqID,$generalID;         ## order should match kmer
#        $negNum ++;
#        $negNum{$generalID} = 1;
        if (!$negNum{$generalID}){
            $negNum{$generalID} = 1;
            push @negSeqID,$generalID;
        }
        $negID2pos{$id} = $negPos;
        $negPos ++;
        push @{$negGeneralID2id{$generalID}},$id;
    }
    my $negNumLen = keys %negNum;
    warn "========= read kmer ========\n";
    open IN,$trainKmer;
    my @posKmer = ();
    my @negKmer = ();
    while (my $line = <IN>){
        chomp($line);
        my ($lab, $feat) = split(/\s+/,$line,2);
        if ($lab == 1){
            push @posKmer,$line;    ## order should match id
        }
        else{
            push @negKmer,$line;    ## order should match id
        }
    }
    close IN;
    warn "========= read IMM predict scores ========\n";
    my %seqID2score = ();
    my %seqID2pos = ();
    my $count = 0;
    open IN,$trainIMMpred;
    while (my $line = <IN>){
        chomp($line);
        next if $line =~ /^crm/;
        my ($seqID, $score) = split(/\s+/,$line,2);
        $seqID2score{$seqID} = $line;
        $seqID2pos{$seqID} = $count;
        $count ++;
    }
    close IN;
    
    warn "========= read IMM predict scores lib ========\n";
    my @lib2score = ();
    open IN,$trainIMMpredLib;
    while (my $line = <IN>){
        chomp($line);
        push @lib2score,$line;
    }
    close IN;
    
    warn "======= generate 5-fold cross validation sets =======\n";
    my $crmStart = 0;
    my $crmEnd = 0;
    
    my $negStart = 0;
    my $negEnd = 0;
    
    `mkdir $homeDir/trainAUC` unless (-d "$homeDir/trainAUC");
    for (my $j=1;$j<=5;$j++){
        `mkdir $tmpModelDir/fold$j` unless (-d "$tmpModelDir/fold$j");
        $crmEnd = floor(($j/5) * ($#crmSeqID + 1)) - 1;
        $negEnd = floor(($j/5) * ($#negSeqID + 1)) - 1;
        
        open OUT1,">$tmpModelDir/fold$j/test.crm.fasta";
        open OUT2,">$tmpModelDir/fold$j/train.crm.fasta";

        open OUT3,">$tmpModelDir/fold$j/test.label";
        open OUT4,">$tmpModelDir/fold$j/train.label";

        open OUT5,">$tmpModelDir/fold$j/train.crm.fasta.andNeg.rc.6mer";
        open OUT6,">$tmpModelDir/fold$j/test.crm.fasta.andNeg.rc.6mer";

        open OUT7,">$tmpModelDir/fold$j/train.ensembFeat.filGroup2";
        open OUT8,">$tmpModelDir/fold$j/test.ensembFeat.filGroup2";

        open OUT9,">$tmpModelDir/fold$j/train.ensembFeat.filGroup2.lib";
        open OUT10,">$tmpModelDir/fold$j/test.ensembFeat.filGroup2.lib";

        
        warn "======= output pos data ============\n";
        for (my $n=0; $n<=$#crmSeqID;$n++){
            my $generalID = $crmSeqID[$n];
            if ($n <= $crmEnd && $n >= $crmStart){
                ##======= output test data ========##
                for my $id (@{$crmGeneralID2id{$generalID}}){
                    print OUT1 ">$id\n$crmSeq{$id}\n";
                    print OUT3 "$id\t1\n";
                    my $pos = $crmID2pos{$id};
                    my $libPos = $seqID2pos{$id};
                    print OUT6 "$posKmer[$pos]\n";
                    print OUT8 "$seqID2score{$id}\n";
                    print OUT10 "$lib2score[$libPos]\n";
                }
#                print OUT1 ">$crmSeqID[$n]\n$crmSeq{$crmSeqID[$n]}\n";
#                print OUT3 "$crmSeqID[$n]\t1\n";
#                print OUT6 "$posKmer[$n]\n";
#                print OUT8 "$seqID2score{$crmSeqID[$n]}\n";
            }
            else{
                ##======== output train data ========##
                for my $id (@{$crmGeneralID2id{$generalID}}){
                    print OUT2 ">$id\n$crmSeq{$id}\n";
                    print OUT4 "$id\t1\n";
                    my $pos = $crmID2pos{$id};
                    my $libPos = $seqID2pos{$id};
                    print OUT5 "$posKmer[$pos]\n";
                    print OUT7 "$seqID2score{$id}\n";
                    print OUT9 "$lib2score[$libPos]\n";
                }
#                print OUT2 ">$crmSeqID[$n]\n$crmSeq{$crmSeqID[$n]}\n";
#                print OUT4 "$crmSeqID[$n]\t1\n";
#                print OUT5 "$posKmer[$n]\n";
#                print OUT7 "$seqID2score{$crmSeqID[$n]}\n";
            }
        }
        close OUT1;
        close OUT2;

        open OUT1,">$tmpModelDir/fold$j/test.neg.fasta";
        open OUT2,">$tmpModelDir/fold$j/train.neg.fasta";
        warn "========= output neg data =============\n";
        for (my $n=0;$n<=$#negSeqID;$n++){
            my $generalID = $negSeqID[$n];
            if ($n <= $negEnd && $n >= $negStart){
                ##========= output test data =======##
                for my $id (@{$negGeneralID2id{$generalID}}){ 
                    print OUT1 ">$id\n$negSeq{$id}\n";
                    print OUT3 "$id\t-1\n";
                    my $pos = $negID2pos{$id};
                    my $libPos = $seqID2pos{$id};
                    print OUT6 "$negKmer[$pos]\n";
                    print OUT8 "$seqID2score{$id}\n";
                    print OUT10 "$lib2score[$libPos]\n";
                }
#                print OUT1 ">$negSeqID[$n]\n$negSeq{$negSeqID[$n]}\n";
#                print OUT3 "$negSeqID[$n]\t-1\n";
#                print OUT6 "$negKmer[$n]\n";
#                print OUT8 "$seqID2score{$negSeqID[$n]}\n";
            }
            else{
                ##======== output train data ==========##
                for my $id (@{$negGeneralID2id{$generalID}}){ 
                    print OUT2 ">$id\n$negSeq{$id}\n";
                    print OUT4 "$id\t-1\n";
                    my $pos = $negID2pos{$id};
                    my $libPos = $seqID2pos{$id};
                    print OUT5 "$negKmer[$pos]\n";
                    print OUT7 "$seqID2score{$id}\n";
                    print OUT9 "$lib2score[$libPos]\n";
                }
#                print OUT2 ">$negSeqID[$n]\n$negSeq{$negSeqID[$n]}\n";
#                print OUT4 "$negSeqID[$n]\t-1\n";
#                print OUT5 "$negKmer[$n]\n";
#                print OUT7 "$seqID2score{$negSeqID[$n]}\n";
            }
        }
        close OUT1;
        close OUT2;
        close OUT3;
        close OUT4;
        close OUT5;
        close OUT6;
        close OUT7;
        close OUT8;
        close OUT9;
        close OUT10;
        
        $crmStart = $crmEnd + 1;
        $negStart = $negEnd + 1;
        `cat $tmpModelDir/fold$j/test.crm.fasta $tmpModelDir/fold$j/test.neg.fasta > $tmpModelDir/fold$j/test.crm.and.neg.fasta`;
    }

    warn "========== IMM prediction ==============\n";
    `perl $Bin/msIMM.model.trainAUC.pl $tmpModelDir`;
   #
    warn "========== kmerSVM prediction ==========\n";
    `perl $Bin/kmerSVM.trainAUC.pl $tmpModelDir`;
   #
    warn "========== RF+IMM prediction =============\n";
    `perl $Bin/svm.v5.RFmodel.trainAUC.pl $tmpModelDir`;

    warn "========== SVM+IMM prediction ==========\n";
    `perl $Bin/svm.v5.SVMmodel.libsvm.trainAUC.pl $tmpModelDir`;

    `tar -zcvf $tmpDir/trainAUC.tar.gz -C $tmpModelDir .`;
    if (-e "$homeDir/trainAUC/"){
        `rm -rf $homeDir/trainAUC/*`;
#        `mv $tmpModelDir/* $homeDir/trainAUC/`;
        `mv $tmpDir/trainAUC.tar.gz $homeDir/trainAUC/`;
        `tar -zxvf $homeDir/trainAUC/trainAUC.tar.gz -C $homeDir/trainAUC/`;
    }
    else{
        `mv $tmpDir/trainAUC.tar.gz $homeDir/trainAUC/`;
        `tar -zxvf $homeDir/trainAUC/trainAUC.tar.gz -C $homeDir/trainAUC/`;
#        `mv $tmpModelDir/* $homeDir/trainAUC/`;
    }
}
}
