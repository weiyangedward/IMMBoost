use strict;
use warnings;
use File::Basename;
use Bio::SeqIO;

die "Usage: perl $0 CRMDir OutDir CRMGroupTable\n" unless @ARGV==4;

my $CRMDir = $ARGV[0];
my $outdir = $ARGV[1];
my $CRMGroup = $ARGV[2];
my $times = $ARGV[3];

my %CRMsets = ();
my %crm2seqID = ();
while (<$CRMDir/*>)
{
    chomp(my $crm = $_);
    my $crmName = basename($crm);

    $CRMsets{$crmName} = 1;
    my $crmFa = "$CRMDir/$crmName/fasta/CRM.fasta";
    my $fa1 = Bio::SeqIO->new(-file=>$crmFa,-format=>'Fasta');
    while (my $curSeq = $fa1->next_seq())
    {
        my $id = $curSeq->id();
        my $seq = $curSeq->seq();
        $crm2seqID{$crmName}{$id} = 1;
    }
}

my %crm2group = ();
open IN,$CRMGroup or die "cannot open $CRMGroup";
while (<IN>){
    chomp(my $line = $_);
    my @group = split /\s+/,$line;
    if ($group[1]){
        for my $crmID (@group[1..$#group]){
            $crm2group{$group[0]}{$crmID} = 1;
        } 
    }
}
close IN;

my %groupCRM2seqID = ();
for my $crm (keys %crm2group){
    for my $member (keys %{$crm2group{$crm}}){
        for my $seqID (keys %{$crm2seqID{$member}}){
            $groupCRM2seqID{$crm}{$seqID} = 1;
        }
    }
}

# while (<$outdir/*>){
for my $crmName (keys %CRMsets)
{
    # chomp(my $crm = $_);
    # my $crmName = basename($crm);
    warn "$crmName\n";
    for (my $k=1;$k<=$times;$k++)
    {
        warn "time $k ...\n";
        for (my $i=1;$i<=5;$i++)
        {
            my $labelFileTrain = "$outdir/$crmName/time$k/fold$i/train.label";
            my $labelFileTest = "$outdir/$crmName/time$k/fold$i/test.label";

            open LAB,$labelFileTrain or die "cannot open $labelFileTrain";
            my %seq2labTrain = ();
            while (<LAB>){
                chomp(my $line = $_);
                my ($id, $lab) = split /\s+/,$line;
                $seq2labTrain{$id} = $lab;
            }
            close LAB;

            open LAB,$labelFileTest or die "cannot open $labelFileTest";
            my %seq2labTest = ();
            while (<LAB>){
                chomp(my $line = $_);
                my ($id, $lab) = split /\s+/,$line;
                $seq2labTest{$id} = $lab;
            }
            close LAB;
            ##====== Create lib file for training data ========##
            my $trainFile = "$outdir/$crmName/time$k/fold$i/train.ensembFeat";
            open IN,$trainFile or die "cannot open $trainFile";
            open OUT1,">$outdir/$crmName/time$k/fold$i/train.ensembFeat.filGroup2.Dmel";
            open OUT2,">$outdir/$crmName/time$k/fold$i/train.ensembFeat.filGroup2.Dmel.lib";
            while (<IN>){
                chomp(my $line = $_);
                if ($line =~ /^crm/){
                    print OUT1 "$line\n";
                }
                else{
                    if ($line =~ /^Dmel/){
                        my @array = split /\s+/,$line;
                        my $seqID = (split(/\_/,$array[0],2))[1];
                        if ($seq2labTrain{$array[0]} == 1){
                            print OUT1 "$line\n";
                            print OUT2 "1 ";
                            for (my $j=1; $j<=$#array; $j++){
                                print OUT2 "$j\:$array[$j] ";
                            }
                            print OUT2 "\n";
                        }
                        else{
                            if (!exists $groupCRM2seqID{$crmName}{$seqID}){
                                print OUT1 "$line\n";
                                print OUT2 "-1 ";
                                for (my $j=1; $j<=$#array; $j++){
                                    print OUT2 "$j\:$array[$j] ";
                                }
                                print OUT2 "\n";
                            }
                        }
                    }
                }
            }
            close OUT1;
            close OUT2;
            close IN;
 
        }
    }
}
