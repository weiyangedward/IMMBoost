=head1 Description
  
  This script uses CRM grouping to filter out closely
   related seq from training data.

=head1 Usage

  perl filterGroupCRM.pl CRMDir OutDir CRMGroupTable

=cut

use strict;
use warnings;
use File::Basename;
use Bio::SeqIO;

# die "Usage: perl $0 CRMDir OutDir CRMGroupTable\n" unless @ARGV==3;
die `pod2text $0` if (@ARGV!=3);

my $CRMDir = $ARGV[0];
my $outdir = $ARGV[1];
my $CRMGroup = $ARGV[2];

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

    # open IN,"$crm/fasta/crms.fasta.lenGC.newID.txt";
    # while (<IN>){
    #     chomp(my $line = $_);
    #     my ($id, $len, $gc) = split /\s+/,$line;
    #     $crm2seqID{$crmName}{$id} = 1;
    # }
    # close IN;
}

my %crm2group = ();
open IN,$CRMGroup;
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
    for (my $k=1;$k<=10;$k++)
    {
        warn "trial $k\n";
        for (my $i=1;$i<=5;$i++)
        {
            my $labelFileTrain = "$outdir/$crmName/time$k/fold$i/train.label";
            my $labelFileTest = "$outdir/$crmName/time$k/fold$i/test.label";

            open LAB,$labelFileTrain;
            my %seq2labTrain = ();
            while (<LAB>){
                chomp(my $line = $_);
                my ($id, $lab) = split /\s+/,$line;
                $seq2labTrain{$id} = $lab;
            }
            close LAB;

            open LAB,$labelFileTest;
            my %seq2labTest = ();
            while (<LAB>){
                chomp(my $line = $_);
                my ($id, $lab) = split /\s+/,$line;
                $seq2labTest{$id} = $lab;
            }
            close LAB;
            ##====== Create lib file for training data ========##
            my $trainFile = "$outdir/$crmName/time$k/fold$i/train.ensembFeat";
            open IN,$trainFile;
            open OUT1,">$outdir/$crmName/time$k/fold$i/train.ensembFeat.filGroup2";
            open OUT2,">$outdir/$crmName/time$k/fold$i/train.ensembFeat.filGroup2.lib";
            while (<IN>){
                chomp(my $line = $_);
                if ($line =~ /^crm/){
                    print OUT1 "$line\n";
                }
                else{
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
            close OUT1;
            close OUT2;
            close IN;
            ##====== create lib format for test data =========##
            my $testFile = "$outdir/$crmName/time$k/fold$i/test.ensembFeat";
            open IN,$testFile;
            open OUT1,">$outdir/$crmName/time$k/fold$i/test.ensembFeat.filGroup2";
            open OUT2,">$outdir/$crmName/time$k/fold$i/test.ensembFeat.filGroup2.lib";
            while (<IN>){
                chomp(my $line = $_);
                if ($line =~ /^crm/){
                    print OUT1 "$line\n";
                }
                else{
                    my @array = split /\s+/,$line;
                   #  my $seqID = (split(/\_/,$array[0],2))[1];
                    if ($seq2labTest{$array[0]} == 1){
                        print OUT1 "$line\n";
                        print OUT2 "1 ";
                        for (my $j=1; $j<=$#array; $j++){
                            print OUT2 "$j\:$array[$j] ";
                        }
                        print OUT2 "\n";
                    }
                    else{
                    #    if (!exists $groupCRM2seqID{$crmName}{$seqID}){
                            print OUT1 "$line\n";
                            print OUT2 "-1 ";
                            for (my $j=1; $j<=$#array; $j++){
                                print OUT2 "$j\:$array[$j] ";
                            }
                            print OUT2 "\n";
                     #   }
                    }
                }
            }
            close OUT1;
            close OUT2;
            close IN;
        }
    }
}
