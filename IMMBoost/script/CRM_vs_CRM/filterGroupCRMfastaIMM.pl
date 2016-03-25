use strict;
use warnings;
use lib '/shared-mounts/sinhas/lib/';
use File::Basename;
use Bio::SeqIO;

die "Usage: perl $0 CRMDir OutDir CRMGroupTable\n" unless @ARGV==3;

my $CRMDir = $ARGV[0];
my $outdir = $ARGV[1];
my $CRMGroup = $ARGV[2];

##======== store enhancer sequence ID for each CRM set ========##
my %crm2seqID = ();
while (<$CRMDir/*>){
    chomp(my $crm = $_);
    my $crmName = basename($crm);
    open IN,"$crm/fasta/crms.fasta.lenGC.newID.txt";
    while (<IN>){
        chomp(my $line = $_);
        my ($id, $len, $gc) = split /\s+/,$line;
        $crm2seqID{$crmName}{$id} = 1;
    }
    close IN;
}

##======== store CRM IDs that are in the same group ========##
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

##======== store enhancer sequence ID for a CRM group =========##
my %groupCRM2seqID = ();
for my $crm (keys %crm2group){
    for my $member (keys %{$crm2group{$crm}}){
        for my $seqID (keys %{$crm2seqID{$member}}){
            $groupCRM2seqID{$crm}{$seqID} = 1;
        }
    }
}

while (<$outdir/*>){
    chomp(my $crm = $_);
    my $crmName = basename($crm);
    warn "======= $crmName ========\n";
    for (my $k=1;$k<=10;$k++){
        warn "$k ....\n";
        for (my $i=1;$i<=5;$i++){
            open OUT,">$crm/time$k/fold$i/train.neg.filGroup.fasta";
            my $negFa = "$crm/time$k/fold$i/train.neg.fasta"; 
            my $fa = Bio::SeqIO->new(-file=>$negFa,-format=>'Fasta');
            while (my $curSeq = $fa->next_seq()){
                my $id = $curSeq->id();
        # warn "$id\n";
                my $seqID = (split(/\_/,$id,2))[1];
        #  warn "$seqID\n";
                my $seq = $curSeq->seq();
                if (!exists $groupCRM2seqID{$crmName}{$seqID}){
                    print OUT ">$id\n$seq\n";
                }
            }
            close OUT;
        }
    }
}
