=head1 Description
  
  This script generates negative background seq
  with similar length and GC content to query CRMset.

=head1 Usage

  perl prepareModelAndData.pl CRM_seq acc_seq how_many_neg_seq > neg_bkg_seq

=cut

use strict;
use warnings;
use Bio::SeqIO;
use List::Util qw(shuffle);
use POSIX;

die `pod2text $0` if (@ARGV!=3);

my $CRM_seq = $ARGV[0];
my $acc_seq = $ARGV[1];
my $how_many_neg_seq = $ARGV[2];

my %crm_len = ();
my %crm_gc = ();
my %acc_len = ();
my %acc_gc = ();
my %acc_seq = ();
my @acc_seq = ();
my $acc_amount = 0;
my $crm_amount = 0;

sub readin_CRM {
	my $fa = Bio::SeqIO->new(-file=>$CRM_seq,-format=>'Fasta');
	while (my $nextSeq = $fa->next_seq()){
		$crm_amount +=1;
        my $id = $nextSeq->id();
        my $seq = $nextSeq->seq();
        my $len = length($seq);
        my $N = $seq =~ tr/Nn//; # gaps
        my $len_no_gap = $len - $N;
        my $gc = $seq =~ tr/GCgc//;
        $gc = $gc / $len_no_gap;
        $crm_len{$id} = $len;
        $crm_gc{$id} = $gc;
    }
}

sub readin_acc_seq {
	my $fa = Bio::SeqIO->new(-file=>$acc_seq,-format=>'Fasta');
	while (my $nextSeq = $fa->next_seq()){
		$acc_amount +=1;
        my $id = (split /\s+/,($nextSeq->id()))[0];
        my $seq = $nextSeq->seq();
        $acc_seq{$id} = $seq;
        my $len = length($seq);
        my $N = $seq =~ tr/Nn//; # gaps
        my $len_no_gap = $len - $N;
        my $gc = $seq =~ tr/GCgc//;
        $gc = $len_no_gap > 0 ? $gc / $len_no_gap : 0;
        $acc_len{$id} = $len;
        $acc_gc{$id} = $gc;
        push(@acc_seq,$id);
    }
}

sub get_bgk_seq {
	my $count = 0;
    my $iter = 0;
    my %seen = ();
    my @acc_seq_shuf = shuffle @acc_seq;

    while ($count < $how_many_neg_seq && $iter < $how_many_neg_seq * $crm_amount)
    {
	    for my $id1 (keys %crm_len)
	    {
	        $iter +=1;
	        last if $count >= $how_many_neg_seq; # exit when 100 seq is found

	        for my $id2 (@acc_seq_shuf)
	        {
	            if ( ($acc_len{$id2} < $crm_len{$id1}*1.05 && $acc_len{$id2} > $crm_len{$id1}*0.95) && 
	            	($acc_gc{$id2} < $crm_gc{$id1}*1.05 && $acc_gc{$id2} > $crm_gc{$id1}*0.95) && 
	            	!$seen{$id2} )
	            {
	                # warn "$id2\t$acc_len{$id2}\t$acc_gc{$id2}\n";
	                print ">$id2\n$acc_seq{$id2}\n";
	                $seen{$id2} = 1;
	                $count +=1;
	                last; # exit searching for neg given query CRM
	            }
	        }
	    }
	}
}

readin_CRM();
readin_acc_seq();
get_bgk_seq();