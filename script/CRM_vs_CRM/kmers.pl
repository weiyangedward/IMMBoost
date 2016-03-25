#!/usr/bin/perl

=head1 Description

	This script inputs a FASTA file and calculates 
	the k-mer counts for each sequence in the file. 
	k-mers are interpreted as tertiary codes (A=0,C=1,G=2,T=3) 
	and are converted to a base 10 integer. In the k-mer 
	file format, the column indicates the particular 
	k-mer and the value its corresponding count. 
	There are 3 columns at the beginning of each line
	(an id, sub-id, and count of base pairs in sequence) 
	before columns correspond to k-mers. A FASTA file 
	with multiple sequences will have a line of counts 
	for each sequence.

=head1 Usage

	perl kmers.pl [OPTIONS] FILE K

	FILE    Input FASTA file.
	K       k-mer length

	--help       Print help and exit
	--rc          Do reverse complement

=cut

use strict;
use FindBin qw($Bin $Script);
use lib '$Bin/../lib/';
# use lib '/shared-mounts/sinhas/lib/';
use File::Basename;
use Bio::SeqIO;
use Getopt::Long;

## Option triggers
my $help	= undef;			##	Help option
my $rc		= undef;			##	Do reverse complement. No by default.

## Get Options
GetOptions(
	"h|help" => \$help,
	"r|rc"   => \$rc
);

$rc	= defined($rc) ? 1 : 0;	##	Cast reverse complement to 1(Yes) or (0) No.
my $ifile = $ARGV[0]; 			##	FASTA file
my $k    	= $ARGV[1];				##	K mer length

## Check exit conditions
die `pod2text $0` if (@ARGV != 2 || $help);

## Iterate over each sequence in the FASTA file.
my $fasta = Bio::SeqIO->new(-file=>$ifile, -format=>'Fasta');
my $i 		= 0;
while (my $fh = $fasta->next_seq()) {
    my $seq 	= $fh->seq();
		my $len		=	CountBases($seq);
    my $distr = KDistr($seq, $k, $rc);
		print "$i\t-\t$len\t$distr\n";
    $i++;
}

sub CountBases {
	my $seq 	= shift;	
	my $bases = $seq;
	$bases 		=~ s/N//g;
	return length($bases);
}

sub KDistr {
	## Variables
	my $seq 	= lc(shift);	## Lowercase FASTA.
	my $k 		= shift;			## K-mer length.
	my $rc		= shift;			## Do rc.
	my $kmer 	= 0;					## K-mer code.
	my $f  		= 0; 					## Will be set to i if the ith position from right is a non 0/1/2/3.
	my @bin 	= ();					## Binary sequence. [acgt-] = [01234]
	my %hash	= ();					## hash[kmer] -> count
	
	## Convert DNA sequence to string of numbers.
	foreach my $c(split("",$seq)) {
		my $n	= ($c eq "a") ? 0 : 
						($c eq "c") ? 1 : 
						($c eq "g") ? 2 : 
						($c eq "t") ? 3 : 4;
		push(@bin, $n);
	}

	## Scan for 1st kmer id.  
	for (my $i = 0; $i < $k; $i++) {
		my $b = $bin[$i];
		$f 		=	($b > 3) ? 1 : 
						($f > 0) ? ($f + 1)%($k + 1) : $f;
		$kmer = $kmer * 4 + ($b % 4); ## append 0 for 4
	}
	if ($f == 0) {	$hash{$kmer}++;}	## Increment kmer count if legit.

	## Scan for 1st reverse complement kmer.
	my $word 		= $kmer;	##	A copy of the current kmer id.
	my $kmerrc 	= 0;			##	Reverse complement id.
	for (my $i = 0; $i < $k; $i++) {
		my $w = $word % 4;							##	Get next 3' base of kmer
		if ($w <= 3) { $w = 3 - $w; }		##	Complement base.
		$kmerrc = $kmerrc*4 + $w;				##	Add to kmer.
		$word = int($word/4);						##	Get rid of 3' base of kmer.
	}
	if ($f == 0 && $rc == 1) {	$hash{$kmerrc} ++;}	## Count kmerrc if legit.

	## Scan for next kmers and reverse complements.
	my $max = 4**$k;					## maximum k-mer id
	my $sig = 4**($k-1);			## most sig bit
	my $len = length($seq);		## length of sequence
	for (my $i = $k; $i < ($len - $k + 1); $i++) {
		my $b = $bin[$i];
		## check if the next char is valid or not
		$f 		=	($b > 3) ? 1 : 
						($f > 0) ? ($f + 1)%($k + 1) : $f;

		## update kmers
		$b = $b % 4; ## 0 for 4
		$kmer = ($kmer*4)%$max + $b;
		$kmerrc = (3-$b)*$sig + int($kmerrc/4);
		## update hash
		if ($f == 0) {	$hash{$kmer}++;}	## Increment kmer count if legit.
		if ($f == 0 && $rc == 1) {	$hash{$kmerrc} ++;}	## Count kmerrc if legit.
	}

	## Print kmer counts out to buffer.
	my $c 	= defined($hash{0}) ? $hash{0} : 0;
	my $buf = "$c";
	for (my $i = 1; $i < $max; $i++) {
		$c = defined($hash{$i}) ? $hash{$i} : 0;
		$buf = "$buf\t$c";
	}
	return $buf;
}

sub Base4 {
	my $num = shift;
	my @val = ();
	while ($num > 0) {
		my $ch 	= $num%4;
		$num 	 	= int($num/4);
		push @val, $ch;
	}
	my $v = join("",reverse(@val)); 
	return $v;
}

