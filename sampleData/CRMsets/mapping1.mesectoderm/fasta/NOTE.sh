less negmsCRM.fasta | grep ">" | sed 's/>//g' > negmsCRM.fasta.list
cat negmsCRM.fasta.list | perl -MList::Util=shuffle -e 'print shuffle(<STDIN>);' > negmsCRM.fasta.list.shul
head -100 negmsCRM.fasta.list.shul > negmsCRM.fasta.list.shul.100
fishInWinter.pl --bformat table --fformat fasta negmsCRM.fasta.list.shul.100 negmsCRM.fasta  > negmsCRM.fasta.100
mv negmsCRM.fasta.100 negmsCRM.fasta
