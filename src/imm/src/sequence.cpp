// Majid Kazemian
// Version:  1.00  15 Apr 2010
// general routines for DNA sequences.

#include "sequence.h"


Sequence::~Sequence()
{
	if(_fp!=NULL)
		fclose(_fp);
	if(_indexedseq)
		delete [] _indexedseq;
	
}


Sequence::Sequence(char *flname)
{
	_indexedseq=NULL;
	_fp = fopen(flname,"r");
	if (_fp == NULL) {
		printf("Error reading sequence file %s\n",flname);
		printf("File should be in FASTA Format\n");
		exit(1);
	}
	
	if(!Next_Seq()){
		printf("Error reading sequence file %s\n",flname);
		printf("File should be in FASTA Format\n");
		exit(1);
	}
}

Sequence::Sequence(FILE *fpp)
{
        _indexedseq=NULL;
        _fp = fpp;
        if (_fp == NULL) {
                printf("Error reading in stream\n");
                exit(1);
        }

        if(!Next_Seq()){
                printf("Error reading in stream\n");
                exit(1);
        }
}


Sequence::Sequence(string seq, string name)
{
	_seq= seq;
	_fp=NULL;
	_name=name;
	_length = seq.length();
	_indexedseq = new int[_length];
	Index_Seq();
}


char Sequence::Complement(char Ch){
	
	/* Returns the DNA complement of  Ch . */
	
	switch  (toupper(Ch)){
	case  'A' :
        return  'T';
	case  'C' :
        return  'G';
	case  'G' :
        return  'C';
	case  'T' :
        return  'A';
	case  'R' :          // a or g
        return  'Y';
	case  'Y' :          // c or t
        return  'R';
	case  'S' :          // c or g
        return  'S';
	case  'W' :          // a or t
        return  'W';
	case  'M' :          // a or c
        return  'K';
	case  'K' :          // g or t
        return  'M';
	case  'B' :          // c, g or t
        return  'V';
	case  'D' :          // a, g or t
        return  'H';
	case  'H' :          // a, c or t
        return  'D';
	case  'V' :          // a, c or g
        return  'B';
	default :            // anything
        return  'N';
	}
}


char Sequence::Filter(char Ch){
	
	//  Return a single  a, c, g or t  for  Ch .
	switch  (toupper(Ch)){
	case  'A' :
	case  'C' :
	case  'G' :
	case  'T' :
        return  Ch;
	case  'R' :     // a or g
        return  'G';
	case  'Y' :     // c or t
        return  'C';
	case  'S' :     // c or g
        return  'C';
	case  'W' :     // a or t
        return  'T';
	case  'M' :     // a or c
        return  'C';
	case  'K' :     // g or t
        return  'T';
	case  'B' :     // c, g or t
        return  'C';
	case  'D' :     // a, g or t
        return  'G';
	case  'H' :     // a, c or t
        return  'C';
	case  'V' :     // a, c or g
        return  'C';
	default :       // anything
        return  'N';
    }
}

int Sequence::Next_Seq(){
// Read next string from  fp  (assuming FASTA format) into  _seq[1 ..]
// Return  TRUE  if successful,  FALSE otherwise (e.g., EOF).
	
	string T;
	char *P, Line [MAX_LINE];
	long int  Len, Lo, Hi;
	int  Ch, Ct;
	
	while((Ch = fgetc (_fp)) != EOF && Ch != '>');
	
	if(Ch == EOF) return false;
	
	fgets (Line, MAX_LINE, _fp);
	
	Len = strlen (Line);
	assert (Len > 0 && Line [Len - 1] == '\n');
	
	P = strtok (Line, " \t\n");
	
	if(P != NULL)
		_name = P;
	else
		_name = "noname";
	
	Lo = 0;  Hi = LONG_MAX;
	
	Ct = 0;
	Len = 0;
	while((Ch = fgetc (_fp)) != EOF && Ch != '>'){
		if  (isspace (Ch)) continue;
		Ct ++;
		if(Ct < Lo || Ct > Hi) continue;
		Ch = toupper (Ch);
		switch  (Ch){
		case  'A' :
		case  'C' :
		case  'G' :
		case  'T' :
		case  'S' :
		case  'W' :
		case  'R' :
		case  'Y' :
		case  'M' :
		case  'K' :
		case  'B' :
		case  'D' :
		case  'H' :
		case  'V' :
		case  'N' :
			break;
		default :
			fprintf (stderr, "Unexpected character `%c\' in string %s\n",Ch, _name.c_str());
			Ch = 'N';
        	}
		
		T.append(1, Filter(Ch));
		Len++;
	}
	
	if  (Ch == '>') ungetc (Ch, _fp);
	
	_seq = T;
	_length = Len;
	
	if(_indexedseq!=NULL){
		delete [] _indexedseq;
	}
	_indexedseq = new int[Len+1];  
	Index_Seq();
	
	return  true;
}


void Sequence::Reverse_Complement(){
	
	//  get seq reverse complement.
	char  Ch;
	long int  i, j;

	for(i = 0, j = _length-1;  i < j;  i ++, j --){
		Ch = _seq[j];
		_seq[j] = Complement (_seq[i]);
		_seq[i] = Complement (Ch);
	}
	
	if  (i == j)
		_seq[i] = Complement (_seq[i]);

	// Indexes should be complemented too
	Index_Seq();
}


void Sequence::Fasta_Print(char *flname, char * hdr){
	
	//  Print string  s  in fasta format to  flname .  Put string  hdr
	//  on header line.
	
	FILE* fp = fopen(flname,"w");
	if (fp == NULL) {
		printf("Error writting file %s\n",flname);
		exit(1);
	}
	
	int  ct = 0;
	if(hdr != NULL)
		fprintf (fp, ">%s\n", hdr);
	else
		fprintf (fp, ">%s\n", _name.c_str());
	
	for(int i = 0;  i < _length;  i++){
		if(ct == FASTA_WIDTH){
			fputc ('\n', fp);
			ct = 0;
		}
		fputc (_seq[i], fp);
		ct ++;
	}
	fputc ('\n', fp);
	fclose(fp);
	
	return;
}


void Sequence::Fasta_Print_N(char *flname, int n, char * hdr){
	
	//  Print first  n  bytes of  string  s  in fasta format to  fp .
	//  Put string  hdr  on header line.
	
	FILE* fp = fopen(flname,"w");
	if (fp == NULL) {
		printf("Error writting file %s\n",flname);
		exit(1);
	}
	
	int  ct = 0, i;
	if(hdr != NULL)
		fprintf (fp, ">%s\n", hdr);
	else
		fprintf (fp, ">%s\n", _name.c_str());
	
	
	for(i = 0;  i < n;  i ++){
		if(ct == FASTA_WIDTH){
			fputc ('\n', fp);
			ct = 0;
		}
		fputc (_seq[i], fp);
		ct ++;
	}
	fputc ('\n', fp);
	fclose(fp);
	
	return;

}


void Sequence::Index_Seq()
{
	for (int i=0; i<_length; i++) {
		switch(_seq[i]) {
		case 'a':
		case 'A': _indexedseq[i] = 0; break;
		case 'c':
		case 'C': _indexedseq[i] = 1; break;
		case 'g': 
		case 'G': _indexedseq[i] = 2; break;
		case 't':
		case 'T': _indexedseq[i] = 3; break;
		default:  _indexedseq[i] = 4;
		} 
	}
}
