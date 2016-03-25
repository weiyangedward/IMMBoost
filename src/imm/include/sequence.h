// Majid Kazemian 
// Version:  1.00  15 Apr 2010
// general routines for DNA sequences.


#ifndef  __SEQUENCE_H
#define  __SEQUENCE_H

#ifndef FALSE
#define FALSE               0
#endif

#ifndef TRUE
#define TRUE                1
#endif

#include  <stdio.h>
#include  <stdlib.h>
#include  <iostream>
#include  <fstream>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <limits.h>

const long int  INCR_SIZE 			= 10000;
const long int  INIT_SIZE 			= 10000;
const int  			MAX_LINE 				= 1000;
const int  			FASTA_WIDTH 		= 60;  // Max number of characters to print on a FASTA data line
static const int MAX_SEQ_LENGTH = 40000000;

using namespace std;

class Sequence {
public:
	string _seq;
	string _name;
	long int _length;
	int* _indexedseq;

private:
	FILE* _fp;

public:
	Sequence(char *flname);
	Sequence(FILE *fpp);
	Sequence(string seq, string name);
	int  Next_Seq();
	void Reverse_Complement();
	void Fasta_Print(char *flname, char * hdr);
	void Fasta_Print_N(char *flname, int n, char * hdr);
	virtual ~Sequence();
private:
	void Index_Seq();
	char Complement(char);
	char Filter(char Ch);

};

#endif
