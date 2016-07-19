#ifndef __SCORE_SEQ_H
#define __SCORE_SEQ_H

#include "sequence.h"
#include <ctype.h>
#include "tree.h"
#include <map>
#include <stdlib.h>
#include <math.h>

void  procOpts(int argc, char * argv []);
static void  usage ();
map<string,double> readModel(char* filename, int& depth);
void scoreSequences(char* filename);
void getSequenceScore(string seq, string seqname);
double scoreSequence(string seq, string seqname);
double scoreSequenceEX(string seq, string seqname);
#endif
