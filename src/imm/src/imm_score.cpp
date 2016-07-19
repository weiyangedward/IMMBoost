#include "imm_score.h"
#include <getopt.h>
map<string, double> pos_hash;
map<string, double> neg_hash;

int pos_depth = 6;			// Depth of positive model.
int neg_depth = 6;			// Depth of negative model.
int col_idx 	= 3; 			// Column index of model file. IMM = 3 and Fixed-Model=2
bool rc 			=	false;	// Use reverse complement.
bool use_neg	=	false;	// Use negative model.
char* neg_file;					// Negative model file.

int main(int argc, char * argv []){
	procOpts (argc, argv);
	pos_hash = readModel(argv[argc-2], pos_depth);
	if(use_neg){
		neg_hash = readModel(neg_file, neg_depth);
	}
	scoreSequences((char*)argv[argc-1]);
	return 0;
}

map<string, double> readModel(char* filename, int& depth){
	FILE* fpmodel = fopen (filename, "r");
	if(fpmodel == NULL){
		cout<<"ERROR: Could not open model-file "<<filename<<endl;
    exit(1);
  }
	map<string, double> rethash;
	char* p;
	char line [MAX_LINE];

	while(fgets(line, MAX_LINE, fpmodel)!=NULL){
		int idx    = 0;
		double imm = 0;	
		string str;
		p = strtok (line, " \t\n");
		while (p != NULL){
			if (idx==0) {
				str = p;
			}
			if (idx==col_idx) {
				imm = atof(p);
			}
			p = strtok (NULL, " \t\n");
			idx++;
		}
		rethash.insert(pair<string,double>(str,imm));
		depth = str.length();
	}
	fclose(fpmodel);
	return rethash;
}

void scoreSequences(char* filename){
	Sequence seq (filename);
	getSequenceScore(seq._seq, seq._name);
	while(seq.Next_Seq()){
		getSequenceScore(seq._seq, seq._name);
	}
}

void getSequenceScore(string seq, string seqname){
	double retval = scoreSequence(seq,seqname);
	if(rc){
		Sequence tmp(seq,"tmp");
		tmp.Reverse_Complement();
		retval+= scoreSequence(tmp._seq,seqname);
	}
	cout<<seqname<<" "<<retval<<endl;
} 

// Score a single IMM Sequence
double scoreSequence(string seq, string seqname){
	double psum = 0;	
	for(int i=0; i <= seq.length()-pos_depth; i++){
		string str;
		str.assign(seq, i, pos_depth);
		map<string,double>::iterator it;
		it=pos_hash.find(str);
		double tsum =0; 
		if(it!=pos_hash.end())
			tsum=it->second;
		if(tsum>0) {
			psum+=log(tsum);
		} else {
			psum+=-2; //log(0.01)
		}
	}

	// If not using a negative model, return
	if(!use_neg) {
		return psum;
	}

	double nsum = 0;
	for(int i=0; i <= seq.length()-neg_depth; i++){
		string str;
		str.assign(seq, i, neg_depth);
		map<string,double>::iterator it;
		it = neg_hash.find(str);

		double tsum=0;
		if(it!=neg_hash.end()) {
			tsum=it->second;
		}
		if(tsum>0) {
			nsum+=log(tsum);
		} else {
			nsum+=-2; //log(0.01)
		}
	}
	return (psum-nsum);
}

double scoreSequenceEX(string seq, string seqname){
	double sum = 0;	
	for(int i = 0; i <= seq.length()-pos_depth; i++){
		string str;
		str.assign(seq, i, pos_depth);
		map<string,double>::iterator it;
		it = pos_hash.find(str);
		double psum =0; 
		if(it != pos_hash.end()) {
			psum=it->second;
		}	
		if(neg_hash.size()){ // negative model is provided
			double nsum = 0;
			it = neg_hash.find(str);
			if(it!=neg_hash.end()) {
				nsum=it->second;
			}
			if(psum && nsum){
				sum += log(psum/nsum);
			}
		} else if(psum>0) {
			sum+=log(psum);
		}
	}
	return sum;
}

//  Get options and parameters from command line.
void  procOpts (int argc, char * argv []){
	int errflg = FALSE;
	int ch;
	char* p;
	char* q;
	FILE* fp;
	optarg = NULL;
	while (! errflg && ((ch = getopt (argc, argv, "n:fr")) != EOF)){
  	switch(ch){
    	case 'n':    // read negative model file
  			use_neg		= true;
  			neg_file 	= optarg;
     		break;
			case 'f':
  			col_idx = 2;
  			break;
			case 'r':
  			rc = true;
  			break;
      case '?':
        fprintf (stderr, "Unrecognized option -%c\n", optopt);
      default :
        errflg = TRUE;
   }
 }
	if (errflg || argc < 3 || argc < optind + 2){
		usage ();
		exit (EXIT_FAILURE);
	}
	return;
}


//  Print the description of options and command line for this program.
static void  usage (){
fprintf (stderr,
	"USAGE:  Score [options] <model-file> <sequence-file>\n"
	"Reads sequences from sequence-file and scores every Fasta sequence\n"
	"using the given model. The scores are printed in standard output.\n\n"
	"Options:\n"
	" -n <filename>\n"
	"    Reads the model on negative sequences from <filename> and outputs the LL score\n"
	" -f \n"
	"    for using fixed kmer model instead of IMM use this option\n"
	" -r \n"
	"    score sequence and its reverse complement\n"
	);
	return;
}

