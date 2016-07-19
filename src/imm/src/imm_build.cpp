#include  <ctype.h>
#include "sequence.h"
#include "imm_build.h"
#include "tree.h"
#include "vector"
#include <getopt.h>

int depth			= 6;			// Depth of model.
bool rc 			= false;	// Use reverse complement
double lambda = 0;			// Interpolation parameter.

int main(int argc, char * argv []){
	procOpts (argc, argv);
	Sequence seq (stdin);
	Tree t(depth, lambda);
	t.BuildTree(seq._indexedseq, seq._seq.length());
	if(rc){
		seq.Reverse_Complement();
		t.BuildTree(seq._indexedseq, seq._seq.length());
	}
	while(seq.Next_Seq()){
		t.BuildTree(seq._indexedseq, seq._seq.length());
	  if(rc){
     seq.Reverse_Complement();
     t.BuildTree(seq._indexedseq, seq._seq.length());
		}
	}
	t.BuildFullModel();
	return 0;
}

//  Get options and parameters from command line.  
void  procOpts (int argc, char * argv []){
  int  errflg = FALSE;
  int  ch;
	char  * p, * q;
    optarg = NULL;
	while(! errflg && ((ch = getopt (argc, argv, "k:rl:")) != EOF)){
  	switch(ch){
			case 'k':    // model depth
      	depth = int (strtol (optarg, & p, 10));
        if (p == optarg || depth <= 0 || depth >20){
        	cerr<<"Bad model depth value" << optarg<<endl;
          errflg = TRUE;
        }	
        break;
			case 'r':   // include model reverse complement 
	  		rc = true;
	  		break;
			case 'l':
				lambda = strtod(optarg, &p);
	 	 		if(lambda <=0 || lambda>1){
					cerr<<"\nlambda should be between 0 and 1"<<endl;
					exit(0);
	  		}
	  		break;
  		case '?':
   			fprintf (stderr, "Unrecognized option -%c\n", optopt);
  		default :
    		errflg = TRUE;
		}
	}
  if (errflg || optind != argc){
	 	usage (argv[0]);
		exit (EXIT_FAILURE);
  }
	return;
}

// Print the description of options and command line for this program. 
static void  usage (char * command){
	fprintf (stderr,
		"\nUSAGE:  %s [options] \n\n"
		"Read sequences from standard input and build an Interpolated\n"
		"Markov Model (IMM) for them\n\n"
		"Options:\n"
		" -k <num>\n"
		"    Set the depth of model to <num>\n"
		" -r \n"
		"    Add reverse complement to the model\n"
		" -l <double>\n"
		"    for fixed lambda interpolation enter a value between 0 and 1\n"
		,command);
  return;
}
