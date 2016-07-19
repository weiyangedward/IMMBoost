#include "tree.h"
#include "tree_node.h"
#define CHI2_PNTS 15
#define MAGIC_NUMBER	40
#include <cstdlib>

Tree::Tree(int layers, double lambda)
{
	_nlayer=layers;
	_lambda=lambda;
	_root = new TreeNode(0, layers);
}

Tree::~Tree()
{
	//recursively delete Tree starting from the root
	DeleteTree(_root);	
}

void Tree::DeleteTree(TreeNode* cNode){
	if(cNode!=NULL){
		for (int i=0; i<DEGREE; i++){
			DeleteTree(cNode->_child[i]);
		}
		delete cNode;
	}
}

void Tree::Print()
{
	//BFS
	queue<TreeNode*> q;
	q.push(_root);
	while(!q.empty()){
		TreeNode*cNode =q.front();		
		q.pop();
		if(cNode->_parent!=NULL){
			cout<<cNode->_name <<"("<< cNode->_parent->_name<<") "<<cNode->_count<<" "<<cNode->_level<<endl;
		} else{
                        cout<<cNode->_name <<"( ) "<<cNode->_count<<" "<<cNode->_level<<endl;
		}
		for (int i=0; i<DEGREE; i++){
			if(cNode->_child[i]!=NULL){
				q.push(cNode->_child[i]);
			}
		}
	}
}

void Tree::BuildTree(int* stridx, int arrSize){
// build tree from given indexed sequence
	for(int i=0; i<=arrSize-_nlayer; i++){
		//cout<<"kmer: "<<i<<endl;
		AddKmer(stridx, i);
	}
}

void Tree::AddKmer(int* stridx, int start){
// add kmer stridx[start:start+_nlayer] to the tree
	TreeNode* cNode = _root;
	char nuc[] = { 'A', 'C', 'G', 'T', 'N'};
	_root->_count++;

	for(int i=start+_nlayer-1;i>=start;i--){
		int chrIdx = stridx[i];

		if(cNode->_child[chrIdx]==NULL){
			cNode->_child[chrIdx] = new TreeNode(1, i-start);
			cNode->_child[chrIdx]->_name = nuc[chrIdx];
			cNode->_child[chrIdx]->_idx= chrIdx;
			cNode->_child[chrIdx]->_parent=cNode;
		}
		else{
			cNode->_child[chrIdx]->_count++;
		}
		cNode = cNode->_child[chrIdx];
	}
}

int Tree::getKmerCount(vector<int> kmer){ 
// returns the number of kmers from tree
        TreeNode* cNode = _root; 
        for(int i=kmer.size()-1;i>=0;i--){
                int chrIdx = kmer[i];
                if(cNode->_child[chrIdx]==NULL)	
			return 0;
                cNode = cNode->_child[chrIdx];
        }
	return cNode->_count;
}

int Tree::getKM1merCount(vector<int> kmer){
// returns the sum of all kmers that share the same prefix (i.e. (k-1)mer)
        TreeNode* cNode = _root;
	int sum=0;
        for(int i=0; i<DEGREE-1; i++){ // (ignore N)
		kmer[kmer.size()-1] = i;
                sum+=getKmerCount(kmer);
        }
        return sum;
}

double Tree::get_IMM(vector<int> kmer){
	if(kmer.size()==0) 
		return 0;
	if(kmer.size()==1)
		return getKmerProb(kmer);

        map<vector<int>,double>::iterator it;
        it=_immhash.find(kmer);
        if(it!=_immhash.end()){
                //cout<<"haha2"<<endl;
                return (double)it->second;
        }

	vector<int> vec1 = kmer;
	vec1.erase(vec1.begin());

	double retval = 0;
	vector<int> prefix = kmer;
	prefix.pop_back();
	double lambda = _lambda;

	if(_lambda<=0){ // if not a fixed lambda provided
		lambda = get_lambda(prefix);
	}

	
	if(lambda==1){
		retval= lambda*getKmerProb(kmer);
	}else{
		double timm=get_IMM(vec1);
		retval= lambda*getKmerProb(kmer)+(1-lambda)*timm;
	}
	
//	cerr<<"size "<<kmer.size()<<" lambda "<<lambda<<" prob "<< getKmerProb(kmer) << " imm "<<timm<< " retval " <<retval<<endl; 

	
	_immhash.insert(pair<vector<int>,double>(kmer,retval));
	return retval;
}

double Tree::get_lambda (vector<int> kmer){

	map<vector<int>,double>::iterator it;
	it=_lambdahash.find(kmer);
	if(it!=_lambdahash.end()){
		//cout<<"haha"<<endl;
		return (double)it->second;
	}
	
	vector <double> obs(DEGREE-1);
	double totalsum = 0.0;
	for(int i=0;i<DEGREE-1; i++){ // for A,C,G,T, (ignore N)
		vector<int> v=kmer;
		v.push_back(i);
		obs[i] = (double)getKmerCount(v);
		totalsum+=obs[i];
	}

/*	if(totalsum>0.75*MAGIC_NUMBER)
		return 0.75;
	else
		return (totalsum/MAGIC_NUMBER);
*/

	vector <double> exp(DEGREE-1);
	for(int i=0;i<DEGREE-1; i++){ // for A,C,G,T, (ignore N)
		vector<int> v=kmer;
		v.push_back(i);
		v.erase(v.begin());
		exp[i]=get_IMM(v);
	}

	for(int i=0;i<DEGREE-1; i++){
		exp[i]*=totalsum;
	}

	double chi = chi2_test(obs,exp);
	double retval=0;

/*
	if(totalsum/MAGIC_NUMBER>0.5){
		retval = 0.5*chi;
	}else{
		retval = chi*totalsum/MAGIC_NUMBER;
	}
*/
	if(totalsum/MAGIC_NUMBER>0.75){
		retval = chi*0.75;
	}else{
		retval = chi*totalsum/MAGIC_NUMBER;
	}

	_lambdahash.insert(pair<vector<int>,double>(kmer,retval));

/*	if(retval>1){
		// this will never happen
		cerr <<"lambda > 1: it has been set to 1"<<retval<<" "<<totalsum<<endl;
		retval=1;
	}
*/
	return retval;
}

double Tree::getKmerProb (vector<int> kmer){
// returns the probability of kmer given the (k-1) previous bases
	if(kmer.size()==0) return 0;
	int km1mer = getKM1merCount(kmer);
	if(km1mer==0)
		return 0;
	return (double)getKmerCount(kmer)/km1mer;
}

double Tree::chi2_test(vector<double> obs, vector<double> exp){

	double ret_val=0;
	float Chi2_Val [CHI2_PNTS] = {0.11,0.58,1.01,1.42,1.87,2.37, 2.95, 4.11, 5.1 ,6.25, 7.81, 9.35, 11.34, 12.8,16.27};
	float Chi2_Sig [CHI2_PNTS] = {0.01,0.1,0.2,0.3,0.4,0.50, 0.6 ,0.75, 0.83 ,0.90, 0.95, 0.975, 0.99, 0.995,0.999};

	double chi2=0;
	for(int i=0; i<obs.size();i++){
		if(exp[i]>0.0){
			chi2+=pow(obs[i]-exp[i],2)/exp[i];
		}
	}

	int i=0;
	for (i = 0; i < CHI2_PNTS && Chi2_Val[i] < chi2; i++);
	if (i == 0) ret_val = 0.0;
	else if (i == CHI2_PNTS) ret_val = 1.0;
	else ret_val = Chi2_Sig[i-1] + ((chi2-Chi2_Val[i-1])/(Chi2_Val[i]-Chi2_Val[i-1]))*(Chi2_Sig[i]-Chi2_Sig[i-1]);

	//cerr<<"\nrtv:"<<ret_val<<" Chi2:"<< chi2 <<endl;

	return ret_val;
}

void Tree::BuildModel(){

	//DFS
	string nuc[] = { "A", "C", "G", "T", "N"};
	stack<TreeNode*> stk;
	stk.push(_root);
	while(!stk.empty()){
		TreeNode* cNode=stk.top();
		stk.pop();
		bool nochild =0;
		for(int i=0; i<DEGREE-1; i++){ // for A,C, G, T, (ignore N)
			if(cNode->_child[i]!=NULL){
				nochild=1;
				stk.push(cNode->_child[i]);
			}
		}
		if(!nochild){
			vector<int> kmer;
			kmer.push_back(cNode->_idx);
			while(cNode->_parent!=NULL){
				cNode = cNode->_parent;
				kmer.push_back(cNode->_idx);
			}
			kmer.pop_back();
			// handle the NNNNNNN kmers
			for(int i=0; i<kmer.size();i++){
				cout <<nuc[kmer[i]];
			}
			double imm = get_IMM(kmer);			
			cout << " " <<getKmerCount(kmer)<<" "<< getKmerProb(kmer)<<" " <<imm <<endl;
		}
	}
}

void Tree::BuildFullModel(){
	string nuc[] = { "A", "C", "G", "T", "N"};
	vector<int> kmer (_nlayer,0);
	for(int i=0; i<pow(DEGREE-1, _nlayer);i++){
		toBase(i,kmer, DEGREE-1,_nlayer);
		for(int i=0; i<kmer.size();i++){
			cout <<nuc[kmer[i]];
		}
		double imm = get_IMM(kmer);			
		cout << " " <<getKmerCount(kmer)<<" "<< getKmerProb(kmer)<<" " <<imm <<endl;
	}
}

void Tree::toBase(int value, vector<int> &kmer, int base, int fixedDigitCount){
  int digitCount = 1;
  int temp = value;
  if(value< 0){
	cerr<< "No negative conversion is supported"<<endl;
	exit(0);
  }
  // Find digit count
  while ( temp /= base )
    digitCount++;
  // Compare with fixeed number of digits, use highest
  int i = max(fixedDigitCount, digitCount);
  // Convert
  for (int off;i;i--) {
    off = (value % base);
    value /= base;
    kmer[i-1] = off;
  }
}
