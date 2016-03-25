#ifndef __TREE_H
#define __TREE_H

#include "tree_node.h"
#include <math.h>
#include <string.h>
#include <vector>
#include <queue>
#include <stack>
#include <map>

class Tree{
	public:
		Tree(int layers=5, double lambda=0);
		~Tree();
		void Print();
		TreeNode* _root;
		int _nlayer;
		double _lambda;
		void BuildTree(int* stridx, int arrSize);
		double getKmerProb(vector<int> kmer);
		double get_IMM(vector<int> kmer);
		void BuildModel();
		void BuildFullModel();
	private:
		void AddKmer(int* stridx, int start=0); // from stridx[start..start+_nlayer]
		void DeleteTree(TreeNode* cNode);
		int getKM1merCount(vector<int> kmer);
		int getKmerCount(vector<int> kmer);
		double get_lambda (vector<int> kmer);
		double chi2_test(vector<double> obs, vector<double> exp);
		map<vector<int>, double> _lambdahash;
		map<vector<int>, double> _immhash;
		void toBase(int value, vector <int> &kmer, int base, int fixedDigitCount = 0);
};

#endif

