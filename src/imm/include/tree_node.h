#ifndef __TREENODE_H
#define __TREENODE_H

#define DEGREE	5 // (A, C, G, T, N)// in most calculations, N is ignored
#include <stdio.h>
#include <iostream>
#include <string>

using namespace std;

class TreeNode{
	public:
		TreeNode(int cnt=0, int lvl=0);
		~TreeNode();
		void PrintNode();
		TreeNode* _child[DEGREE];
		TreeNode* _parent;
		int _count;
		int _level;
		char _name; //(A, C, G, T, N)
		int _idx;     //(0, 1, 2, 3, 4)
	private:
};

#endif

