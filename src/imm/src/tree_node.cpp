#include "tree_node.h"

TreeNode::TreeNode(int cnt, int lvl)
{
	_count = cnt;
	_level = lvl;
	_name = '-';
	for(int i=0; i<DEGREE; i++){
		_child[i] = NULL;
	}
	_parent = NULL;
	_idx = -1;
}

TreeNode::~TreeNode()
{
}

void TreeNode::PrintNode()
{
	cout << "Node : "<<_name<<" "<<_count<<endl;
	if(_parent!=NULL){
		cout<<"parent: "<<_parent->_name<<" "<<_parent->_count<<endl;
	}
}

