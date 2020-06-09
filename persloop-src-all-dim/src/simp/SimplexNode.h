/*
(c) 2015 Fengtao Fan, Dayu Shi
*/
#ifndef _SIMPLICIAL_TREE_NODE_H_
#define _SIMPLICIAL_TREE_NODE_H_

#include <iostream>
#include <list>
#include <unordered_map>
#include <unordered_set> 
#include <map>
#include <vector>
#include <set>
#include <queue>
#include <stack>
#include <algorithm>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/pointer_cast.hpp>
#include <boost/functional/hash.hpp>

using namespace std;

extern std::vector<std::unordered_map<int, pair<int, int>>> persistences;
extern float fThreshold;
extern int filtration_step;
extern int max_dimension;
extern vector<float> vecFiltrationScale;

/*-----Simplicial Tree---------*/
/* declaration of simpicial tree node*/
class SimplicialTreeNode;
/* shared pointer to simplicial tree node */
typedef boost::shared_ptr<SimplicialTreeNode> SimplicialTreeNode_ptr;

/*-----Annotation Matrix classes---------*/
class ListNode;
// boost shared pointer to the list node 
typedef boost::shared_ptr<ListNode> ListNodePtr;
//typedef boost::weak_ptr<ListNode> ListNodeWeakPtr;
/*------------------------------------------------*/

/*--------union find deletion classes--------*/
class TreeNode;
class TreeRootNode;

typedef boost::shared_ptr<SimplicialTreeNode> ElementNodePtr;
typedef boost::shared_ptr<TreeNode> TreeNodePtr;
typedef boost::shared_ptr<TreeRootNode> TreeRootNodePtr;
/*------------------------------------------------*/


class SimplicialTreeNode
{
public:
	// constructor
	SimplicialTreeNode() : label(-1), index_in_filtration(-1)
	{}
	SimplicialTreeNode(const int v_index) : label(v_index), index_in_filtration(-1)
	{}
	SimplicialTreeNode(const int v_index, const int idx_filtratioin) : label(v_index), index_in_filtration(idx_filtratioin)
	{}
	SimplicialTreeNode(const int v_index, TreeNodePtr in_tree_node) : label(v_index), tree_node(in_tree_node)
	{}
public:
	int label;
	int index_in_filtration;
	SimplicialTreeNode_ptr parent;
	std::unordered_map<int, SimplicialTreeNode_ptr> children;
	SimplicialTreeNode_ptr next_circular_ptr;
	SimplicialTreeNode_ptr prev_circular_ptr;
	/*By convention, a simplex with zero annotation has a nullptr tree_node */
	TreeNodePtr tree_node; // linked to the element used in the union-find-delete data struture;
	int iStatus;	//0 means the simplex is uncovered by spanning tree
};

/*************Union Find Deletion*********************************/
class TreeNode
{
public:
	TreeNode() : cListSize(0), rank(0)
	{}
public:
	ElementNodePtr elem;
	TreeNodePtr parent;
	TreeNodePtr cListHead; // children list head pointer
	TreeNodePtr cPrev; // previous pointer in its siblings list 
	TreeNodePtr cNext; // next pointer in its siblings list
	TreeNodePtr nPrev;
	TreeNodePtr nNext;
	TreeNodePtr dfsPrev;
	TreeNodePtr dfsNext;
	int cListSize;
	int rank;
};

class TreeRootNode : public TreeNode
{
public:
	TreeRootNode() : TreeNode()
	{}
public:
	TreeNodePtr nListHead;
	boost::shared_ptr<ListNode> attribute;
};


class UnionFindDeletion
{
public:
	~UnionFindDeletion()
	{
		//std::cout << "Deconstructor" << std::endl; 
	}
	TreeRootNodePtr MakeSet(ElementNodePtr &elem);
	TreeRootNodePtr Union(TreeRootNodePtr & a, TreeRootNodePtr & b);
	TreeRootNodePtr Find(TreeNodePtr curr);
	void Delete(TreeNodePtr  a);
	bool Is_singleton(TreeNodePtr a) {
		TreeRootNodePtr root = Find(a);
		return !(root->cListHead);
	}
	// 
	//public:
private:
	TreeRootNodePtr find_root(TreeNodePtr a){
		while (a->parent != a) {
			a = a->parent;
		}
		return boost::static_pointer_cast<TreeRootNode>(a);
	}
	void link_tree_node_and_element(TreeNodePtr a, ElementNodePtr b)
	{//a<------>b
		a->elem = b;
		b->tree_node = a;
		return;
	}
	void unlink_tree_node_and_element(TreeNodePtr a, ElementNodePtr b)
	{//a<---x--->b
		a->elem.reset();
		b->tree_node.reset();
		return;
	}
	void LinkSingleNodeToTreeRoot(TreeNodePtr a, TreeRootNodePtr &T);
	void insert_into_nlist(TreeNodePtr a, TreeRootNodePtr &T);
	void insert_into_clist_and_dfs_list(TreeNodePtr elem, TreeNodePtr parent);
	void Relink(TreeNodePtr a);
	void LocalRebuild(TreeNodePtr &p);
	void RemoveFromDFSList(TreeNodePtr a);
	void RemoveFromCList(TreeNodePtr a);
	void RemoveFromNList(TreeNodePtr a);
	void InsertIntoCList(TreeNodePtr p, TreeNodePtr ref, TreeNodePtr a, bool right);
	void InsertIntoNList(TreeRootNodePtr p, TreeNodePtr ref, TreeNodePtr a);
	void DeleteFromReducedTree(TreeNodePtr a);
};
/**********Annotation matrix***************************/
// Each column in the annotation matrix is a circular list with a dummy head
//// dummy head with row index -1
// such that all row elements are doublely linked 
class ListNode
{
public:
	ListNode() : row(-1), val(0)
	{// dummy head with row index -1
	}
	ListNode(const ListNode &rhs);
	ListNode& operator=(const ListNode &rhs);
	ListNode(const int r, const int v) : row(r), val(v)
	{}
	//
	void RowConnections(ListNodePtr &r, ListNodePtr &l){
		right = r;
		left = l;
		return;
	}
	void ResetRowPtrs() {
		right.reset();
		left.reset();
	}
public:
	int row; // the row 
	int val; // entry value
	ListNodePtr next; // pointer to next element in the column
	ListNodePtr right; // pointer to its right element in the row
	ListNodePtr left; // pointer to its left element int the row
};
// hash function for ListNodePtr
struct hash_ListNodePtr
{
	std::size_t operator() (const ListNodePtr &  nodePtr) const
	{
		std::size_t s = 0;
		/*ignore the dummy node [nodePtr]*/
		ListNodePtr p(nodePtr->next);
		// zero annotation has hash value 0
		while (p != nodePtr)
		{
			boost::hash_combine(s, p->row);
			boost::hash_combine(s, p->val);
			//
			p = p->next;
		}
		return s;
	}
};
struct equal_ListNodePtr {
	bool operator() (const ListNodePtr   lhs, const ListNodePtr   rhs) const {
		/*ignore the dummy node [lhs] [rhs]*/
		ListNodePtr p(lhs->next);
		ListNodePtr q(rhs->next);
		if (p == lhs && q == rhs) {
			// both are zero annotation
			return true;
		}
		if (p == lhs && q != rhs || q == rhs && p != lhs) {
			// one is zero and the other is not zero
			return false;
		}
		while (p != lhs && q != rhs) {
			if (p->row != q->row) {
				return false;
			}
			p = p->next;
			q = q->next;
		}
		return ((p == lhs) && (q == rhs));
	}
};
//
class AnnotationMatrix
{
	/*ListNodePtr == nullptr indicates the zero annotation*/
	/* each column is a list with a dummy head node */
public:
	AnnotationMatrix() : timeStamp(0)
	{}
	// only copy empty annotation matrix
	AnnotationMatrix(const AnnotationMatrix & rhs)
	{}
	~AnnotationMatrix()
	{

	}
	ListNodePtr DeepCopyAnnotationColumn(const ListNodePtr &head);
	// 
	bool search(ListNodePtr &ptr) {
		return (ann_mat.find(ptr) != ann_mat.end());
	}
	TreeRootNodePtr tree_root(ListNodePtr &ptr) {
		if (ann_mat.find(ptr) != ann_mat.end()) {
			return ann_mat[ptr];
		}
		return TreeRootNodePtr();
	}
	void update_tree_root(ListNodePtr &ptr, TreeRootNodePtr &root) {
		ann_mat[ptr] = root;
		return;
	}
	ListNodePtr make_zero_annotation() {
		ListNodePtr p(boost::make_shared<ListNode>());
		p->next = p;
		return p;
	}
	void Insert(ListNodePtr &ptr, const TreeRootNodePtr x, UnionFindDeletion &ufd);
	TreeRootNodePtr Delete(ListNodePtr &ptr);
	void clearNode(ListNodePtr &ptr, bool bUpdatePers = true);
	ListNodePtr extract_column(ListNodePtr & head) {
		if (head) {
			ListNodePtr col_head = DeepCopyAnnotationColumn(head);
			//Delete(head);
			clearNode(head, false);
			return col_head;
		}
		return ListNodePtr();
	}
	int sum_dst_with_one_column(ListNodePtr &out_dst, ListNodePtr &in_src) {
		if (search(in_src)) {
			return sum_two_annotation_with_changed_dst(out_dst, in_src);
		}
		return -1;
	}
	int genus() {
		return (int)row_ptr.size();
	}
	ListNodePtr create_cocycle(TreeRootNodePtr &root, UnionFindDeletion &ufd, bool zero_elem = false);
	int sum_two_annotation_with_changed_dst(ListNodePtr & out_dst, ListNodePtr & in_src);
	int lowest_one(ListNodePtr & head);
	void kill_cocycle_last_nonzero_bit(const int u, ListNodePtr &ext_src, UnionFindDeletion &ufd);
	bool empty() {
		return ann_mat.empty();
	}
public:
	std::unordered_map<ListNodePtr, TreeRootNodePtr, hash_ListNodePtr, equal_ListNodePtr> ann_mat;
	std::unordered_map<int, ListNodePtr> row_ptr;
	long long timeStamp;		//always stores the next coming time stamp
	int annoDim;  //dimension of the simplices this annotation matrix is for
private:
	void delete_from_doubly_linked_list(const int row, ListNodePtr & p, bool bUpdatePers = true) {
		if (p->left == p) {
			// it is the last element
			row_ptr[row].reset();
			row_ptr.erase(row);
			if (bUpdatePers && annoDim <= max_dimension)
			{
				//update persistences in dim annoDim
				unordered_map<int, pair<int, int>>::iterator iterPers = persistences[annoDim].find(row);
				(iterPers->second).second = filtration_step;
				float diff = vecFiltrationScale[filtration_step] - vecFiltrationScale[(iterPers->second).first];
				if (diff <= fThreshold && diff >= 0)
					persistences[annoDim].erase(iterPers);
			}
		}
		else {
			// not last element
			if (p == row_ptr[row]) {
				// the head 
				row_ptr[row] = p->left;
			}
			p->left->right = p->right;
			p->right->left = p->left;
		}
		//
		p->left.reset();
		p->right.reset();
		return;
	}
	void insert_into_doubly_linked_list(const int row, ListNodePtr & p) {
		std::unordered_map<int, ListNodePtr>::iterator findIter = row_ptr.find(row);
		if (findIter != row_ptr.end())
		{// insert behind the head
			p->left = findIter->second->left;
			p->right = findIter->second;
			p->left->right = p;
			findIter->second->left = p;
		}
		else
		{// the first one
			row_ptr[p->row] = p;
			p->left = p;
			p->right = p;
		}
		return;
	}
};

#endif // _SIMPLICIAL_TREE_NODE_H_
