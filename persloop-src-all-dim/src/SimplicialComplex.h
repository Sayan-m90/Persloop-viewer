/*
(c) 2015 Fengtao Fan, Dayu Shi
*/
#ifndef _SIMPLICIAL_TREE_H_
#define _SIMPLICIAL_TREE_H_

#include "SimplexNode.h" 

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <cstring>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <queue>
#include <ctime>


//template<class T, class Allocator = std::allocator<T>>
//class SimplicialTree;

typedef boost::shared_ptr<AnnotationMatrix> AnnotationMatrixPtr;
typedef boost::shared_ptr<std::unordered_map<int, SimplicialTreeNode_ptr> > LabelsDictionaryPtr;

extern vector<float> vecFiltrationScale;
extern int time_in_each_filtration_step;
extern std::clock_t start, timer1, timer2;
extern double dFuncTimeSum;
extern double dInsertTime;
extern double dCollapseTime;
// extern TreeLoopTracker_ptr tlt_ptr;
// extern TreeLoopTracker_ptr tlt_ptr_filt;

class SimpleGraph 
{
public:
	void addEdge(int a, int b)
	{
		std::vector<int>& nbs = adj[a];
		for (int nb : nbs)
			if (nb == b) return;
		nbs.push_back(b);
		adj[b].push_back(a);
		maxIndex = std::max(maxIndex, a);
		maxIndex = std::max(maxIndex, b);
		return;
	}
	int getNumE() { return this->numE; }
	bool findInVec(vector<int>& vec, int num) {
		for (int n : vec)
			if (n == num) return true;
		return false;
	}
	//void printGraph();
	std::vector<int> BFS(int src, int dest)
	{
		std::vector<int> vecRe;
		std::vector<int> parent(maxIndex + 1, -1);
		parent[src] = src;
		parent[dest] = dest;
		vector<int> curSet = { src }, nextSet = { dest };
		while (!curSet.empty()) {
			vector<int> newSet;
			for (int nd : curSet) {
				for (int nb : adj[nd]) {
					if (findInVec(nextSet, nb)) {
						//recover path and return
						int p = nb;
						while (p != dest && p != src) {
							vecRe.push_back(p);
							p = parent[p];
						}
						if (p == src) 
							vecRe.push_back(src);
						else 
							vecRe.push_back(dest);
						std::reverse(vecRe.begin(), vecRe.end());
						p = nd;
						while (p != dest && p != src) {
							vecRe.push_back(p);
							p = parent[p];
						}
						if (p == src)
							vecRe.push_back(src);
						else
							vecRe.push_back(dest);
						return vecRe;
					}
					if (parent[nb] == -1) {
						parent[nb] = nd;
						newSet.push_back(nb);
					}
				}
			}
			curSet = nextSet;
			nextSet = newSet;
		}
		return vecRe;
	}
private:
	std::unordered_map<int, std::vector<int>> adj;
	int numE;
	int maxIndex;
};

/* declaration of simpicial tree*/
template<class T>
class SimplicialTree
{
public:
	/*constructor */
	// default constructor
	SimplicialTree() : dim(-1), EuclideanDataPtr(NULL), accumulativeSimplexSize(0), bGenerator(false)
	{
		vecTS.resize(max_dimension + 5);
		//counts.resize(100);
	}
	//
	void UpdateAnnotationArray(const int simplex_dim);
	ListNodePtr find_annotation(SimplicialTreeNode_ptr p);
	SimplicialTreeNode_ptr InsertSimplexWithAnnotation(std::vector<int> &simplex_vertices, ListNodePtr &simplex_ann);
	SimplicialTreeNode_ptr ElementaryInsersion(std::vector<int> &simplex_vertices);
	// SimplicialTreeNode_ptr ElementaryInsersionPostShortLoop(std::vector<int> &simplex_vertices);
	std::vector<int> AddTwoCycles(std::vector<int>, std::vector<int>);
	/*deconstructor */
	~SimplicialTree()
	{}
	// size statistics of the complex
	int SimplexDim(const SimplicialTreeNode_ptr simplex);
	int ComplexSize(int dim = max_dimension);
	int dDimSimplicesSize(const int d);
	void Report();
	// 
	void retrieve_vertex_indices(const SimplicialTreeNode_ptr &simplex, vector<int> &simplex_vertices);
	SimplicialTreeNode_ptr find(vector<int> &simplex_vertices);
	void insert_into_circular_list(const int simplex_dim, const SimplicialTreeNode_ptr & simplex);
	SimplicialTreeNode_ptr insert_into_simplicial_tree(std::vector<int> &simplex_vertices);
	// bool CheckBoundary (std::vector<int> simplex_vertices);
	// retrieve boundary
	bool Boundary(const SimplicialTreeNode_ptr& sigma, std::vector<SimplicialTreeNode_ptr> &bdries);
	// retrieve coboundary
	bool is_link_condition_satisfied(SimplicialTreeNode_ptr &edge, boost::unordered_set<SimplicialTreeNode_ptr> &intersectedLinks);
	bool is_upto_p_link_condition_satisfied(SimplicialTreeNode_ptr &edge, boost::unordered_set<SimplicialTreeNode_ptr> &intersectedLinks, int p);
	bool is_coface(SimplicialTreeNode_ptr &sigma, vector<int> &simplex_vertices);
	bool is_coface_set(SimplicialTreeNode_ptr &sigma, vector<int> &simplex_vertices, unordered_set<int>& hset);
	bool is_coface1(SimplicialTreeNode_ptr &sigma, vector<int> &simplex_vertices);
	void retrieve_simplices_from_subtrees(const SimplicialTreeNode_ptr &sigma, vector<SimplicialTreeNode_ptr> &coBdries, bool dfs_visit = false);
	void retrieve_simplices_from_subtrees_set(const SimplicialTreeNode_ptr &sigma, vector<SimplicialTreeNode_ptr> &coBdries, unordered_set<int>& hset, int p = INT_MAX);
	void CoDimensionOneFaces(const SimplicialTreeNode_ptr& sigma, std::vector<SimplicialTreeNode_ptr> &codim_one_faces);
	bool CoBoundary(const SimplicialTreeNode_ptr &sigma, std::vector<SimplicialTreeNode_ptr> &coBdries, bool dfs_visit = false);
	void CoBoundaryWithRestricSet(const SimplicialTreeNode_ptr &sigma, std::vector<SimplicialTreeNode_ptr> &coBdries, vector<int>& commonNB, int p = INT_MAX, bool dfs_visit = false);
	void LinkSubcomplex(const SimplicialTreeNode_ptr &sigma, std::vector<SimplicialTreeNode_ptr> &simplexLink);
	void LinkSubcomplexSet(const SimplicialTreeNode_ptr &sigma, std::vector<SimplicialTreeNode_ptr> &simplexLink, vector<int>& commonNB, int p = INT_MAX);
	void LinsubcomplexIntersection(SimplicialTreeNode_ptr & a, SimplicialTreeNode_ptr & b, boost::unordered_set<SimplicialTreeNode_ptr> &intersectedLinks, int p = INT_MAX);
	void FindCommonNeighbors(int a, int b, vector<int>& commonNB);
	// cone over the subcomplexes
	void AddExtraSimplicesToSatisfyLinkCondition(SimplicialTreeNode_ptr &edge, boost::unordered_set<SimplicialTreeNode_ptr> &intersectedLinks);
	void EdgeContraction(SimplicialTreeNode_ptr &edge, const int preserve_vertex_label);
	void remove_simplex_from_both_complex_and_ufdForest(SimplicialTreeNode_ptr sigma, bool bUpdatePers = true);
	void remove_simplex_from_complex(SimplicialTreeNode_ptr sigma);
	void rename_simplices_in_subtree(SimplicialTreeNode_ptr sigma, const int new_label);
	// annotation transfer 
	void AnnotationTransfer(SimplicialTreeNode_ptr &edge, const int preserve_vertex_label);
	// construct the complex
	//virtual bool Construction() = 0; 
	void ElementaryCollapse(int remove_label,  int preserve_label);
	/* ------------- simplicial complex I/O ---------------- */
	string annotation_to_string(SimplicialTreeNode_ptr sigma);
	ListNodePtr string_to_annotation(string &line);
	void ReadComplex(const char* pFileName);
	void ReadComplexWithAnnotation(const char* pFileName);
	void ReadSimplicialMap(const char* pFileName, vector<pair<int, int> > &vertex_map);
	void WriteComplex(const char* pFileName);
	void WriteComplexWithAnnotation(const char* pFileName);
	void PrintComplexWithAnnotation();
	//
	void WriteStatisticsToFile(const char* pFileName);
	//
	/* ----------------persistence related functions ---------------*/
	void SnapshotHomologicalFeatures(vector<std::unordered_set<int> > &hom_info);
	void PerformSimplicialCollapse(vector<pair<int, int> > &vertex_map, unordered_map<int, int> &updated_vertex_map);
	void InitializeByRenamingIncomingComplex(SimplicialTree<T> &src, unordered_map<int, int> &vertex_map);
	//relabeling instead of copying
	void RelabelingVertices(std::unordered_map<int, int> &vertex_map);
	void AddRemainingSimpliciesFromFile(const char* pFileName);
	void CheckPersistence(vector<std::unordered_set<int> > &homo_info);
	//

public:
	T* EuclideanDataPtr; // data used for constructing this complex, like the underlying graph for Rips complex
	int dim; // dimension of the simpicial complex  
	std::vector<int> simplex_sizes;
	//
	std::vector<LabelsDictionaryPtr> labels_dict_in_each_dim; // circular list for each 
	// union-find-deletion data structure
	UnionFindDeletion ufd;
	// Annotation in each dimension
	std::vector<AnnotationMatrixPtr> annotations;
	//max index of image vertices in vertex map
	int maxImageVertex;
	//time stamps in each dimension
	std::vector<long long> vecTS;
	//accumulative simplex size
	int accumulativeSimplexSize;

	//spanning tree
	SimpleGraph spanTree;

	bool bGenerator;
	//generators for dim 1 <timeStamp, generating cycle vertices>
	std::unordered_map<int, std::vector<int>> gen1;

	std::unordered_map<int, int> reindex;

	//std::vector<long long> counts;

	int findRoot(int a) {
		int parent = reindex[a];
		if (parent == a) return a;
		else return reindex[a] = findRoot(parent);
	}

	//clear function
	void clearData()
	{
		this->dim = 0;
		this->simplex_sizes.clear();
		this->labels_dict_in_each_dim.clear();
		//ufd
		this->annotations.clear();
		this->vecTS.resize(max_dimension + 3);
	}
	void clearMemory();		//clear memory of annotation matrices, simplicial trees and udf trees

	bool check_parent_children(SimplicialTreeNode_ptr sigma) {
		for (std::unordered_map<int, SimplicialTreeNode_ptr>::iterator mIter = sigma->children.begin();
			mIter != sigma->children.end(); ++mIter) {
			if (mIter->first != mIter->second->label) {
				cout << "label is not consistent " << endl;
				return false;
			}
			if (mIter->second->parent != sigma) {
				cout << "parent is not consistent" << endl;
				return false;
			}
		}
		//int simplex_dim = SimplexDim(sigma);
		//if(annotations[simplex_dim]->ann_mat.find(find_annotation(sigma)) == annotations[simplex_dim]->ann_mat.end()) {
		//	cout << "sigma annotation is not found" << endl;
		//	return false;
		//}
		return true;
	}
	bool check_ufd_elem(SimplicialTreeNode_ptr sigma) {
		if (!sigma->tree_node) {
			cout << "empty tree_node" << endl;
			return false;
		}
		else {
			if (sigma->tree_node->elem != sigma) {
				while (sigma) {
					cout << sigma->label << " ";
					sigma = sigma->parent;
				}
				cout << endl;
				cout << "tree node is not consistent" << endl;
				return false;
			}
		}
		return true;
	}
	bool check_ufd_anno() {
		for (int i = 0; i < annotations.size(); ++i) {
			std::unordered_map<ListNodePtr, TreeRootNodePtr, hash_ListNodePtr, equal_ListNodePtr>::iterator  mIter = annotations[i]->ann_mat.begin();
			for (; mIter != annotations[i]->ann_mat.end(); ++mIter) {
				if (mIter->second->attribute != mIter->first) {
					cout << "ufd and anno is not consistent" << endl;
					return false;
				}
			}
		}
		return true;
	}
	int check_status() {
		// check children and parent
		// check tree_node and elem
		// check root and annotation
		cout << "check status" << endl;
		check_ufd_anno();
		vector<int> simplex_vertices;
		simplex_vertices.reserve(dim + 1);
		// visit each simplex through the labels_diction variables
		for (int i = 0; i < labels_dict_in_each_dim.size(); ++i) {

			std::unordered_map<int, SimplicialTreeNode_ptr> & labels_dict = *(labels_dict_in_each_dim[i]);
			//
			std::unordered_map<int, SimplicialTreeNode_ptr>::iterator mIter = labels_dict.begin();
			if (i == 0) {
				// handle vertices
				for (; mIter != labels_dict.end(); ++mIter) {
					//
					if (!check_parent_children(mIter->second) || !check_ufd_elem(mIter->second)) {
						exit(0);
					}
				}
			}
			else {
				for (; mIter != labels_dict.end(); ++mIter) {
					SimplicialTreeNode_ptr trav = mIter->second;
					do{
						if (!check_parent_children(trav) || !check_ufd_elem(trav)) {
							exit(0);
						}
						//
						trav = trav->next_circular_ptr;
					} while (trav != mIter->second);
				}
			}
		}
		return 1;

	}
private:
	struct SimplicialTreeNodePtrLessThan {
		bool operator()(const SimplicialTreeNode_ptr lhs, const SimplicialTreeNode_ptr rhs) const {
			int lhs_size = 0, rhs_size = 0;
			SimplicialTreeNode_ptr trav = lhs;
			while (trav) {
				++lhs_size;
				trav = trav->parent;
			}
			trav = rhs;
			while (trav) {
				++rhs_size;
				trav = trav->parent;
			}
			return lhs_size < rhs_size;
		}
	};
	void merge_two_sorted_arrays(vector<int> &A, vector<int> &B, vector<int> &res) {
		res.resize(A.size() + B.size());
		int a = (int)A.size() - 1, b = (int)B.size() - 1, c = (int)res.size() - 1;
		while (a >= 0 && b >= 0) {
			res[c--] = max(A[a], B[b]);
			A[a] > B[b] ? --a : --b;
		}
		while (a >= 0) {
			res[c--] = A[a--];
		}
		while (b >= 0) {
			res[c--] = B[b--];
		}
		return;
	}
	void delete_from_circular_list(SimplicialTreeNode_ptr &sigma, const int simplex_dim) {
		std::unordered_map<int, SimplicialTreeNode_ptr> &head = (*labels_dict_in_each_dim[simplex_dim]);
		if (simplex_dim == 0 || sigma->next_circular_ptr == sigma) {
			// the last label in the this row 
			head.erase(sigma->label);
		}
		else {
			SimplicialTreeNode_ptr prev = sigma->prev_circular_ptr;
			SimplicialTreeNode_ptr next = sigma->next_circular_ptr;
			/*while (prev->next_circular_ptr != sigma) {
			prev = prev->next_circular_ptr;
			}*/
			prev->next_circular_ptr = next;
			next->prev_circular_ptr = prev;
			if (head[sigma->label] == sigma) {
				// 
				head[sigma->label] = sigma->next_circular_ptr;
			}
		}
		sigma->next_circular_ptr.reset();
		sigma->prev_circular_ptr.reset();
	}
	//
	bool is_zero_annotation(ListNodePtr &ptr) {
		if (!ptr) {
			return false; // empty 
		}
		return (ptr->next == ptr);
	}
};

template<typename T>
inline int SimplicialTree<T>::SimplexDim(const SimplicialTreeNode_ptr simplex) {
	int simplex_dim = -1;
	SimplicialTreeNode_ptr trav = simplex;
	while (trav) {
		++simplex_dim;
		trav = trav->parent;
	}
	return simplex_dim;
}
template<typename T>
inline int SimplicialTree<T>::ComplexSize(int dim)
{
	int total_size = 0;
	for (int i = 0; i < simplex_sizes.size() && i <= dim; i++){
		total_size += simplex_sizes[i];
	}
	return total_size;
}

template<typename T>
inline int SimplicialTree<T>::dDimSimplicesSize(const int d)
{
	if (d > dim || d < 0)
	{
		std::cout << "Querying dimension exceeding the Complex Dimension " << std::endl;
		return  -1;
	}
	return simplex_sizes[d];
}

template<typename T>
inline void SimplicialTree<T>::Report()
{
	for (int i = 0; i < dim + 1; i++)
	{
		std::cout << "DIM " << i << " : " << simplex_sizes[i] << std::endl;
	}
	std::cout << "total size : " << ComplexSize() << std::endl;
	return;
}
template<typename T>
ListNodePtr SimplicialTree<T>::find_annotation(SimplicialTreeNode_ptr p){
	if (p->tree_node) {
		TreeRootNodePtr root = ufd.Find(p->tree_node);
		return root->attribute;
	}
	else {
		if (SimplexDim(p) <= max_dimension)
		{
			cout << "Simplex has no annotation available" << endl;
			exit(0);
		}
		else
			return NULL;
	}
	return ListNodePtr();
}
template<typename T>
ListNodePtr SimplicialTree<T>::string_to_annotation(string &line) {
	int start = 0;
	while (start < line.size() && line[start] != '[') {
		++start;
	}
	if (start == line.size()) {
		// invalid annoation
		return ListNodePtr();
	}
	int nonzero_bits = 0;
	// make dummy head
	ListNodePtr dummy_head = boost::make_shared<ListNode>();
	ListNodePtr trav = dummy_head;
	//
	stringstream sstr(stringstream::in);
	//
	int row = 0, val = 0;
	for (int i = start; i < line.size(); ++i) {
		if (line[i] == '<') {
			start = i;
		}
		else {
			if (line[i] == '>') {

				sstr.str(string(line.begin() + start + 1, line.begin() + i));
				//cout << sstr.str() << endl;
				sstr >> row >> val;
				//
				trav->next = boost::make_shared<ListNode>(row, val);
				trav = trav->next;
				//
				sstr.str("");
				sstr.clear();
			}
		}
	}
	trav->next = dummy_head;
	//
	return dummy_head;
}
template<typename T>
string SimplicialTree<T>::annotation_to_string(SimplicialTreeNode_ptr sigma) {
	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	ListNodePtr ann = find_annotation(sigma);
	if (!ann) {
		/*cout << "annotation is not valid " << endl;
		exit(0);*/
		return "";
	}
	ListNodePtr trav(ann->next);
	sstr << "[";
	while (trav != ann) {
		sstr << "<" << trav->row << " " << 1 << ">" << (trav->next == ann ? "" : " ");
		trav = trav->next;
	}
	sstr << "]";
	return sstr.str();
}
template<typename T>
void SimplicialTree<T>::retrieve_vertex_indices(const SimplicialTreeNode_ptr &simplex, vector<int> &simplex_vertices) {
	if (simplex) {
		SimplicialTreeNode_ptr trav = simplex;
		while (trav) {
			simplex_vertices.push_back(trav->label);
			trav = trav->parent;
		}
		reverse(simplex_vertices.begin(), simplex_vertices.end());
	}
	return;
}
template<typename T>
SimplicialTreeNode_ptr SimplicialTree<T>::find(vector<int> &simplex_vertices){
	// preq: numbers in simplex_vertices are sorted either in increasing order
	if (simplex_vertices.empty() || labels_dict_in_each_dim.empty()) {
		// no vertices in input
		// or the complex is empty 
		return SimplicialTreeNode_ptr();
	}
	int label = simplex_vertices.front();
	LabelsDictionaryPtr vertices = labels_dict_in_each_dim.front();
	SimplicialTreeNode_ptr pIter = (vertices->find(label) == vertices->end() ? SimplicialTreeNode_ptr() : (*vertices)[label]);
	if (pIter) {
		for (int i = 1; i < simplex_vertices.size(); ++i) {
			label = simplex_vertices[i];
			std::unordered_map<int, SimplicialTreeNode_ptr>::iterator pNext;
			if (pIter->children.empty())
				return SimplicialTreeNode_ptr();
			else
			{
				pNext = pIter->children.find(label);
				if (pNext == pIter->children.end())
					return SimplicialTreeNode_ptr();
			}
			pIter = pNext->second;
		}
	}
	return pIter;
}
template<typename T>
SimplicialTreeNode_ptr SimplicialTree<T>::insert_into_simplicial_tree(std::vector<int> &simplex_vertices) {
	// preq:	1) simplex_vertices is sorted
	//			2) all boundaries are inserted already
	//			3) this is a new simplex
	const int simplex_dim = (int)simplex_vertices.size() - 1;
	//
	if (simplex_sizes.size() <= simplex_dim) {
		simplex_sizes.push_back(1);
		// creater the dictonary (hash) when insert the first element
		labels_dict_in_each_dim.push_back(boost::make_shared<unordered_map<int, SimplicialTreeNode_ptr> >());
		++dim;
	}
	else {
		++simplex_sizes[simplex_dim];
	}
	//
	unordered_map<int, SimplicialTreeNode_ptr> & vertices = (*labels_dict_in_each_dim.front());
	//
	if (simplex_dim == 0) {
		int ver_num = simplex_vertices[0];
		// it is a new vertex
		//SimplicialTreeNode_ptr simplex(boost::make_shared<SimplicialTreeNode>(ver_num, ver_num));
		SimplicialTreeNode_ptr simplex(boost::make_shared<SimplicialTreeNode>(ver_num, ++time_in_each_filtration_step));
		vertices[ver_num] = simplex;
		//
		return vertices[ver_num];
	}
	SimplicialTreeNode_ptr pIter = vertices[simplex_vertices[0]];
	int i;
	for (i = 1; i < simplex_vertices.size() - 1 && pIter != NULL; ++i) {
		pIter = pIter->children[simplex_vertices[i]];
	}
	if (pIter == NULL)
	{
		cerr << "one of the sub-simplecies of";
		for (int k = 0; k < simplex_vertices.size(); ++k)
		{
			cerr << " " << simplex_vertices[k];
		}
		cerr << " is missing.\n";
		exit(0);
	}
	if (pIter->children.find(simplex_vertices.back()) == pIter->children.end()) {
		// this simplex is not existed
		SimplicialTreeNode_ptr simplex(boost::make_shared<SimplicialTreeNode>(simplex_vertices.back(), ++time_in_each_filtration_step));
		pIter->children[simplex_vertices.back()] = simplex;
		// 
		simplex->parent = pIter;
		//
		insert_into_circular_list(simplex_dim, simplex);
		//
	}
	//
	int sz = pIter->children.size();
	//counts[sz / 10]++;
	return pIter->children[simplex_vertices.back()];
}
template<typename T>
void SimplicialTree<T>::UpdateAnnotationArray(const int simplex_dim) {
	if (annotations.size() < simplex_dim + 1) {
		annotations.push_back(boost::make_shared<AnnotationMatrix>());
		annotations.back()->timeStamp = this->vecTS[simplex_dim];
		annotations.back()->annoDim = simplex_dim;
	}
	return;
}

template<typename T>
SimplicialTreeNode_ptr SimplicialTree<T>::InsertSimplexWithAnnotation(std::vector<int> &simplex_vertices, ListNodePtr &simplex_ann) {
	SimplicialTreeNode_ptr simplex = find(simplex_vertices);
	if (!simplex) {
		// new simplex
		// 1) insert it into the simplicial tree
		simplex = insert_into_simplicial_tree(simplex_vertices);
		//
		int simplex_dim = (int)simplex_vertices.size() - 1;
		UpdateAnnotationArray(simplex_dim);
		//
		TreeRootNodePtr root = ufd.MakeSet(simplex);
		annotations[simplex_dim]->Insert(simplex_ann, root, ufd);
		//  
		//spanning tree...
	}
	return simplex;
}

template<typename T>
std::vector<int> SimplicialTree<T>::AddTwoCycles(std::vector<int> v1, std::vector<int> v2)
{
	std::vector<int> vecRe;
	std::unordered_set<int> setV1;
	std::unordered_set<int> setCommon;
	int a1 = -1, a2 = -1;
	for (int i = 0; i < v1.size(); ++i)
	{
		setV1.insert(v1[i]);
	}
	//int iSize2 = v2.size();
	for (int i = 0; i < v2.size(); ++i)
	{
		if (setV1.find(v2[i]) != setV1.end())
		{
			setCommon.insert(v2[i]);
			if (a1 == -1)	//normal case
			{
				a1 = i;
			}
			else if (a1 == 0 && a2 != -1)	//common part is broken into two parts, begin and end
			{
				a1 = i;
			}
		}
		else
		{
			if (a1 != -1 && a2 == -1)
				a2 = (i - 1) % v2.size();
		}
	}
	if (a2 == -1)
	{
		if (a1 == 0)
			return vecRe;
		else
			a2 = v2.size() - 1;
	}
	int b1, b2;
	for (int i = 0; i < v1.size(); ++i)
	{
		if (v1[i] == v2[a1])
			b1 = i;
		if (v1[i] == v2[a2])
			b2 = i;
	}
	int iCommon = (a2 - a1 + 1 + v2.size()) % v2.size();
	//a2 to a1
	for (int i = 0; i < v2.size() - iCommon + 2; ++i)
	{
		int index = (a2 + i) % v2.size();
		vecRe.push_back(v2[index]);
	}
	//b1 to b2
	for (int i = 0; i < v1.size() - iCommon; ++i)
	{
		int index = (b1 - i - 1 + v1.size()) % v1.size();
		vecRe.push_back(v1[index]);
	}
	return vecRe;
}


// Find annotation of each cycle in the given filtration
// Do reduction in the matrix form, see which ones survive


template<typename T>
SimplicialTreeNode_ptr SimplicialTree<T>::ElementaryInsersion(std::vector<int> &simplex_vertices) {
	int a, b;
	if (simplex_vertices.size() == 2) {
		a = simplex_vertices[0];
		b = simplex_vertices[1];
	}
	//reindex simplex
	if (simplex_vertices.size() == 1 && reindex.find(simplex_vertices[0]) == reindex.end())
		reindex[simplex_vertices[0]] = simplex_vertices[0];
	else {
		for (int i = 0; i < simplex_vertices.size(); i++) {
			simplex_vertices[i] = findRoot(simplex_vertices[i]);
		}
		std::sort(simplex_vertices.begin(), simplex_vertices.end());
	}

	// Preq: its faces are inserted
	// ensure that the annotations in required dimensions are allocated
	if (simplex_vertices.size() > max_dimension + 2)
		return SimplicialTreeNode_ptr();
	SimplicialTreeNode_ptr simplex = find(simplex_vertices);

	if (!simplex) {
		// new simplex
		this->accumulativeSimplexSize += 1;
		// 1) insert it into the simplicial tree
		simplex = insert_into_simplicial_tree(simplex_vertices);
		//
		int simplex_dim = (int)simplex_vertices.size() - 1;
		if (simplex_dim <= max_dimension && annotations.size() < simplex_dim + 1)
			UpdateAnnotationArray(simplex_dim);
		//
		if (simplex_dim == 0) {

			// put it into the union find data struture;
			TreeRootNodePtr root = ufd.MakeSet(simplex); // connect simplex with tree_node
			ListNodePtr newCycleAnno;
			newCycleAnno = annotations[0]->create_cocycle(root, ufd); // connect tree_root_node with annotation
			//update persistences in dim simplex_dim
			int newTS = annotations[simplex_dim]->lowest_one(newCycleAnno);
			if (persistences.size() < simplex_dim + 1)
				persistences.resize(simplex_dim + 1);
			persistences[simplex_dim].insert(std::make_pair(newTS, std::make_pair(filtration_step, -1)));
			//   
			//check_status();
		}
		else {
			// 2) take the boundary of this simplex;
			vector<SimplicialTreeNode_ptr> boundaries;
			Boundary(simplex, boundaries);

			// 3) check the annotation sum of these boundary simplices
			//find_annotation returns listnode pointer
			// cout<<"Before deep copy\n";
			ListNodePtr sum = annotations[simplex_dim - 1]->DeepCopyAnnotationColumn(find_annotation(boundaries.front()));
			int dead_bit = -1;
			for (int i = 1; i < boundaries.size(); ++i){
				if (boundaries[i] == NULL)
				{
					cerr << "one of the sub-simplecies of";
					for (int k = 0; k < simplex_vertices.size(); ++k)
					{
						cerr << " " << simplex_vertices[k];
					}
					cerr << " is missing.\n";
					exit(0);
				}
				ListNodePtr tempLNP = find_annotation(boundaries[i]);
				dead_bit = annotations[simplex_dim - 1]->sum_two_annotation_with_changed_dst(sum, tempLNP);
			}

			//empty annotation
			if (dead_bit == -1) {
				// 4.1) the sum of boundary annotations is zero
				//		create a new cocyle
				if (simplex_dim <= max_dimension)
				{
					TreeRootNodePtr root = ufd.MakeSet(simplex); // connect simplex with tree_node
					ListNodePtr newCycleAnno;
					newCycleAnno = annotations[simplex_dim]->create_cocycle(root, ufd, annotations[simplex_dim]->empty()); // connect tree_root_node with annotation
					//update persistences in dim simplex_dim
					int newTS = annotations[simplex_dim]->lowest_one(newCycleAnno);
					if (persistences.size() < simplex_dim + 1)
						persistences.resize(simplex_dim + 1);
					persistences[simplex_dim].insert(std::make_pair(newTS, std::make_pair(filtration_step, -1)));
					//find generator for this persistence homology class
					if (bGenerator && simplex_dim == 1 )
						this->gen1[newTS] = spanTree.BFS(a, b);
				}

			}
			else {
				// 4.2) the sum of boundary annotations is not zero
				//		kill the cocycle represented by the last nonzero bit  

				annotations[simplex_dim - 1]->kill_cocycle_last_nonzero_bit(dead_bit, sum, ufd);
				//update persistences in dim simplex_dim - 1
				std::unordered_map<int, pair<int, int>>::iterator itPer = persistences[simplex_dim - 1].find(dead_bit);

				(itPer->second).second = filtration_step;
				float diff = vecFiltrationScale[filtration_step] - vecFiltrationScale[(itPer->second).first];
				if (diff <= fThreshold && diff >= 0)
					persistences[simplex_dim - 1].erase(itPer);
				if (bGenerator) {
					//add an edge in spanning tree
					if (simplex_dim == 1)
						this->spanTree.addEdge(a, b);
				}
				// assign zero annotation to this simplex
				if (simplex_dim <= max_dimension)
				{
					TreeRootNodePtr root = ufd.MakeSet(simplex);
					ListNodePtr zero_ann = annotations[simplex_dim]->make_zero_annotation();
					annotations[simplex_dim]->Insert(zero_ann, root, ufd);
				}
			}
			//check_status();
			//break the cycle of "sum" annotation for deletion
			sum->next.reset();
		}
	}
	return simplex;
}

template <class T>
bool SimplicialTree<T>::Boundary(const SimplicialTreeNode_ptr& sigma, std::vector<SimplicialTreeNode_ptr> &bdries)
{
	// traverse back to the root to get all boundary faces
	SimplicialTreeNode_ptr pIter(sigma);
	SimplicialTreeNode_ptr pBdryFace;
	//
	unordered_map<int, SimplicialTreeNode_ptr> & vertices = (*labels_dict_in_each_dim.front());
	//
	// vertex simplex doesn't have boundary
	// the half vertex set
	std::vector<int> tail_ver_index_set;
	if (pIter->parent)
	{// simplex of dim > 0
		do
		{// the boundary face without current vertex 
			// go from its parent to visit all tail vertex in tail_ver_index_set
			// at the end of visit, the boundary face pointer is reached.
			if (pIter->parent)
			{
				pBdryFace = pIter->parent;
			}
			else
			{// at the moment of removing the first vertex
				// require tail_ver_index_set NOT empty
				pBdryFace = vertices[tail_ver_index_set[0]];
				tail_ver_index_set.erase(tail_ver_index_set.begin());
			}
			//
			for (int child : tail_ver_index_set)
			{//
				pBdryFace = (pBdryFace->children)[child];
			}
			//
			bdries.push_back(pBdryFace);
			//
			if (pIter->parent == NULL) break;
			tail_ver_index_set.insert(tail_ver_index_set.begin(), pIter->label);
			// traverse back one more step
			pIter = pIter->parent;
		} while (pIter);
	}
	/******************************************/
	return true;
}
template <class T>
void SimplicialTree<T>::insert_into_circular_list(const int simplex_dim, const SimplicialTreeNode_ptr & simplex)
{
	if (simplex_dim > 0) {
		unordered_map<int, SimplicialTreeNode_ptr> & dict = (*labels_dict_in_each_dim[simplex_dim]);
		if (dict.find(simplex->label) == dict.end()) {
			// it is first one
			dict[simplex->label] = simplex;
			simplex->next_circular_ptr = simplex;
			simplex->prev_circular_ptr = simplex;
		}
		else {
			// insert just behind the head element
			simplex->next_circular_ptr = dict[simplex->label]->next_circular_ptr;
			dict[simplex->label]->next_circular_ptr = simplex;
			simplex->prev_circular_ptr = dict[simplex->label];
			simplex->next_circular_ptr->prev_circular_ptr = simplex;
		}
	}
	// no need to handle the vertices
	return;
}

template <typename T>
void SimplicialTree<T>::retrieve_simplices_from_subtrees_set(const SimplicialTreeNode_ptr &sigma, vector<SimplicialTreeNode_ptr> &coBdries, unordered_set<int>& hset, int p) {
	queue<SimplicialTreeNode_ptr> Q;
	int curDim = SimplexDim(sigma);
	Q.push(sigma);
	while (!Q.empty() && curDim < p + 1) {
		int num = Q.size();
		while (num--) {
			SimplicialTreeNode_ptr curr = Q.front();
			Q.pop();
			coBdries.push_back(curr);
			for (auto cd : curr->children) {
				if (hset.find(cd.first) != hset.end())
					Q.push(cd.second);
			}
		}
		curDim++;
	}
}

template <typename T>
void SimplicialTree<T>::retrieve_simplices_from_subtrees(const SimplicialTreeNode_ptr &sigma,
	vector<SimplicialTreeNode_ptr> &coBdries, bool dfs_visit) {
	//
	if (dfs_visit) {
		// do a bfs visit to get all cofaces
		boost::unordered_set<SimplicialTreeNode_ptr> visited;
		stack<SimplicialTreeNode_ptr> S;
		S.push(sigma);
		visited.insert(sigma);
		while (!S.empty()) {
			SimplicialTreeNode_ptr curr = S.top();
			bool to_be_visited = true;
			if (!curr->children.empty()) {
				for (std::unordered_map<int, SimplicialTreeNode_ptr>::iterator mIter = curr->children.begin();
					mIter != curr->children.end(); ++mIter) {
					if (visited.find(mIter->second) == visited.end()) {
						to_be_visited = false;
						S.push(mIter->second);
						visited.insert(mIter->second);
					}
				}
			}
			if (to_be_visited) {
				coBdries.push_back(curr);
				S.pop();
			}
		}
	}
	else {
		queue<SimplicialTreeNode_ptr> Q;
		Q.push(sigma);
		while (!Q.empty()) {
			SimplicialTreeNode_ptr curr = Q.front();
			Q.pop();
			coBdries.push_back(curr);
			for (std::unordered_map<int, SimplicialTreeNode_ptr>::iterator mIter = curr->children.begin();
				mIter != curr->children.end(); ++mIter) {
				Q.push(mIter->second);
			}
		}
	}
	return;
}
template<typename T>
bool SimplicialTree<T>::is_coface_set(SimplicialTreeNode_ptr &sigma, vector<int> &simplex_vertices, unordered_set<int>& hset) {
	// labels in simplex_vertices and sigma are sorted
	SimplicialTreeNode_ptr trav = sigma;
	int i = (int)simplex_vertices.size() - 1;
	while (trav) {
		if (i >= 0 && trav->label == simplex_vertices[i]) {
			trav = trav->parent;
			--i;
		}
		else {
			if ((i >= 0 && trav->label < simplex_vertices[i]) || (hset.find(trav->label) == hset.end())) {
				// the lagest remaining labels is smaller than the current traget
				//or there is a vertex not in common neighbor set
				// failed.
				return false;
			}
			trav = trav->parent;
		}
	}
	return (i < 0);
}
template<typename T>
bool SimplicialTree<T>::is_coface(SimplicialTreeNode_ptr &sigma, vector<int> &simplex_vertices) {
	// labels in simplex_vertices and sigma are sorted
	SimplicialTreeNode_ptr trav = sigma;
	int i = (int)simplex_vertices.size() - 1;
	while (trav && i >= 0) {
		if (trav->label == simplex_vertices[i]) {
			trav = trav->parent;
			--i;
		}
		else {
			if (trav->label < simplex_vertices[i]) {
				// the lagest remaining labels is smaller than the current traget
				// failed.
				return false;
			}
			trav = trav->parent;
		}
	}
	return (i < 0);
}

template<typename T>
bool SimplicialTree<T>::is_coface1(SimplicialTreeNode_ptr &sigma, vector<int> &simplex_vertices) {
	// labels in simplex_vertices and sigma are sorted
	int iMiss = 0;
	SimplicialTreeNode_ptr trav = sigma;
	int i = (int)simplex_vertices.size() - 1;
	while (trav && i >= 0) {
		if (trav->label == simplex_vertices[i]) {
			trav = trav->parent;
			--i;
		}
		else {
			if (iMiss == 0)
				iMiss += 1;
			else	//can't miss twice while still be a co-dimension-1 face
				return false;
			if (trav->label < simplex_vertices[i]) {
				// the lagest remaining labels is smaller than the current traget
				// failed.
				return false;
			}
			trav = trav->parent;
		}
	}
	return (i < 0);
}
template <typename T>
void SimplicialTree<T>::CoDimensionOneFaces(const SimplicialTreeNode_ptr& sigma, std::vector<SimplicialTreeNode_ptr> &codim_one_faces) {
	vector<int> simplex_vertices;
	retrieve_vertex_indices(sigma, simplex_vertices);
	int simplex_dim = (int)simplex_vertices.size() - 1;
	//
	for (std::unordered_map<int, SimplicialTreeNode_ptr>::iterator mIter = sigma->children.begin();
		mIter != sigma->children.end(); ++mIter) {
		codim_one_faces.push_back(mIter->second);
	}
	//
	int label = simplex_vertices.back();
	unordered_map<int, SimplicialTreeNode_ptr> & labels_dict = (*labels_dict_in_each_dim[simplex_dim + 1]);
	if (labels_dict.find(label) != labels_dict.end()) {
		// visit all simplices in current dimension containing label
		SimplicialTreeNode_ptr trav = labels_dict[label];
		SimplicialTreeNode_ptr start = trav;
		do {
			// check it contains all labels in the input simplex or not
			if (is_coface1(trav, simplex_vertices)) {
				codim_one_faces.push_back(trav);
			}
			trav = trav->next_circular_ptr;
		} while (trav != start);
	}
	return;
}

template <typename T>
void SimplicialTree<T>::CoBoundaryWithRestricSet(const SimplicialTreeNode_ptr &sigma, std::vector<SimplicialTreeNode_ptr> &coBdries, vector<int>& commonNB, int p, bool dfs_visit) {
	unordered_set<int> hset(commonNB.begin(), commonNB.end());
	// visit subtress using bfs
	vector<int> simplex_vertices;
	retrieve_vertex_indices(sigma, simplex_vertices);
	int simplex_dim = (int)simplex_vertices.size() - 1;
	//
	for (auto cd : sigma->children)
		if(hset.find(cd.first) != hset.end())
			retrieve_simplices_from_subtrees_set(cd.second, coBdries, hset, p);
	///////////
	int label = simplex_vertices.back();
	for (int i = simplex_dim + 1; i < labels_dict_in_each_dim.size() && i < p + 1; ++i) {
		unordered_map<int, SimplicialTreeNode_ptr> & labels_dict = (*labels_dict_in_each_dim[i]);
		if (labels_dict.find(label) != labels_dict.end()) {
			// visit all simplices in current dimension containing label
			SimplicialTreeNode_ptr trav = labels_dict[label];
			SimplicialTreeNode_ptr start = trav;
			do {

				// check it contains all labels in the input simplex or not
				if (is_coface_set(trav, simplex_vertices, hset)) {
					retrieve_simplices_from_subtrees_set(trav, coBdries, hset, p);
				}
				trav = trav->next_circular_ptr;
			} while (trav != start);
		}
	}
	return;
}

template <typename T>
bool SimplicialTree<T>::CoBoundary(const SimplicialTreeNode_ptr& sigma, std::vector<SimplicialTreeNode_ptr> &coBdries, bool dfs_visit){
	// visit subtress using bfs
	vector<int> simplex_vertices;
	retrieve_vertex_indices(sigma, simplex_vertices);
	int simplex_dim = (int)simplex_vertices.size() - 1;
	//
	for (std::unordered_map<int, SimplicialTreeNode_ptr>::iterator mIter = sigma->children.begin();
		mIter != sigma->children.end(); ++mIter) {
		retrieve_simplices_from_subtrees(mIter->second, coBdries, dfs_visit);
	}
	//
	int label = simplex_vertices.back();
	for (int i = simplex_dim + 1; i < labels_dict_in_each_dim.size(); ++i) {
		unordered_map<int, SimplicialTreeNode_ptr> & labels_dict = (*labels_dict_in_each_dim[i]);
		if (labels_dict.find(label) != labels_dict.end()) {
			// visit all simplices in current dimension containing label
			SimplicialTreeNode_ptr trav = labels_dict[label];
			SimplicialTreeNode_ptr start = trav;
			do {
				// check it contains all labels in the input simplex or not
				if (is_coface(trav, simplex_vertices)) {
					retrieve_simplices_from_subtrees(trav, coBdries, dfs_visit);
				}
				trav = trav->next_circular_ptr;
			} while (trav != start);
		}
	}
	return true;
}
template <typename T>
void SimplicialTree<T>::LinkSubcomplexSet(const SimplicialTreeNode_ptr &sigma, std::vector<SimplicialTreeNode_ptr> &simplexLink, vector<int>& commonNB, int p){
	//std::vector<SimplicialTreeNode_ptr> sl;
	CoBoundaryWithRestricSet(sigma, simplexLink, commonNB, p);
	//CoBoundary(sigma, sl);
	vector<int> simplex_vertices;
	retrieve_vertex_indices(sigma, simplex_vertices);
	vector<int> link_simplex_vertices;
	for (int i = 0; i < simplexLink.size(); ++i) {
		SimplicialTreeNode_ptr trav = simplexLink[i];
		int index = (int)simplex_vertices.size() - 1;
		while (trav) {
			if (index >= 0 && trav->label == simplex_vertices[index]) {
				trav = trav->parent;
				--index;
			}
			else {
				if (index < 0 || trav->label > simplex_vertices[index]) {
					link_simplex_vertices.push_back(trav->label);
					trav = trav->parent;
				}
				else {
					cout << "Not a coface in getting link simplex" << endl;
					exit(0);
				}
			}
		}
		//
		reverse(link_simplex_vertices.begin(), link_simplex_vertices.end());
		//simplicies.push_back(link_simplex_vertices);
		simplexLink[i] = find(link_simplex_vertices);
		link_simplex_vertices.clear();
	}
	return;
}
template <typename T>
void SimplicialTree<T>::LinkSubcomplex(const SimplicialTreeNode_ptr &sigma, std::vector<SimplicialTreeNode_ptr> &simplexLink) {
	CoBoundary(sigma, simplexLink);
	vector<int> simplex_vertices;
	retrieve_vertex_indices(sigma, simplex_vertices);
	vector<int> link_simplex_vertices;
	for (int i = 0; i < simplexLink.size(); ++i) {
		SimplicialTreeNode_ptr trav = simplexLink[i];
		int index = (int)simplex_vertices.size() - 1;
		while (trav) {
			if (index >= 0 && trav->label == simplex_vertices[index]) {
				trav = trav->parent;
				--index;
			}
			else {
				if (index < 0 || trav->label > simplex_vertices[index]) {
					link_simplex_vertices.push_back(trav->label);
					trav = trav->parent;
				}
				else {
					cout << "Not a coface in getting link simplex" << endl;
					exit(0);
				}
			}
		}
		//
		reverse(link_simplex_vertices.begin(), link_simplex_vertices.end());
		simplexLink[i] = find(link_simplex_vertices);
		link_simplex_vertices.clear();
	}
	return;
}
template <typename T>
void SimplicialTree<T>::FindCommonNeighbors(int a, int b, vector<int>& commonNB) {
	vector<int> aNB, bNB;
	//find all neighbors in <a, *>
	SimplicialTreeNode_ptr trav = (*labels_dict_in_each_dim[0])[a];
	for (auto p : trav->children)
		aNB.push_back(p.first);
	//find all neighbors in <*, a>
	if (labels_dict_in_each_dim[1]->find(a) != labels_dict_in_each_dim[1]->end()) {
		trav = (*labels_dict_in_each_dim[1])[a];
		SimplicialTreeNode_ptr start = trav;
		do {
			aNB.push_back(trav->parent->label);
			trav = trav->next_circular_ptr;
		} while (trav != start);
	}
	//find all neighbors in <b, *>
	trav = (*labels_dict_in_each_dim[0])[b];
	for (auto p : trav->children)
		bNB.push_back(p.first);
	//find all neighbors in <*, b>
	if (labels_dict_in_each_dim[1]->find(b) != labels_dict_in_each_dim[1]->end()) {
		trav = (*labels_dict_in_each_dim[1])[b];
		SimplicialTreeNode_ptr start = trav;
		do {
			bNB.push_back(trav->parent->label);
			trav = trav->next_circular_ptr;
		} while (trav != start);
	}
	if (aNB.size() < 20) {
		for (int m : bNB)
			for (int n : aNB)
				if (m == n)
					commonNB.push_back(m);
	}
	else if (bNB.size() < 20) {
		for (int m : aNB)
			for (int n : bNB)
				if (m == n)
					commonNB.push_back(m);
	}
	else if (aNB.size() < bNB.size()) {
		unordered_set<int> hset(aNB.begin(), aNB.end());
		for (int n : bNB)
			if (hset.find(n) != hset.end())
				commonNB.push_back(n);
	}
	else {
		unordered_set<int> hset(bNB.begin(), bNB.end());
		for (int n : aNB)
			if (hset.find(n) != hset.end())
				commonNB.push_back(n);
	}
}
template <typename T>
void SimplicialTree<T>::LinsubcomplexIntersection(SimplicialTreeNode_ptr & a, SimplicialTreeNode_ptr & b, boost::unordered_set<SimplicialTreeNode_ptr> &intersectedLinks, int p) {
	// no simplices in the subtree of sigma are in the intersection
	// 
	vector<int> commonNB;
	FindCommonNeighbors(a->label, b->label, commonNB);
	
	vector<SimplicialTreeNode_ptr> a_link_subcomplex;
	vector<SimplicialTreeNode_ptr> b_link_subcomplex;
	LinkSubcomplexSet(a, a_link_subcomplex, commonNB, p);
	LinkSubcomplexSet(b, b_link_subcomplex, commonNB, p);
	vector<SimplicialTreeNode_ptr> * smaller = &a_link_subcomplex;
	vector<SimplicialTreeNode_ptr> * larger = &b_link_subcomplex;
	if (a_link_subcomplex.size() > b_link_subcomplex.size()) {
		smaller = &b_link_subcomplex;
		larger = &a_link_subcomplex;
	}
	for (int i = 0; i < smaller->size(); ++i) {
		intersectedLinks.insert((*smaller)[i]);
	}
	int left = 0, right = (int)larger->size();
	while (left < right) {
		if (intersectedLinks.find((*larger)[left]) != intersectedLinks.end()) {
			++left;
		}
		else {
			swap((*larger)[left], (*larger)[right - 1]);
			--right;
		}
	}
	intersectedLinks.clear();
	for (int i = 0; i < left; ++i) {
		intersectedLinks.insert((*larger)[i]);
	}
	// 
	/*boost::unordered_set<SimplicialTreeNode_ptr> intersectedLinks2;
	vector<SimplicialTreeNode_ptr> a_link_subcomplex2;
	vector<SimplicialTreeNode_ptr> b_link_subcomplex2;
	LinkSubcomplex(a, a_link_subcomplex2);
	LinkSubcomplex(b, b_link_subcomplex2);
	vector<SimplicialTreeNode_ptr> * smaller2 = &a_link_subcomplex2;
	vector<SimplicialTreeNode_ptr> * larger2 = &b_link_subcomplex2;
	if (a_link_subcomplex2.size() > b_link_subcomplex2.size()) {
		smaller2 = &b_link_subcomplex2;
		larger2 = &a_link_subcomplex2;
	}
	for (int i = 0; i < smaller2->size(); ++i) {
		intersectedLinks2.insert((*smaller2)[i]);
	}
	left = 0, right = (int)larger2->size();
	while (left < right) {
		if (intersectedLinks2.find((*larger2)[left]) != intersectedLinks2.end()) {
			++left;
		}
		else {
			swap((*larger2)[left], (*larger2)[right - 1]);
			--right;
		}
	}
	intersectedLinks2.clear();
	for (int i = 0; i < left; ++i) {
		intersectedLinks2.insert((*larger2)[i]);
	}
	if (intersectedLinks.size() != intersectedLinks2.size())
		int a = 4;*/
	return;
}
//check for up to p-link condition, namely check any up to (p-1)-simplex in intersectedLinks subcomplx
template<typename T>
bool SimplicialTree<T>::is_upto_p_link_condition_satisfied(SimplicialTreeNode_ptr &edge, boost::unordered_set<SimplicialTreeNode_ptr> &intersectedLinks, int p) {
	unordered_map<int, SimplicialTreeNode_ptr> & vertices = (*labels_dict_in_each_dim[0]);
	SimplicialTreeNode_ptr a = vertices[edge->label]; // vertex with larger label
	SimplicialTreeNode_ptr b = vertices[edge->parent->label]; // vertex with smaller label
	LinsubcomplexIntersection(a, b, intersectedLinks, p);
	vector<SimplicialTreeNode_ptr> edgeLinkSubcomplex;
	LinkSubcomplex(edge, edgeLinkSubcomplex);
	//delete any simplex higher than (p-1)-simplex in intersectedLinks
	// by theory, edge link subcomplex is a subset of the intersection subcomplex
	/*boost::unordered_set<SimplicialTreeNode_ptr>::iterator iter = intersectedLinks.begin();
	boost::unordered_set<SimplicialTreeNode_ptr>::iterator itNext;
	for (; iter != intersectedLinks.end(); iter = itNext)
	{
		itNext = iter;
		++itNext;
		if (SimplexDim(*iter) >= p)
			intersectedLinks.erase(iter);
	}*/

	for (int i = 0; i < edgeLinkSubcomplex.size(); ++i) {
		if (SimplexDim(edgeLinkSubcomplex[i]) >= p)
			continue;
		intersectedLinks.erase(edgeLinkSubcomplex[i]);
	}
	return intersectedLinks.empty();
}

template<typename T>
bool SimplicialTree<T>::is_link_condition_satisfied(SimplicialTreeNode_ptr &edge, boost::unordered_set<SimplicialTreeNode_ptr> &intersectedLinks) {
	unordered_map<int, SimplicialTreeNode_ptr> & vertices = (*labels_dict_in_each_dim[0]);
	SimplicialTreeNode_ptr a = vertices[edge->label]; // vertex with larger label
	SimplicialTreeNode_ptr b = vertices[edge->parent->label]; // vertex with smaller label
	LinsubcomplexIntersection(a, b, intersectedLinks);
	vector<SimplicialTreeNode_ptr> edgeLinkSubcomplex;
	LinkSubcomplex(edge, edgeLinkSubcomplex);
	// by theory, edge link subcomplex is a subset of the intersection subcomplex
	for (int i = 0; i < edgeLinkSubcomplex.size(); ++i) {
		if (intersectedLinks.find(edgeLinkSubcomplex[i]) != intersectedLinks.end()) {
			intersectedLinks.erase(edgeLinkSubcomplex[i]);
		}
	}
	return intersectedLinks.empty();
}
template<typename T>
void SimplicialTree<T>::AddExtraSimplicesToSatisfyLinkCondition(SimplicialTreeNode_ptr &edge,
	boost::unordered_set<SimplicialTreeNode_ptr> &intersectedLinks){
	// since simplices in the intersection subcomplex are from the cofaces of the endpoints of the edge
	// only need to add simplicies [a b *]
	vector<int> edge_vertices(2, edge->parent->label);
	edge_vertices.back() = edge->label;
	//
	vector<int> link_simplex_vertices;
	vector<int> simplex_vertices;
	link_simplex_vertices.reserve(dim + 1);
	simplex_vertices.reserve(2 * dim + 2);

	vector<SimplicialTreeNode_ptr> baseSubcomplex(intersectedLinks.begin(), intersectedLinks.end());
	// no need to sort as simplices are ordered by dimension because of bfs visit
	if (baseSubcomplex.size() > 1) {
		sort(baseSubcomplex.begin(), baseSubcomplex.end(), SimplicialTreeNodePtrLessThan());
	}
	// use the fact that the intersection link subcomplexes of two vertices is a subcomplex as well
	for (int i = 0; i < baseSubcomplex.size(); ++i) {
		link_simplex_vertices.clear();
		simplex_vertices.clear();
		retrieve_vertex_indices(baseSubcomplex[i], link_simplex_vertices);
		merge_two_sorted_arrays(edge_vertices, link_simplex_vertices, simplex_vertices);
		//
		//if (simplex_vertices.size() <= HIGH_DIM + 1)
		ElementaryInsersion(simplex_vertices);
	}
	return;
}
template<typename T>
void SimplicialTree<T>::AnnotationTransfer(SimplicialTreeNode_ptr &edge, const int preserve_vertex_label) {
	// remove the cofaces of the edge which contain both vertices
	//cout << "AnnotationTransfer" << endl;
	const int removed_vertex_label = (edge->label == preserve_vertex_label ? edge->parent->label : edge->label);
	vector<SimplicialTreeNode_ptr> edge_cofaces;
	vector<SimplicialTreeNode_ptr> codim_one_faces;
	CoBoundary(edge, edge_cofaces, true);
	edge_cofaces.push_back(edge);
	// cofaces are ordered with dimension larger first 
	vector<int> simplex_vertices;
	simplex_vertices.reserve(dim + 1);
	for (int i = 0; i < edge_cofaces.size(); ++i) {
		// every simplex in the coface is a vanishing simplex
		// if the vanishing simpelx has non-zero annotation
		int simplex_dim = SimplexDim(edge_cofaces[i]);
		if (simplex_dim <= max_dimension)
		{
			ListNodePtr vanishing_simplex_annotation = find_annotation(edge_cofaces[i]);
			if (!is_zero_annotation(vanishing_simplex_annotation)) {
				//check_status();
				vanishing_simplex_annotation = annotations[simplex_dim]->DeepCopyAnnotationColumn(vanishing_simplex_annotation);
				simplex_vertices.clear();
				retrieve_vertex_indices(edge_cofaces[i], simplex_vertices);
				// find the mirror simplex
				for (int j = 0; j < simplex_vertices.size(); ++j) {
					if (simplex_vertices[j] > removed_vertex_label) {
						simplex_vertices[j - 1] = simplex_vertices[j];
					}
				}
				simplex_vertices.resize(simplex_vertices.size() - 1);
				//
				SimplicialTreeNode_ptr mirror_simplex = find(simplex_vertices);
				// obtain all codimension one faces of the mirror simplex
				codim_one_faces.clear();
				CoDimensionOneFaces(mirror_simplex, codim_one_faces);
				//
				// add the annotation of the vanishing simplex to all codimension one faces of the mirror simplex
				for (int j = 0; j < codim_one_faces.size(); ++j) {
					//
					//cout << j << "\t"; check_status(); 
					ListNodePtr codim_one_face_annotation = annotations[simplex_dim]->DeepCopyAnnotationColumn(find_annotation(codim_one_faces[j]));
					// delete the tree and the annotation it is pointed to
					if (ufd.Is_singleton(codim_one_faces[j]->tree_node)) {
						// after deletion also delete the annotation
						TreeRootNodePtr root = ufd.Find(codim_one_faces[j]->tree_node);
						//annotations[simplex_dim]->Delete(root->attribute);
						annotations[simplex_dim]->clearNode(root->attribute, false);
					}
					ufd.Delete(codim_one_faces[j]->tree_node);
					//sigma->tree_node.reset(); reset in ufd.Delete

					annotations[simplex_dim]->sum_two_annotation_with_changed_dst(codim_one_face_annotation, vanishing_simplex_annotation);
					//
					// insert the simplex with new annotation into the union-find-deletion forest and annotation matrix
					TreeRootNodePtr root = ufd.MakeSet(codim_one_faces[j]);
					annotations[simplex_dim]->Insert(codim_one_face_annotation, root, ufd);
					//
				}
			}
		}
	}
	//cout << " out AnnotationTransfer" << endl;
	return;
}
template<typename T>
void SimplicialTree<T>::ElementaryCollapse(int remove_label, int preserve_label) {
	if (remove_label == preserve_label)
		return;
	remove_label = findRoot(remove_label);
	preserve_label = findRoot(preserve_label);
	vector<int> simplex_vertices(2, remove_label);
	remove_label < preserve_label ? simplex_vertices[1] = preserve_label : simplex_vertices[0] = preserve_label;
	SimplicialTreeNode_ptr edge_simplex = find(simplex_vertices);
	cout<<"Collapse0";
	if (!edge_simplex) {
		// the edge is not exist
		cout<<"Collapse1";
		edge_simplex = ElementaryInsersion(simplex_vertices);
	}
	// 
	boost::unordered_set<SimplicialTreeNode_ptr> intersectedLinkSubcomplex;
	if (!is_upto_p_link_condition_satisfied(edge_simplex, intersectedLinkSubcomplex, max_dimension)) {
		// insert necessary simplices to ensure the link condition
		cout<<"Collapse2";
		AddExtraSimplicesToSatisfyLinkCondition(edge_simplex, intersectedLinkSubcomplex);
	}
	//always collpase larger index to smaller index and update reindex
	if (remove_label < preserve_label) {
		cout<<"Collapse3";
		reindex[preserve_label] = remove_label;
		std::swap(remove_label, preserve_label);
	}
	// apply annotation transfer
	AnnotationTransfer(edge_simplex, preserve_label);
	// perform the edge contraction
	EdgeContraction(edge_simplex, preserve_label);
	//
	return;
}
template<typename T>
void SimplicialTree<T>::EdgeContraction(SimplicialTreeNode_ptr &edge, const int preserve_vertex_label) {
	unordered_map<int, SimplicialTreeNode_ptr> & vertices = (*labels_dict_in_each_dim.front());
	SimplicialTreeNode_ptr preserved_vertex = vertices[preserve_vertex_label];
	SimplicialTreeNode_ptr removed_vertex = vertices[edge->label == preserve_vertex_label ? edge->parent->label : edge->label];
	// remove the cofaces of the edge which contain both vertices
	vector<SimplicialTreeNode_ptr> cofaces;
	CoBoundary(edge, cofaces, true);
	// cofaces are ordered with dimension larger first
	cofaces.push_back(edge); // delete this edge as well
	for (int i = 0; i < cofaces.size(); ++i) {
		remove_simplex_from_both_complex_and_ufdForest(cofaces[i]);
		//PrintComplexWithAnnotation();
	}
	// handle the vertex first

	rename_simplices_in_subtree(removed_vertex, preserve_vertex_label);
	//PrintComplexWithAnnotation();
	// 
	int simplex_dim = 0;
	int label = removed_vertex->label;
	for (int i = simplex_dim + 1; i < labels_dict_in_each_dim.size(); ++i) {
		unordered_map<int, SimplicialTreeNode_ptr> & labels_dict = (*labels_dict_in_each_dim[i]);
		while (labels_dict.find(label) != labels_dict.end()) {
			// visit all simplices in current dimension containing the removed label 
			SimplicialTreeNode_ptr trav = labels_dict[label];
			rename_simplices_in_subtree(trav, preserve_vertex_label);
			//PrintComplexWithAnnotation();
		}
	}
	// update the union-find-deletion forest and annotation matrix
	return;
}
template<typename T>
void SimplicialTree<T>::rename_simplices_in_subtree(SimplicialTreeNode_ptr sigma, const int new_label){
	// record the deletion option for the simplices
	// 1) remove only from the simplicial tree -- false
	// 2) remove both from the simplicial tree and union-find-deletion forest --- true
	boost::unordered_map<SimplicialTreeNode_ptr, pair<SimplicialTreeNode_ptr, bool> > flags;
	stack<SimplicialTreeNode_ptr> S;
	vector<int> simplex_vertices;
	retrieve_vertex_indices(sigma, simplex_vertices);
	// sigma->label is replaced with new_label
	simplex_vertices.back() = new_label;
	sort(simplex_vertices.begin(), simplex_vertices.end());
	SimplicialTreeNode_ptr newNode = find(simplex_vertices);
	if (!newNode) {
		// this is a new simplex to be inserted
		//insert the the new one into the simplicial tree
		newNode = insert_into_simplicial_tree(simplex_vertices);
		if (SimplexDim(sigma) <= max_dimension)
		{
			newNode->tree_node = sigma->tree_node;
			sigma->tree_node->elem = newNode;
			sigma->tree_node.reset();
		}
		// 
		flags[sigma] = pair<SimplicialTreeNode_ptr, bool>(newNode, false); // only delete from simplicial tree
	}
	else {
		flags[sigma] = pair<SimplicialTreeNode_ptr, bool>(newNode, true); // delete from both simplicial tree and union-find-deletion forest
	}

	S.push(sigma);
	while (!S.empty()) {
		SimplicialTreeNode_ptr curr = S.top();
		simplex_vertices.clear();
		retrieve_vertex_indices(flags[curr].first, simplex_vertices);
		//
		if (!curr->children.empty() && flags.find(curr->children.begin()->second) == flags.end()) {
			// all childrens are not visited 
			// check the subtrees
			simplex_vertices.push_back(0); // one dimension higher simplex
			vector<int> ordered_children;
			ordered_children.reserve(curr->children.size());
			for (std::unordered_map<int, SimplicialTreeNode_ptr>::iterator mIter = curr->children.begin();
				mIter != curr->children.end(); ++mIter) {
				ordered_children.push_back(mIter->first);
			}
			// 
			int read = 0;
			vector<int> backup_simplex = simplex_vertices;
			///////////
			while (read < ordered_children.size()) {
				SimplicialTreeNode_ptr child = curr->children[ordered_children[read]];
				// update the simplex vertices array
				simplex_vertices = backup_simplex;
				simplex_vertices.back() = child->label;
				if (simplex_vertices.back() < simplex_vertices[simplex_vertices.size() - 2]) {
					swap(simplex_vertices.back(), simplex_vertices[simplex_vertices.size() - 2]);
				}
				newNode = find(simplex_vertices);
				//
				if (!newNode) {
					if (ordered_children[read] > new_label) {
						// no need to generate a new node;
						// move this node and all children to the new position
						SimplicialTreeNode_ptr parent = flags[curr].first;
						parent->children[child->label] = child;
						// remove it from it parent
						child->parent->children.erase(child->label);
						// update its new parent
						child->parent = parent;
						//do not move write pointer
					}
					else {
						// copy the node in union-find-deletion forest and 
						newNode = insert_into_simplicial_tree(simplex_vertices);
						if (SimplexDim(newNode) <= max_dimension)
						{
							newNode->tree_node = child->tree_node;
							child->tree_node->elem = newNode;
							child->tree_node.reset();
						}
						//
						S.push(child);
						flags[child] = pair<SimplicialTreeNode_ptr, bool>(newNode, false); // only delete from simplicial tree 
					}
				}
				else {
					S.push(child);
					flags[child] = pair<SimplicialTreeNode_ptr, bool>(newNode, true); // delete from both simplicial tree and union-find-deletion forest
				}
				++read;
			}
		}
		else {
			// no child or all children are visited before  
			// second visit
			// delete from the simplicial tree
			S.pop();
			// perform the deletion action
			if (flags[curr].second) {
				// remove it from both simplicial tree and union-find-deletion forest
				remove_simplex_from_both_complex_and_ufdForest(curr);
			}
			else {
				remove_simplex_from_complex(curr);
			}
		}
	}
	return;
}
template<typename T>
void SimplicialTree<T>::remove_simplex_from_complex(SimplicialTreeNode_ptr sigma){
	// Preq: all cofaces are deleted already
	// remove the tree node in union-find-deletion data structure
	if (!sigma->children.empty()) {
		cout << "can not be deleted as some cofaces exists" << endl;
		exit(0);
	}
	int simplex_dim = SimplexDim(sigma);
	//sigma->tree_node.reset(); reset in ufd.Delete
	// parent-children relation
	sigma->parent->children.erase(sigma->label);
	sigma->parent.reset();
	//delete from the circular list or the unordered map
	delete_from_circular_list(sigma, simplex_dim);
	// update size info
	--simplex_sizes[simplex_dim];
	if (simplex_sizes[simplex_dim] == 0) {
		simplex_sizes.resize(simplex_dim);
		--dim;
		labels_dict_in_each_dim.resize(simplex_dim);
		//keep the time stamp for annotation matrix of dim: simplex_dim
		this->vecTS[simplex_dim] = annotations[simplex_dim]->timeStamp;
		annotations.resize(simplex_dim);
	}
	return;
}

template<typename T>
void SimplicialTree<T>::remove_simplex_from_both_complex_and_ufdForest(SimplicialTreeNode_ptr sigma, bool bUpdatePers){
	// Preq: all cofaces are deleted already
	// remove the tree node in union-find-deletion data structure
	if (!sigma->children.empty()) {
		cout << "can not be deleted as some cofaces exists" << endl;
		exit(0);
	}
	int simplex_dim = SimplexDim(sigma);
	if (simplex_dim <= max_dimension)
	{
		if (ufd.Is_singleton(sigma->tree_node)) {
			// after deletion also delete the annotation
			TreeRootNodePtr root = ufd.Find(sigma->tree_node);
			annotations[simplex_dim]->clearNode(root->attribute, bUpdatePers);
		}
		ufd.Delete(sigma->tree_node);
	}
	//sigma->tree_node.reset(); reset in ufd.Delete
	// parent-children relation
	if (simplex_dim > 0) {
		// the parent of a vertex is null
		sigma->parent->children.erase(sigma->label);
		sigma->parent.reset();
	}
	//delete from the circular list or the unordered map
	delete_from_circular_list(sigma, simplex_dim);
	//break cycle
	//sigma->tree_node.reset();
	// update size info
	--simplex_sizes[simplex_dim];
	if (simplex_sizes[simplex_dim] == 0) {
		simplex_sizes.resize(simplex_dim);
		--dim;
		//labels_dict_in_each_dim[simplex_dim].reset();
		labels_dict_in_each_dim.resize(simplex_dim);
		//keep the time stamp for annotation matrix of dim: simplex_dim
		if (simplex_dim <= max_dimension)
			this->vecTS[simplex_dim] = annotations[simplex_dim]->timeStamp;
		annotations.resize(simplex_dim);
	}
	return;
}
/*----------------Persistence Related operations-----------------------------------*/
template <typename T>
void SimplicialTree<T>::SnapshotHomologicalFeatures(vector<unordered_set<int> > &hom_info) {
	for (int i = 0; i < annotations.size(); ++i) {
		unordered_set<int> cycles;
		for (unordered_map<int, ListNodePtr>::iterator mIter = annotations[i]->row_ptr.begin();
			mIter != annotations[i]->row_ptr.end(); ++mIter) {
			cycles.insert(mIter->first);
		}
		hom_info.push_back(cycles);
	}
	return;
}
template <typename T>
void SimplicialTree<T>::CheckPersistence(vector<unordered_set<int> > &homo_info) {
	for (int i = (int)annotations.size(); i < homo_info.size(); ++i) {
		// all homology cycles with dimension higher than the dimension of the complex are killed
		homo_info[i].clear();
	}
	for (int i = 0; i < annotations.size(); ++i) {
		if (i < homo_info.size()) {
			vector<int> erase_set;
			erase_set.reserve(homo_info[i].size() + 1);
			for (unordered_set<int>::iterator sIter = homo_info[i].begin(); sIter != homo_info[i].end(); ++sIter) {
				if (annotations[i]->row_ptr.find(*sIter) == annotations[i]->row_ptr.end()) {
					erase_set.push_back(*sIter);
				}
			}
			for (int j = 0; j < erase_set.size(); ++j) {
				homo_info[i].erase(erase_set[j]);
			}
		}
		else {
			break;
		}
	}
	return;
}
template <typename T>
void SimplicialTree<T>::PerformSimplicialCollapse(vector<pair<int, int> > &vertex_map, unordered_map<int,int> &updated_vertex_map) {
	// record the non-singleton pre-images of the vertex_map
	//vector<int> empVec;
	vector<vector<int>> preimages(maxImageVertex + 1);
	for (int i = 0; i < vertex_map.size(); ++i) {
		preimages[vertex_map[i].second].push_back(vertex_map[i].first);
	}
	for (int i = 0; i < preimages.size(); ++i) {
		//cout << mIter->first << endl;
		if (!preimages[i].empty())
		{
			if (preimages[i].size() > 1) {
				//sort(mIter->second.begin(), mIter->second.end());
				//for (int i = 0; i < mIter->second.size() - 1; ++i) {
				//	ElementaryCollapse((mIter->second)[i], mIter->second.back());
				//} 
				for (int j = preimages[i].size() - 1; j > 0; --j) {
					//check_status();
					ElementaryCollapse(preimages[i][j], preimages[i].front());
					//cout << "after " << endl;
					//check_status();
				}
			}
			updated_vertex_map[preimages[i].front()] = i;
		}
	}
	return;
}
template <typename T>
void SimplicialTree<T>::clearMemory()		//clear memory of annotation matrices, simplicial trees and udf trees
{
	//clear simplicial trees, ufd trees and annotation matrices
	for (int i = labels_dict_in_each_dim.size() - 1; i >= 0; --i) {

		unordered_map<int, SimplicialTreeNode_ptr> & labels_dict = *(labels_dict_in_each_dim[i]);
		//
		unordered_map<int, SimplicialTreeNode_ptr>::iterator mIter = labels_dict.begin();
		if (i == 0)
		{
			// handle vertices
			while (true)
			{
				if (labels_dict.size() == 1)
				{
					this->remove_simplex_from_both_complex_and_ufdForest(labels_dict.begin()->second, false);
					break;
				}
				this->remove_simplex_from_both_complex_and_ufdForest(labels_dict.begin()->second, false);
			}
		}
		else {
			while (true)
			{
				SimplicialTreeNode_ptr trav = labels_dict.begin()->second;
				SimplicialTreeNode_ptr next;
				if (labels_dict.size() == 1)
				{
					while (trav != trav->next_circular_ptr)
					{
						next = trav->next_circular_ptr;
						this->remove_simplex_from_both_complex_and_ufdForest(trav, false);
						trav = next;
					}
					this->remove_simplex_from_both_complex_and_ufdForest(trav, false);
					break;
				}
				while (trav != trav->next_circular_ptr)
				{
					next = trav->next_circular_ptr;
					this->remove_simplex_from_both_complex_and_ufdForest(trav, false);
					trav = next;
				}
				this->remove_simplex_from_both_complex_and_ufdForest(trav, false);
			}
		}
	}
}

//relabeling vertices
template <typename T>
void SimplicialTree<T>::RelabelingVertices(unordered_map<int, int> &vertex_map)
{
	std::vector<LabelsDictionaryPtr> new_labels_dict_in_each_dim;
	for (int i = 0; i < this->labels_dict_in_each_dim.size(); ++i)
	{
		unordered_map<int, SimplicialTreeNode_ptr> & labels_dict = *(labels_dict_in_each_dim[i]);
		//add a new dimension in new_labels_dict_in_each_dim
		new_labels_dict_in_each_dim.push_back(boost::make_shared<unordered_map<int, SimplicialTreeNode_ptr> >());
		//relabeling each simplical tree node
		//also relabel its children map
		unordered_map<int, SimplicialTreeNode_ptr>::iterator mIter = labels_dict.begin();
		for (; mIter != labels_dict.end(); ++mIter)
		{
			if (i == 0)
			{
				//handle vertices
				//update current tree node label 
				mIter->second->label = vertex_map[mIter->first];
				//update children map
				if (mIter->second->children.size() != 0)
				{
					std::unordered_map<int, SimplicialTreeNode_ptr> newChildren;
					std::unordered_map<int, SimplicialTreeNode_ptr>::iterator mIterChildren = mIter->second->children.begin();
					for (; mIterChildren != mIter->second->children.end(); ++mIterChildren)
					{
						newChildren[vertex_map[mIterChildren->first]] = mIterChildren->second;
					}
					mIter->second->children = newChildren;
				}
			}
			else
			{
				SimplicialTreeNode_ptr trav = mIter->second;
				do
				{
					//update current tree node label
					trav->label = vertex_map[mIter->first];
					if (trav->children.size() != 0)
					{
						//update children map
						std::unordered_map<int, SimplicialTreeNode_ptr> newChildren;
						std::unordered_map<int, SimplicialTreeNode_ptr>::iterator mIterChildren = trav->children.begin();
						for (; mIterChildren != trav->children.end(); ++mIterChildren)
						{
							newChildren[vertex_map[mIterChildren->first]] = mIterChildren->second;
						}
						trav->children = newChildren;
					}
					//next tree node
					trav = trav->next_circular_ptr;
				} while (trav != mIter->second);
			}
			//copy new labels to new_labels_dict in each dimension
			//add a new dimension in new_labels_dict_in_each_dim
			unordered_map<int, SimplicialTreeNode_ptr> & new_labels_dict = *(new_labels_dict_in_each_dim[i]);
			new_labels_dict[vertex_map[mIter->first]] = labels_dict[mIter->first];
		}
	}
	//relabeling labels_dict_in_each_dim by copying new_labels_dict_in_each_dim back
	//labels_dict_in_each_dim.clear();
	labels_dict_in_each_dim = new_labels_dict_in_each_dim;
	return;
}


template <typename T>
void SimplicialTree<T>::InitializeByRenamingIncomingComplex(SimplicialTree<T> &src, unordered_map<int, int> &vertex_map){
	//clear before use
	clearData();

	vector<int> simplex_vertices;
	simplex_vertices.reserve(dim + 1);
	// visit each simplex through the labels_diction variables
	for (int i = 0; i < src.labels_dict_in_each_dim.size(); ++i) {

		unordered_map<int, SimplicialTreeNode_ptr> & labels_dict = *(src.labels_dict_in_each_dim[i]);
		//
		unordered_map<int, SimplicialTreeNode_ptr>::iterator mIter = labels_dict.begin();
		if (i == 0) {
			// handle vertices
			simplex_vertices.resize(1);
			for (; mIter != labels_dict.end(); ++mIter) {
				ListNodePtr simplex_ann = src.annotations[i]->DeepCopyAnnotationColumn(src.find_annotation(mIter->second));
				simplex_vertices.front() = vertex_map[mIter->first];
				//
				InsertSimplexWithAnnotation(simplex_vertices, simplex_ann);
			}
		}
		else {
			for (; mIter != labels_dict.end(); ++mIter) {
				SimplicialTreeNode_ptr trav = mIter->second;
				do{
					simplex_vertices.clear();
					src.retrieve_vertex_indices(trav, simplex_vertices);
					for (int j = 0; j < simplex_vertices.size(); ++j) {
						simplex_vertices[j] = vertex_map[simplex_vertices[j]];
					}
					sort(simplex_vertices.begin(), simplex_vertices.end());
					// 
					//  
					if (i <= max_dimension)
					{
						ListNodePtr simplex_ann = src.annotations[i]->DeepCopyAnnotationColumn(src.find_annotation(trav));
						InsertSimplexWithAnnotation(simplex_vertices, simplex_ann);
					}
					else
						ElementaryInsersion(simplex_vertices);
					//
					trav = trav->next_circular_ptr;
				} while (trav != mIter->second);
			}
		}
	}
	//copy annotation time stamps
	for (int i = 0; i < src.annotations.size(); ++i)
	{
		this->annotations[i]->timeStamp = src.annotations[i]->timeStamp;
	}
	this->vecTS = src.vecTS;
	return;
}
template <typename T>
void SimplicialTree<T>::AddRemainingSimpliciesFromFile(const char* pFileName) {
	/*
	# each line is a simplex with sorted integer vertex labels
	# ver_i < ver_j if i < j
	# ver_0 ver_1 ... ver_dim
	The following is a triangle (triangle.txt)
	0
	1
	2
	0 1
	1 2
	2 3
	1 2 3
	*/
	ReadComplex(pFileName);
	return;
}
/*----------------simplicial complex I/O-----------------------------------*/
template <typename T>
void SimplicialTree<T>::ReadSimplicialMap(const char* pFileName, vector<pair<int, int> > &vertex_map) {
	/*
	# Each line is a pair of vertices
	# Format: v_i v_j
	# Remark: v_i in the domain complex is mapped to v_j in the range complex
	The following is the map for edge contraction <1, 2> to 2
	1 2
	2 2
	*/
	this->maxImageVertex = -1;
	std::ifstream ifile;
	ifile.open(pFileName, std::ifstream::in | std::ifstream::binary);
	//
	//std::string sBuf;
	//std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	//
	if (ifile.is_open()) {
		//      /*
		//       * Get the size of the file
		//       */ 
		//ifile.seekg(0, ifile.end);
		//      long long iFileSize = ifile.tellg();
		//ifile.seekg(0, ifile.beg);

		////copy the whole file into file buffer
		//char* fBuf = new char[iFileSize + 1];
		//ifile.read(fBuf, iFileSize);
		//// add extra symbol
		//fBuf[iFileSize] = '\n';
		////
		//sBuf.assign(fBuf);
		//cout << sBuf << endl;
		//sstr.str(sBuf);

		//// close file
		////ifile.close();
		//ifile.seekg(0,std::ios::beg);
		//sBuf.clear();
		////deallocate memory
		//delete [] fBuf;
		///*start reading*/
		///*read simplex from each line*/
		string line;
		istringstream iss;
		//first line for current time stamp
		getline(ifile, line);
		iss.str(line);
		string delim;
		float dFilScale;
		iss >> delim >> dFilScale;
		vecFiltrationScale.push_back(dFilScale);
		//following lines of vertex maps
		while (getline(ifile, line)) {
			iss.str(line);
			int src = 0, dst = 0;
			if (iss >> src >> dst) {
				vertex_map.push_back(pair<int, int>(src, dst));
				//update max image vertex
				if (dst > maxImageVertex)
					maxImageVertex = dst;
			}
			line.clear();
			iss.clear();
			iss.str("");
		}
		ifile.close();
	}
	else {
		std::cout << "Can NOT open file: " << pFileName << std::endl;
		exit(0);
	}
}
template <typename T>
void SimplicialTree<T>::ReadComplexWithAnnotation(const char* pFileName) {
	/*
	# The first line contains an integer which is the total number of simplices in this simplicial complex
	# Each later line is a simplex with increasingly sorted integer vertex labels and its annotation []
	# format: ver_0 ver_1 ... ver_dim [<a_1 val_1> <a_2 val_2> ... <a_n val_n>]
	#			ver_i < ver_j if i < j,
	#			<a_i val_i> : means the a_i-th bit has non-zero value val_i,
	#			[] means zero annotation
	# Remark: annotation is assumed to be under Z_2
	The following is an empty triangle (triangle.txt)
	7           --> number of simplices
	0 [<0 1>]
	1 [<0 1>]
	2 [<0 1>]
	0 1 []
	1 2 []
	0 2 [<0 1>]
	*/
	std::ifstream ifile;
	ifile.open(pFileName, std::ifstream::in | std::ifstream::binary);
	//
	std::string sBuf;
	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	//
	if (ifile.is_open()) {
		/*
		* Get the size of the file
		*/
		ifile.seekg(0, std::ios::end);
		long long iFileSize = ifile.tellg();
		ifile.seekg(0, std::ios::beg);

		//copy the whole file into file buffer
		char* fBuf = new char[iFileSize + 1];
		ifile.read(fBuf, iFileSize);
		// add extra symbol
		fBuf[iFileSize] = '\n';
		//
		sBuf.assign(fBuf);
		sstr.str(sBuf);

		// close file
		ifile.close();
		sBuf.clear();
		//deallocate memory
		delete[] fBuf;
		/*start reading*/
		/*read simplex from each line*/
		string line;
		istringstream iss;
		vector<int> simplex_vertices;
		simplex_vertices.reserve(100);
		if (getline(sstr, line)) {
			iss.str(line);
			int simplex_count = 0;
			iss >> simplex_count;
			iss.clear();
			iss.str("");
			while (getline(sstr, line)) {
				simplex_vertices.clear();
				iss.str(line);
				int vertex = 0;
				while (iss >> vertex) {
					simplex_vertices.push_back(vertex);
				}
				if (simplex_count > 0 && !simplex_vertices.empty()) {
					--simplex_count;
					ListNodePtr simplex_ann = string_to_annotation(line);
					InsertSimplexWithAnnotation(simplex_vertices, simplex_ann);
					//initialize time stamp for each annotation matrix
					int simplexDimension = simplex_vertices.size() - 1;
					long long lastBit = annotations[simplexDimension]->lowest_one(simplex_ann);
					if (lastBit >= annotations[simplexDimension]->timeStamp)
						annotations[simplexDimension]->timeStamp = lastBit + 1;
				}
				line.clear();
				iss.clear();
				iss.str("");
			}
		}
	}
	else {
		std::cout << "Can NOT open file: " << pFileName << std::endl;
		exit(0);
	}
	return;
}
template <typename T>
void SimplicialTree<T>::ReadComplex(const char* pFileName) {
	/*
	# The first line contains an integer which is the total number of simplices in this simplicial complex
	# Each later line is a simplex with increasingly sorted integer vertex labels
	# Format: ver_0 ver_1 ... ver_dim
	# Remark:	1) ver_i < ver_j if i < j
	#			2) the simplices in the file are ordered such that the boundary faces of any simplex appear before itself
	#			3) the simple order satisfying this requirement is to sort the simplices by the increasingly order of their dimensions

	The following is a triangle (triangle.txt)
	7		---> number of simplices
	0		---> vertex 0
	1		---> vertex 1
	2		---> vertex 2
	0 1		---> edge <0, 1>
	1 2		---> edge <1, 2>
	0 2		---> edge <0, 2>
	0 1 2	---> triangle <0, 1, 2>
	*/
	std::ifstream ifile;
	ifile.open(pFileName, std::ifstream::in | std::ifstream::binary);
	//
	std::string sBuf;
	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	//
	if (ifile.is_open()) {
		/*
		* Get the size of the file
		*/
		ifile.seekg(0, std::ios::end);
		long long iFileSize = ifile.tellg();
		ifile.seekg(0, std::ios::beg);

		//copy the whole file into file buffer
		char* fBuf = new char[iFileSize + 1];
		ifile.read(fBuf, iFileSize);
		// add extra symbol
		fBuf[iFileSize] = '\n';
		//
		sBuf.assign(fBuf);
		sstr.str(sBuf);

		// close file
		ifile.close();
		sBuf.clear();
		//deallocate memory
		delete[] fBuf;
		/*start reading*/
		/*read simplex from each line*/
		string line;
		istringstream iss;
		vector<int> simplex_vertices;
		simplex_vertices.reserve(100);
		if (getline(sstr, line)) {
			iss.str(line);
			int simplex_count = 0;
			iss >> simplex_count;
			iss.clear();
			iss.str("");
			while (getline(sstr, line)) {
				simplex_vertices.clear();
				iss.str(line);
				int vertex = 0;
				while (iss >> vertex) {
					simplex_vertices.push_back(vertex);
				}
				if (simplex_count > 0 && !simplex_vertices.empty()) {
					--simplex_count;
					std::clock_t timer = std::clock();
					ElementaryInsersion(simplex_vertices);
					dInsertTime += std::clock() - timer;
				}
				line.clear();
				iss.clear();
				iss.str("");
			}
		}
	}
	else {
		std::cout << "Can NOT open file: " << pFileName << std::endl;
		exit(0);
	}
}

template <class T>
void SimplicialTree<T>::WriteStatisticsToFile(const char* pFileName)
{
	std::cout << "Writing statistics < " << pFileName << " >" << std::endl;
	std::ofstream ofile;
	ofile.open(pFileName, std::ifstream::out);
	//
	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	if (ofile.is_open())
	{
		//sstr << pFileName << std::endl;

		for (int i = 0; i < dim + 1; i++)
		{
			sstr << "DIM " << i << " : " << simplex_sizes[i] << std::endl;
		}
		sstr << "Total size : " << ComplexSize() << std::endl;
		//
		//ofile << sstr.rdbuf();
		ofile.write(sstr.str().c_str(), sstr.str().size());
		//
		ofile.close();
		sstr.clear();
	}
	else
	{
		std::cout << "Can NOT open file " << pFileName << std::endl;
		exit(0);
	}
	std::cout << "---- Done ----" << std::endl << std::endl;
	return;
}
template<typename T>
void SimplicialTree<T>::PrintComplexWithAnnotation() {
	vector<int> simplex_vertices;
	simplex_vertices.reserve(dim + 1);
	// visit each simplex through the labels_diction variables
	for (int i = 0; i < labels_dict_in_each_dim.size(); ++i) {

		unordered_map<int, SimplicialTreeNode_ptr> & labels_dict = *(labels_dict_in_each_dim[i]);
		// assigning the bits with new indices
		unordered_map<int, int> newIndices;
		vector<int> nonzero_bits;
		nonzero_bits.reserve(annotations[i]->row_ptr.size());
		for (unordered_map<int, ListNodePtr>::iterator bitIter = annotations[i]->row_ptr.begin();
			bitIter != annotations[i]->row_ptr.end(); ++bitIter) {
			nonzero_bits.push_back(bitIter->first);
		}
		sort(nonzero_bits.begin(), nonzero_bits.end());
		for (int j = 0; j < nonzero_bits.size(); ++j) {
			newIndices[nonzero_bits[j]] = j;
		}
		//
		unordered_map<int, SimplicialTreeNode_ptr>::iterator mIter = labels_dict.begin();
		if (i == 0) {
			// handle vertices
			for (; mIter != labels_dict.end(); ++mIter) {
				cout << mIter->first << " " << annotation_to_string(mIter->second, newIndices) << endl;
			}
		}
		else {
			for (; mIter != labels_dict.end(); ++mIter) {
				SimplicialTreeNode_ptr trav = mIter->second;
				do{
					simplex_vertices.clear();
					retrieve_vertex_indices(trav, simplex_vertices);
					for (int j = 0; j < simplex_vertices.size(); ++j) {
						cout << simplex_vertices[j] << (j == simplex_vertices.size() - 1 ? "" : " ");
					}
					cout << " " << annotation_to_string(trav, newIndices) << endl;
					//
					trav = trav->next_circular_ptr;
				} while (trav != mIter->second);
			}
		}
	}
}
template<typename T>
void SimplicialTree<T>::WriteComplexWithAnnotation(const char* pFileName) {
	/*
	# The first line contains an integer which is the total number of simplices in this simplicial complex
	# Each later line is a simplex with increasingly sorted integer vertex labels and its annotation []
	# format:
	# ver_0 ver_1 ... ver_dim [<a_1 val_1> <a_2 val_2> ... <a_n val_n>]
	# Remark 1:	1) ver_i < ver_j if i < j,
	#			2) <a_i val_i> : means the a_i-th bit has non-zero value val_i,
	#			3) [] means zero annotation
	#			4) the simplices in the file are ordered such that the boundary faces of any simplex appear before itself
	#			5) the simple order satisfying this requirement is to sort the simplices by the increasingly order of their dimensions

	# Remark 2: annotation is assumed to be under Z_2
	The following is an empty triangle (triangle.txt)
	6		--> number of simplices
	0 []
	1 []
	2 []
	0 1 []
	1 2 []
	0 2 [<0 1>]
	*/
	std::ofstream ofile;
	ofile.open(pFileName, std::ifstream::out);
	//
	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	if (ofile.is_open())
	{
		vector<int> simplex_vertices;
		simplex_vertices.reserve(dim + 1);
		// visit each simplex through the labels_diction variables
		sstr << ComplexSize() << endl;
		for (int i = 0; i < labels_dict_in_each_dim.size(); ++i) {

			unordered_map<int, SimplicialTreeNode_ptr> & labels_dict = *(labels_dict_in_each_dim[i]);
			// assigning the bits with new indices
			//unordered_map<int, int> newIndices;
			//vector<int> nonzero_bits;
			//nonzero_bits.reserve(annotations[i]->row_ptr.size());
			//for (unordered_map<int, ListNodePtr>::iterator bitIter = annotations[i]->row_ptr.begin();
			//	bitIter != annotations[i]->row_ptr.end(); ++bitIter) {
			//	nonzero_bits.push_back(bitIter->first);
			//}
			//sort(nonzero_bits.begin(), nonzero_bits.end());
			//for (int j = 0; j < nonzero_bits.size(); ++j) {
			//	newIndices[nonzero_bits[j]] = j;
			//}
			//
			unordered_map<int, SimplicialTreeNode_ptr>::iterator mIter = labels_dict.begin();
			if (i == 0) {
				// handle vertices
				for (; mIter != labels_dict.end(); ++mIter) {
					sstr << mIter->first << " " << annotation_to_string(mIter->second) << endl;
				}
			}
			else {
				for (; mIter != labels_dict.end(); ++mIter) {
					SimplicialTreeNode_ptr trav = mIter->second;
					do{
						simplex_vertices.clear();
						retrieve_vertex_indices(trav, simplex_vertices);
						for (int j = 0; j < simplex_vertices.size(); ++j) {
							sstr << simplex_vertices[j] << (j == simplex_vertices.size() - 1 ? "" : " ");
						}
						sstr << " " << annotation_to_string(trav) << endl;
						//
						trav = trav->next_circular_ptr;
					} while (trav != mIter->second);
				}
			}
		}
		//
		//ofile << sstr.rdbuf();
		ofile.write(sstr.str().c_str(), sstr.str().size());
		//
		ofile.close();
		sstr.clear();
	}
	else
	{
		std::cout << "Can NOT open file " << pFileName << std::endl;
		exit(0);
	}
}
template<typename T>
void SimplicialTree<T>::WriteComplex(const char* pFileName) {
	/*
	# The first line contains an integer which is the total number of simplices in this simplicial complex
	# Each later line is a simplex with increasingly sorted integer vertex labels
	# Format: ver_0 ver_1 ... ver_dim
	# Remark: ver_i < ver_j if i < j
	The following is a triangle (triangle.txt)
	7			-->number of simplices
	0
	1
	2
	0 1
	1 2
	0 2
	0 1 2
	*/
	std::ofstream ofile;
	ofile.open(pFileName, std::ifstream::out);
	//
	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	if (ofile.is_open())
	{
		vector<int> simplex_vertices;
		simplex_vertices.reserve(dim + 1);
		// visit each simplex through the labels_diction variables
		sstr << ComplexSize() << endl;
		for (int i = 0; i < labels_dict_in_each_dim.size(); ++i) {

			unordered_map<int, SimplicialTreeNode_ptr> & labels_dict = *(labels_dict_in_each_dim[i]);
			unordered_map<int, SimplicialTreeNode_ptr>::iterator mIter = labels_dict.begin();
			if (i == 0) {
				// handle vertices
				for (; mIter != labels_dict.end(); ++mIter) {
					sstr << mIter->first << endl;
				}
			}
			else {
				for (; mIter != labels_dict.end(); ++mIter) {
					SimplicialTreeNode_ptr trav = mIter->second;
					do{
						simplex_vertices.clear();
						retrieve_vertex_indices(trav, simplex_vertices);
						for (int j = 0; j < simplex_vertices.size(); ++j) {
							sstr << simplex_vertices[j] << (j == simplex_vertices.size() - 1 ? "" : " ");
						}
						sstr << endl;
						//
						trav = trav->next_circular_ptr;
					} while (trav != mIter->second);
				}
			}
		}
		//
		//ofile << sstr.rdbuf();
		ofile.write(sstr.str().c_str(), sstr.str().size());
		//
		ofile.close();
		sstr.clear();
	}
	else
	{
		std::cout << "Can NOT open file " << pFileName << std::endl;
		exit(0);
	}
}
#endif // _SIMPLICIAL_TREE_H_
