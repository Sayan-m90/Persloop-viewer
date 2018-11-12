/*
(c) 2015 Fengtao Fan, Dayu Shi
*/
#include "SimplexNode.h"  
#include <vector>
#include <queue>
#include <stack>

using namespace std;
TreeRootNodePtr UnionFindDeletion::MakeSet(ElementNodePtr &elem)
{
	//std::unordered_map<int, ElementNodePtr>::iterator findIter = elemSet.find(elem->value);
	//
	//if (findIter == elemSet.end())
	//{ 
	//Create a new tree root node
	TreeRootNodePtr root(boost::make_shared<TreeRootNode>());
	//Set double-links between the element node and the tree node
	link_tree_node_and_element(root, elem);

	/*Initializing tree root data*/
	//Store the tree root in the forest
	//forest[elem->value] = root;
	//Store the element value in the hash of elments
	//elemSet[elem->value] = elem;

	/*Initializing tree node*/
	//Set parent
	root->parent = root;
	//Set rank
	root->rank = 0;
	root->cListSize = 0;
	//Set dfsList
	root->dfsNext = root;
	root->dfsPrev = root;

	//}
	//else
	//{
	//	std::cout << "The element [" << elem->value << "] exists" << std::endl; 
	//	std::cout << "MAKESET operation is NOT performed" << std::endl;
	//}
	//
	//return forest[elem->value];
	return root;
}
void UnionFindDeletion::insert_into_nlist(TreeNodePtr a, TreeRootNodePtr &T)
{
	if (T->nListHead) {
		TreeNodePtr tail = T->nListHead->nPrev;
		// head-->a-->tail
		// head<--a<--tail
		T->nListHead->nPrev = a;
		a->nPrev = tail;
		tail->nNext = a;
		a->nNext = T->nListHead;

		T->nListHead = a;
	}
	else {
		T->nListHead = a;
		a->nNext = a;
		a->nPrev = a;
	}
	return;
}
void UnionFindDeletion::insert_into_clist_and_dfs_list(TreeNodePtr a, TreeNodePtr parent)
{
	//Add to the children list 
	if (parent->cListHead)
	{// insert to the begining of the list
		TreeNodePtr cListTail = parent->cListHead->cPrev;
		// head-->a-->tail
		// head<--a<--tail
		parent->cListHead->cPrev = a;
		a->cPrev = cListTail;
		cListTail->cNext = a;
		a->cNext = parent->cListHead;

		parent->cListHead = a;
	}
	else
	{//The first element of cListHead
		parent->cListHead = a;
		a->cNext = a;
		a->cPrev = a;
	}
	//Add to dfsList of T_keep 
	TreeNodePtr parent_dfsNext = parent->dfsNext;
	TreeNodePtr child_dfsTail = a->dfsPrev;

	// parent-->a----a_tail->parent_next
	// parent<--a----a_tail<--parent_next
	a->dfsPrev = parent;
	parent->dfsNext = a;

	child_dfsTail->dfsNext = parent_dfsNext;
	parent_dfsNext->dfsPrev = child_dfsTail;

	//Update counters
	parent->cListSize++;

	return;
}
void UnionFindDeletion::LinkSingleNodeToTreeRoot(TreeNodePtr a, TreeRootNodePtr &T)
{
	//Link to the root of T_keep
	a->parent = T;
	//Set its rank as zero
	a->rank = 0;
	//Set cListHead and its size
	a->cListHead.reset();
	a->cListSize = 0;
	//prepare for insert into dfs_list
	a->dfsNext = a;
	a->dfsPrev = a;
	//
	insert_into_clist_and_dfs_list(a, T);

	return;

}
TreeRootNodePtr UnionFindDeletion::Union(TreeRootNodePtr & a, TreeRootNodePtr & b)
{
	if (a == b) {
		return a;
	}
	TreeRootNodePtr dead = a; // dead refer to the tree to be destroyed
	TreeRootNodePtr alive = b;
	//
	int a_size = 4;
	if (!a->nListHead) {
		a_size = a->cListSize + 1;
	}
	int b_size = 4;
	if (!b->nListHead) {
		b_size = b->cListSize + 1;
	}

	if (a_size < 4 || b_size < 4)
	{ // one of them is reduced
		if (a_size >= 4)
		{
			dead = b;
			alive = a;
		}
		// remove it from the forest
		//forest.erase(dead->elem->value);
		// 		
		TreeNodePtr trav = dead->cListHead;
		while (trav) {
			if (trav->cNext == dead->cListHead) {
				dead->cListHead.reset();
			}
			else {
				dead->cListHead = trav->cNext;
				trav->cNext->cPrev = trav->cPrev;
				trav->cPrev->cNext = trav->cNext;
			}
			LinkSingleNodeToTreeRoot(trav, alive);
			trav = dead->cListHead;
		}
		// 
		LinkSingleNodeToTreeRoot(dead, alive);
		// update the root's rank
		alive->rank = std::max(alive->rank, 1);
		//
	}
	else
	{
		if (a->rank > b->rank)
		{
			dead = b;
			alive = a;
		}
		dead->parent = alive;
		// update root's rank
		if (dead->rank == alive->rank) {
			++alive->rank;
		}
		//forest.erase(dead->elem->value);
		// set rank of T_keep
		insert_into_nlist(dead, alive);
		insert_into_clist_and_dfs_list(dead, alive);
		// 
		dead->nListHead.reset();
	}
	return alive;
}
void UnionFindDeletion::InsertIntoCList(TreeNodePtr p, TreeNodePtr ref, TreeNodePtr a, bool right)
{
	//
	a->parent = p;
	if (p->cListHead)
	{
		if (right)
		{
			// next--ref---a--prev
			TreeNodePtr prev = ref->cPrev;
			prev->cNext = a;
			a->cPrev = prev;
			ref->cPrev = a;
			a->cNext = ref;
			if (ref == p->cListHead) {
				p->cListHead = a;
			}
		}
		else
		{// insert to the left of ref
			//next--a--ref--prev
			TreeNodePtr next = ref->cNext;
			next->cPrev = a;
			a->cNext = next;
			ref->cNext = a;
			a->cPrev = ref;
		}
	}
	else
	{
		p->cListHead = a;
		a->cNext = a;
		a->cPrev = a;
	}
	// 
	p->cListSize++;
	//

	return;
}
void UnionFindDeletion::RemoveFromCList(TreeNodePtr a)
{// a is not the root
	TreeNodePtr p = a->parent;
	a->parent.reset();
	if (p->cListSize == 1)
	{
		p->cListHead.reset();
		//
		a->cNext.reset();
		a->cPrev.reset();
		// 
	}
	else
	{
		TreeNodePtr prev = a->cPrev;
		TreeNodePtr next = a->cNext;
		//
		prev->cNext = next;
		next->cPrev = prev;
		if (a == p->cListHead) {
			p->cListHead = next;
		}
		//
		a->cNext.reset();
		a->cPrev.reset();
	}
	// 
	p->cListSize--;
	//

	return;
}
void UnionFindDeletion::RemoveFromNList(TreeNodePtr a)
{
	TreeRootNodePtr root = find_root(a);//forest[a->parent->elem->value];
	//
	if (root->nListHead->nNext == root->nListHead)
	{// only one element
		root->nListHead.reset();
	}
	else
	{
		TreeNodePtr head = root->nListHead;
		TreeNodePtr prev = a->nPrev;
		TreeNodePtr next = a->nNext;
		//
		prev->nNext = next;
		next->nPrev = prev;
		if (a == head) {
			root->nListHead = next;
		}
	}
	return;
}
void UnionFindDeletion::InsertIntoNList(TreeRootNodePtr p, TreeNodePtr ref, TreeNodePtr a) {
	if (p->parent == p)
	{// p is also the root
		if (a->cListSize > 0)
		{// insert it into the NList of the root
			TreeNodePtr next = ref->nNext;
			next->nPrev = a;
			a->nNext = next;
			ref->nNext = a;
			a->nPrev = ref;
		}
	}
	return;
}
void UnionFindDeletion::Relink(TreeNodePtr a)
{// p(p(a)) is not the root
	TreeNodePtr p = a->parent;
	if (a->cNext != p->cListHead)
	{// a has a left sibling
		TreeNodePtr L = a->cNext;
		// remove a from CList of p
		RemoveFromCList(a);
		// add it into 
		InsertIntoCList(p->parent, p, a, true);
		if (p->parent->parent == p->parent && a->cListSize > 0) {
			TreeRootNodePtr root = find_root(p->parent);//forest[p->parent->elem->value];
			InsertIntoNList(root, p, a);
		}
		// remove the segment from DFSList
		TreeNodePtr tail = L->dfsPrev;
		// disconnect [a, tail]
		L->dfsPrev = a->dfsPrev;
		a->dfsPrev->dfsNext = L;
		// insert [a, tail] before p in DFSList
		tail->dfsNext = p;
		a->dfsPrev = p->dfsPrev;
		p->dfsPrev->dfsNext = a;
		p->dfsPrev = tail;
	}
	else
	{
		// remove a from CList of p
		RemoveFromCList(a);
		// add it into 
		InsertIntoCList(p->parent, p, a, false);
		if (p->parent->parent == p->parent && a->cListSize > 0) {
			TreeRootNodePtr root = find_root(p->parent);//forest[p->parent->elem->value];
			InsertIntoNList(root, p, a);
		}
	}
	// if p is leaf 
	if (p->cListSize == 0)
	{
		p->rank = 0;
		//
		if (p->parent->parent == p->parent)
		{// p->parent is the root
			RemoveFromNList(p);
			// if NList is empty, the tree is reduced
			TreeRootNodePtr root = find_root(p->parent);//forest[p->parent->elem->value];
			if (!root->nListHead) {
				p->parent->rank = 1;
			}
		}

	}
}
TreeRootNodePtr UnionFindDeletion::Find(TreeNodePtr a)
{
	while (a->parent->parent != a->parent)
	{
		TreeNodePtr p = a->parent;
		Relink(a);
		if (p->cListSize == 2)
		{
			while (p->cListSize > 0)
			{
				Relink(p->cListHead);
			}
		}
		a = p;
	}
	TreeRootNodePtr root = find_root(a->parent);//forest[a->parent->elem->value];
	return root;
}
void UnionFindDeletion::RemoveFromDFSList(TreeNodePtr a)
{// a is not the root unless the tree contains only the root
	if (a->parent == a)
	{
		a->dfsNext.reset();
		a->dfsPrev.reset();
	}
	else
	{
		TreeNodePtr next = a->dfsNext;
		TreeNodePtr prev = a->dfsPrev;
		//
		next->dfsPrev = prev;
		prev->dfsNext = next;
		//
		a->dfsNext.reset();
		a->dfsPrev.reset();
	}
	return;
}
void UnionFindDeletion::DeleteFromReducedTree(TreeNodePtr a)
{// the tree is reduced
	TreeNodePtr to_be_deleted = a;
	if (a->parent == a && a->cListSize > 0)
	{// a is the root 
		to_be_deleted = a->dfsPrev;
		//switch conents
		ElementNodePtr keepElem = to_be_deleted->elem;
		ElementNodePtr rootElem = a->elem;
		//
		unlink_tree_node_and_element(a, rootElem);
		unlink_tree_node_and_element(to_be_deleted, keepElem);
		//
		link_tree_node_and_element(a, keepElem);
		link_tree_node_and_element(to_be_deleted, rootElem);
		//element and treenodes doubly linked 
		//update forest
		//forest[a->elem->value] = forest[to_be_deleted->elem->value];
		//forest.erase(to_be_deleted->elem->value);
	}
	//
	if (to_be_deleted->parent != to_be_deleted) {
		RemoveFromCList(to_be_deleted);
	}
	//
	RemoveFromDFSList(to_be_deleted);
	//
	to_be_deleted->elem->tree_node.reset();
	//elemSet.erase(to_be_deleted->elem->value);
	//
	//if (to_be_deleted->parent == to_be_deleted)
	//{// this is the root in the tree
	// destroy this class 
	//
	//forest.erase(to_be_deleted->elem->value);
	// 
	//} 
	//
	//break cycles
	to_be_deleted->elem.reset();
	to_be_deleted->parent.reset();
	to_be_deleted->cListHead.reset();
	to_be_deleted->cPrev.reset();
	to_be_deleted->cNext.reset();
	to_be_deleted->nPrev.reset();
	to_be_deleted->nNext.reset();
	to_be_deleted->dfsPrev.reset();
	to_be_deleted->dfsNext.reset();

	return;
}
void UnionFindDeletion::LocalRebuild(TreeNodePtr &p)
{
	if (p->parent == p)
	{// p is the root
		// relink the three left most children of a non-leaf node c.
		TreeRootNodePtr root = find_root(p);//forest[p->elem->value];
		TreeNodePtr c = root->nListHead->nPrev;
		//
		if (c->cListSize < 6) {
			while (c->cListSize > 0) {
				Relink(c->cListHead->cPrev);
			}
		}
		else {
			for (int i = 0; i < 3; i++) {
				Relink(c->cListHead->cPrev);
			}
		}
	}
	else {
		if (p->cListSize < 5) {
			while (p->cListSize > 0) {
				Relink(p->cListHead->cPrev);
			}
		}
		else {
			Relink(p->cListHead->cPrev);
			Relink(p->cListHead->cPrev);
		}
	}
	return;
}
void UnionFindDeletion::Delete(TreeNodePtr a)
{
	bool bReducedTree = true;
	if (a->parent == a)
	{// a is the root
		TreeRootNodePtr root = find_root(a);//forest[a->elem->value]; 
		if (root->nListHead){
			bReducedTree = false;
		}
	}
	else
	{
		if (a->parent->parent == a->parent)
		{// a->parent is the root
			TreeRootNodePtr root = find_root(a);//forest[a->parent->elem->value];
			if (root->nListHead) {
				bReducedTree = false;
			}
		}
		else
		{// the height is at least 2
			bReducedTree = false;
		}
	}
	if (bReducedTree)
	{
		DeleteFromReducedTree(a);
	}
	else
	{
		TreeNodePtr del = a;
		// find leaf
		if (a->cListSize > 0)
		{
			if (a->parent == a)
			{// a is the root
				del = a->dfsPrev;
				//switch conents
				ElementNodePtr keepElem = del->elem;
				ElementNodePtr delElem = a->elem;

				unlink_tree_node_and_element(a, delElem);
				unlink_tree_node_and_element(del, keepElem);

				link_tree_node_and_element(a, keepElem);
				link_tree_node_and_element(del, delElem);
				//
				// update forest
				//forest[a->elem->value] = forest[del->elem->value];
				//forest.erase(del->elem->value);
			}
			else
			{
				if (a->cNext == a->parent->cListHead)
				{// no left sibling
					del = a->dfsPrev;
				}
				else
				{
					del = a->cNext->dfsPrev;
				}
				//switch contents
				ElementNodePtr keepElem = del->elem;
				ElementNodePtr delElem = a->elem;

				unlink_tree_node_and_element(a, delElem);
				unlink_tree_node_and_element(del, keepElem);

				link_tree_node_and_element(a, keepElem);
				link_tree_node_and_element(del, delElem);
			}
		}
		//
		TreeNodePtr delParent = del->parent;
		// remove from CList
		RemoveFromCList(del);
		// remove from DFSList
		RemoveFromDFSList(del);
		//
		del->elem->tree_node.reset();
		//elemSet.erase(del->elem->value);
		//  
		LocalRebuild(delParent);
		//break cycles
		del->elem.reset();
		del->parent.reset();
		del->cListHead.reset();
		del->cPrev.reset();
		del->cNext.reset();
		del->nPrev.reset();
		del->nNext.reset();
		del->dfsPrev.reset();
		del->dfsNext.reset();
	}
	return;
}
