/*
(c) 2015 Fengtao Fan, Dayu Shi
*/
#include "SimplexNode.h"

ListNode::ListNode(const ListNode &rhs)
{
	row = rhs.row;
	val = rhs.val;
	next = rhs.next;
	right = rhs.right;
	left = rhs.left;
}
ListNode& ListNode::operator=(const ListNode &rhs) {
	if (this != &rhs) {
		row = rhs.row;
		val = rhs.val;
		next = rhs.next;
		right = rhs.right;
		left = rhs.left;
	}
	return *this;
}
ListNodePtr AnnotationMatrix::DeepCopyAnnotationColumn(const ListNodePtr &head) {
	/*create the dummy node*/
	ListNodePtr new_head(boost::make_shared<ListNode>());
	std::unordered_map<ListNodePtr, TreeRootNodePtr, hash_ListNodePtr, equal_ListNodePtr>::iterator findIter = ann_mat.find(head);
	if (findIter != ann_mat.end())
	{
		if (head->next == head) {
			// zero annotation
			new_head->next = new_head;
			return new_head;
		}
		/*pointing to the dummy head*/
		ListNodePtr trav = head;
		/*point to the first element */
		trav = trav->next;
		//copy the first element
		new_head->next = boost::make_shared<ListNode>(trav->row, trav->val);
		new_head->next->next = new_head;
		//
		ListNodePtr tail = new_head->next;
		trav = trav->next;
		while (trav != head) {
			//copy the element
			tail->next = boost::make_shared<ListNode>(trav->row, trav->val);

			trav = trav->next;
			tail = tail->next;
			tail->next = new_head;
		}
	}
	else {
		// zero annotation;
		cout << "The annotation is not found" << endl;
		ListNodePtr trav = head;
		trav = trav->next;
		while (trav != head) {
			cout << trav->row << " ";
			trav = trav->next;
		}
		cout << endl;
		for (unordered_map<int, ListNodePtr>::iterator mIter = row_ptr.begin();
			mIter != row_ptr.end(); ++mIter) {
			cout << mIter->first << endl;
		}
		exit(0);
		return ListNodePtr();  // return nullptr
	}
	//
	return new_head;
}
ListNodePtr AnnotationMatrix::create_cocycle(TreeRootNodePtr &root, UnionFindDeletion &ufd, bool zero_elem){
	if (zero_elem) {
		ListNodePtr p(boost::make_shared<ListNode>());
		p->next = p;
		//
		Insert(p, root, ufd);
		// 
		return p;
	}
	// not the zero element

	// generate the dummy node first
	ListNodePtr p(boost::make_shared<ListNode>());
	p->next = boost::make_shared<ListNode>(timeStamp, 1);
	timeStamp += 1;
	p->next->next = p; // make it circular
	// link root to p
	Insert(p, root, ufd);
	//
	return p;
}
void AnnotationMatrix::Insert(ListNodePtr &ptr, const TreeRootNodePtr root, UnionFindDeletion &ufd)
{
	if (!search(ptr))
	{
		// insert the list
		ann_mat[ptr] = root;
		// link p to root 
		root->attribute = ptr;
		// link into row lists
		// move to the first element
		ListNodePtr p(ptr->next);
		while (p != ptr)
		{
			insert_into_doubly_linked_list(p->row, p);
			//
			p = p->next;
		}
	}
	else {
		// merge the two clusters
		TreeRootNodePtr root1 = ann_mat[ptr];
		TreeRootNodePtr root2 = root;
		root1->attribute.reset(); //unlink tree with annotation
		TreeRootNodePtr newRoot = ufd.Union(root1, root2); // merge two trees
		// update cluster associated with the annotation 
		newRoot->attribute = ann_mat.find(ptr)->first;
		ann_mat[ptr] = newRoot;
		//clear listNode "ptr"
		ptr->next.reset();
		ptr = newRoot->attribute;
	}
	return;
}

void AnnotationMatrix::clearNode(ListNodePtr &ptr, bool bUpdatePers)
{//

	std::unordered_map<ListNodePtr, TreeRootNodePtr, hash_ListNodePtr, equal_ListNodePtr>::iterator findIter = ann_mat.find(ptr);
	if (findIter != ann_mat.end())
	{
		// delete from the row list
		// move to the first element
		ListNodePtr p(findIter->first->next);
		while (p != findIter->first)
		{
			delete_from_doubly_linked_list(p->row, p, bUpdatePers);
			p = p->next;
		}
		// delete from the list
		TreeRootNodePtr ret = findIter->second;
		ret->attribute.reset(); // unlink the root with annotation
		//
		ann_mat.erase(findIter);
		//break cycle
		p->next.reset();
		return;
	}
	//
	return;
}

TreeRootNodePtr AnnotationMatrix::Delete(ListNodePtr &ptr)
{//

	std::unordered_map<ListNodePtr, TreeRootNodePtr, hash_ListNodePtr, equal_ListNodePtr>::iterator findIter = ann_mat.find(ptr);
	if (findIter != ann_mat.end())
	{
		// delete from the row list
		// move to the first element
		ListNodePtr p(findIter->first->next);
		while (p != findIter->first)
		{
			//don't update persistences
			delete_from_doubly_linked_list(p->row, p, false);
			//
			p = p->next;
		}
		// delete from the list
		TreeRootNodePtr ret = findIter->second;
		ret->attribute.reset(); // unlink the root with annotation
		//
		ann_mat.erase(findIter);
		return ret;
	}
	//
	return boost::make_shared<TreeRootNode>();
}
int AnnotationMatrix::sum_two_annotation_with_changed_dst(ListNodePtr & out_dst, ListNodePtr & in_src) {
	// change dst and keep src unchanged 
	if (!in_src) {
		// in_src is empty
		return lowest_one(out_dst);
	}
	if (!out_dst) {
		// out_dst is empty but in_src is nonzero
		out_dst = DeepCopyAnnotationColumn(in_src);
		return lowest_one(out_dst);
	}
	ListNodePtr dst_prev(out_dst);
	ListNodePtr dst_trav(out_dst->next);
	ListNodePtr src_trav(in_src->next);

	while (dst_trav != out_dst && src_trav != in_src) {
		if (dst_trav->row == src_trav->row) {
			// thest two should be canceled
			dst_prev->next = dst_trav->next;
			//
			dst_trav = dst_trav->next;
			src_trav = src_trav->next;
		}
		else {
			if (dst_trav->row > src_trav->row) {
				ListNodePtr temp(boost::make_shared<ListNode>(src_trav->row, 1));
				dst_prev->next = temp;
				temp->next = dst_trav;

				dst_prev = temp;
				src_trav = src_trav->next;
			}
			else {
				dst_prev = dst_trav;
				dst_trav = dst_trav->next;
			}
		}
	}
	while (src_trav != in_src) {
		dst_prev->next = boost::make_shared<ListNode>(src_trav->row, 1);
		dst_prev = dst_prev->next;
		dst_prev->next = out_dst;

		src_trav = src_trav->next;
	}
	while (dst_prev->next != out_dst) {
		dst_prev = dst_prev->next;
	}
	return (out_dst->next != out_dst ? dst_prev->row : -1);
}
int AnnotationMatrix::lowest_one(ListNodePtr & head) {
	if (head) {
		ListNodePtr trav(head->next);
		if (trav == head) {
			return -1;
		}
		while (trav->next != head) {
			trav = trav->next;
		}
		return trav->row;
	}
	return -1;
}
void AnnotationMatrix::kill_cocycle_last_nonzero_bit(const int u, ListNodePtr &ext_src, UnionFindDeletion &ufd){
	// add every annotation with nonzero u-th bit by exteranl annotation ext_src
	while (row_ptr.find(u) != row_ptr.end()) {
		ListNodePtr mid_node = row_ptr[u];
		// the dummy head has row of -1
		while (mid_node->row != -1) {
			mid_node = mid_node->next;
		}
		// mide_node is the dummy head now
		// remove it from the annotation matrix
		TreeRootNodePtr x = Delete(mid_node);
		sum_two_annotation_with_changed_dst(mid_node, ext_src);
		// the simplex has new annotation 
		// insert the simplex back to the forest
		Insert(mid_node, x, ufd);
	}
	return;
}