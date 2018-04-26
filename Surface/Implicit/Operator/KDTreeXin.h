//*****************************************************************************
//	KDTree.h: Header file for the KDTree class.
//	Xinlai Ni
//	10/25/2005
//*****************************************************************************


#ifndef	XNTOOLS_DS_KDTREE_H
#define	XNTOOLS_DS_KDTREE_H

#include <math.h>
#include <vector>

namespace xntools {
namespace ds {

//	templated KDTree class
//	Traits: the traits for the point, must have the following defined:
//		- point_type: type of the point to be searched
//		- coord_type: type of the coordinate of the points
//		- get_component: function object that returns a coord_type, the i-th component of the point
//		- dim: dimension of the point
template <typename Traits>
class KDTree
{
public:
	//	type of the point
	typedef		typename Traits::point_type		point_type;
	//	type of the coordinates
	typedef		typename Traits::coord_type		coord_type;
	//	type of the get_component function
	typedef		typename Traits::get_component	get_component;

	//	dimension of the point
	static const int dim = Traits::dim;

	//	construct the KD-tree from a vector of points
	void construct(const std::vector<point_type>& points);
	//	search for the nearest point from the query point within the optional radius r (r = -1 for global search)
	//	the index of the nearest point is returned, -1 for empty tree or none found within r
	//	upon return, 'dist' (if not NULL) is the distance to the nearest point
	int search_nearest(const coord_type query_point[], const std::vector<point_type>& points, coord_type* dist = NULL, coord_type r = -1) const;
	//	search for the n nearest points from the query point within the optional radius r (r = -1 for global search)
	//	the number of the points found is returned
	//	'index' should be a preallocated array of size n
	//	upon return, 'index' contains the indices of the qualified points, sorted from near to far
	//	'dist' contains the corresponding distance
	int search_nearest(const coord_type query_point[], const std::vector<point_type>& points, int n, int* index, coord_type* dist = NULL, coord_type r = -1) const;

private:
	//	recursive function to balance the tree
	void balance(std::vector<int>& holder, const std::vector<point_type>& points, int root, int start, int end, coord_type box_max[], coord_type box_min[]);
	//	partition the nodes around the root
	void partition(std::vector<int>& holder, const std::vector<point_type>& points, int start, int end, int root, int split_dim);
	//	recursive function to search the nearest point
	void search_nearest(const coord_type query_point[], const std::vector<point_type>& points, coord_type& r, int tree, int& result) const;
	//	recursive function to search the n nearest points
	void search_nearest(const coord_type query_point[], const std::vector<point_type>& points, int n, coord_type rad[], int tree, int& acc, int* heap) const;
	//	array of the point indices
	std::vector<int>		_nodes;
	//	splitting dimensions
	std::vector<int>		_splitters;
	//	function object
	static get_component	_component;
};

//	initialize the function object
template <typename Traits>
typename KDTree<Traits>::get_component KDTree<Traits>::_component;

//	construct the KD-tree from a vector of points
template <typename Traits>
void KDTree<Traits>::construct(const std::vector<point_type>& points)
{
	int num_points = (int)points.size();
	//	singleton leaf
	if( num_points <= 1 ) {
		_nodes.clear();
		_splitters.clear();
		_nodes.push_back(-1);
		_nodes.push_back(0);
		_splitters.push_back(-1);
		_splitters.push_back(0);
		return;
	}
	//	initialize the array
	_nodes.clear();
	_splitters.clear();
	_nodes.push_back(-1);
	_splitters.push_back(-1);
	//	temporary holder for the indices
	std::vector<int> holder;
	holder.push_back(-1);
	int i;
	for(i = 0; i < num_points; i++) {
		holder.push_back(i);
		_nodes.push_back(-1);
		_splitters.push_back(-1);
	}
	//	calculate the bounding box of the points
	int d;
	coord_type box_min[dim];
	coord_type box_max[dim];
	for(d = 0; d < dim; d++)
		box_min[d] = box_max[d] = _component(points[0], d);
	for(i = 1; i < num_points; i++) {
		for(d = 0; d < dim; d++) {
			coord_type x = _component(points[i], d);
			if( x < box_min[d] ) box_min[d] = x;
			if( x > box_max[d] ) box_max[d] = x;
		}
	}

	//	make a balanced binary tree
	if( num_points > 1 )
		balance(holder, points, 1, 1, num_points, box_max, box_min);
}

//	search for the nearest point from the query point within the optional radius r (r = -1 for global search)
template <typename Traits>
int KDTree<Traits>::search_nearest(const coord_type query_point[], const std::vector<point_type>& points, coord_type* dist, coord_type r) const
{
	//	is the tree empty?
	if( _nodes.size() <= 1 ) return -1;
	//	if r < 0, do global search, otherwise, make it square
	if( r > 0 ) r *= r;
	else r = -1;
	int result = -1;
	search_nearest(query_point, points, r, 1, result);
	if( dist != NULL )
		*dist = sqrt(r);
	return result;
}

//	search for the nearest n points from the query point
template <typename Traits>
int KDTree<Traits>::search_nearest(const coord_type query_point[], const std::vector<point_type>& points, int n, int* index, coord_type* dist, coord_type r) const
{
	//	is the tree empty?
	if( _nodes.size() <= 1 || index == NULL || n <= 0 ) return 0;
	int acc = 0;
	coord_type* rad = new coord_type[n + 1];
	if( r > 0 )	rad[0] = r * r;
	else rad[0] = -1;
	int* heap = new int[n + 1];
	search_nearest(query_point, points, n, rad, 1, acc, heap);
	//	export the heap to the sorted index array
	int size = acc;
	while( size > 0 ) {
		//	pop the current largest
		index[size - 1] = heap[1];
		if( dist != NULL )
			dist[size - 1] = sqrt(rad[1]);
		int bubble = 1;
		int current = 2;
		//	move the last item up
		while( current <= size ) {
			//	find greater child
			if( current < size && rad[current + 1] > rad[current] )
				current++;
			if( rad[size] >= rad[current] ) break;
			rad[bubble] = rad[current];
			heap[bubble] = heap[current];
			bubble = current;
			current <<= 1;
		}
		rad[bubble] = rad[size];
		heap[bubble] = heap[size];
		size--;
	}
	delete[] rad;
	delete[] heap;
	return acc;
}

//	recursive function to balance the tree
template <typename Traits>
void KDTree<Traits>::balance(std::vector<int>& holder, const std::vector<point_type>& points, int tree, int start, int end, coord_type box_max[], coord_type box_min[])
{
	//	find the subtree root index in [start, end]
	//	if n = end - start + 1 is the number of nodes in the range, the index of the root should keep the tree balanced
	//	let k = log n, the formula is:
	//	if n + 1 < 3*2^(k-1) then m = n - 2^(k-1) + 1 (left subtree is not full)
	//	else m = 2 ^ k (left subtree is full)
	int root;
	int n = end - start + 1;
	int k = 0;
	while( n != 1 && n != 3 ) { k++; n >>= 1; }
	if( n == 1 )
		root = end - (1 << (k - 1)) + 1;
	else
		root = start - 1 + (1 << (k + 1));

	//	find the splitting dimension, i.e., the longest dimension of the span
	int split_dim = 0;
	coord_type max_span = box_max[0] - box_min[0];
	for(int d = 1; d < dim; d++) {
		if( box_max[d] - box_min[d] > max_span ) {
			max_span = box_max[d] - box_min[d];
			split_dim = d;
		}
	}
	//	partition the holder around the root
	partition(holder, points, start, end, root, split_dim);
	//	copy to the tree
	_nodes[tree] = holder[root];
	_splitters[tree] = split_dim;
	//	recursively balance the left and right subtree
	if( root > start ) {
		if( start < root - 1 ) {
			//	the bounding box of the left subspace
			coord_type temp = box_max[split_dim];
			box_max[split_dim] = _component(points[_nodes[tree]], split_dim);
			//	balance the left subtree
			balance(holder, points, tree << 1, start, root - 1, box_max, box_min);
			box_max[split_dim] = temp;
		}
		else
			_nodes[tree << 1] = holder[start];
	}
	if( root < end ) {
		if( root + 1 < end ) {
			//	the bounding box of the right subspace
			coord_type temp = box_min[split_dim];
			box_min[split_dim] = _component(points[_nodes[tree]], split_dim);
			//	balance the right subtree
			balance(holder, points, (tree << 1) + 1, root + 1, end, box_max, box_min);
			box_min[split_dim] = temp;
		}
		else
			_nodes[(tree << 1) + 1] = holder[end];
	}
}

//	partition the nodes around the root
template <typename Traits>
void KDTree<Traits>::partition(std::vector<int>& holder, const std::vector<point_type>& points, int start, int end, int root, int split_dim)
{
	int left = start;
	int right = end;
	coord_type pivot;
	int i, j;
	while( left < right ) {
		pivot = _component(points[holder[right]], split_dim);
		i = left - 1;
		j = right;
		while( true ) {
			//	find the first point from the left that's larger than pivot
			while( _component(points[holder[++i]], split_dim) < pivot && i < right);
			//	find the first point from the right that's smaller than pivot
			while( _component(points[holder[--j]], split_dim) > pivot && j > left );
			if( i >= j ) break;
			//	swap i and j
			int temp = holder[i];
			holder[i] = holder[j];
			holder[j] = temp;
		}
		//	swap i and right, so that the pivor is on the i-th position
		int temp = holder[i];
		holder[i] = holder[right];
		holder[right] = temp;
		if( i >= root ) right = i - 1;
		if( i <= root ) left = i + 1;
	}
}

//	recursive function to search the nearest point
template <typename Traits>
void KDTree<Traits>::search_nearest(const coord_type query_point[], const std::vector<point_type>& points, coord_type& r, int tree, int& result) const
{
	coord_type dist;
	int split_dim = _splitters[tree];
	int num_nodes = (int)_nodes.size() - 1;
	int subtree;
	//	if 'tree' is not leaf, search its subtrees first
	if( tree <= (num_nodes >> 1) ) {
		//	distance from the query point to the splitting hyperplane
		dist = query_point[split_dim] - _component(points[_nodes[tree]], split_dim);
		if( dist > 0 ) {
			//	search the right subspace first
			if( (subtree = (tree << 1) | 1) <= num_nodes )
				search_nearest(query_point, points, r, subtree, result);
			//	search the left subspace only when splitting hyperplane is within the radius
			if( r < 0 || dist * dist < r )
				search_nearest(query_point, points, r, tree << 1, result);
		}
		else {
			//	search the left subspace first
			search_nearest(query_point, points, r, tree << 1, result);
			if( (r < 0 || dist * dist < r) && (subtree = (tree << 1) | 1) <= num_nodes )
				search_nearest(query_point, points, r, subtree, result);
		}
	}
	//	current node
	coord_type rad = 0;
	for(int d = 0; d < dim; d++) {
		dist = query_point[d] - _component(points[_nodes[tree]], d);
		rad += dist * dist;
	}
	//	if r < 0, doing global search
	if( r < 0 || rad < r ) {
		r = rad;
		result = _nodes[tree];
	}
}

//	recursive function to search the n nearest points
template <typename Traits>
void KDTree<Traits>::search_nearest(const coord_type query_point[], const std::vector<point_type>& points, int n, coord_type rad[], int tree, int& acc, int* heap) const
{
	coord_type dist;
	int split_dim = _splitters[tree];
	int num_nodes = (int)_nodes.size() - 1;
	int subtree;
	//	if 'tree' is not leaf, search its subtrees first
	if( tree <= (num_nodes >> 1) ) {
		//	distance from the query point to the splitting hyperplane
		dist = query_point[split_dim] - _component(points[_nodes[tree]], split_dim);
		if( dist > 0 ) {
			//	search the right subspace first
			if( (subtree = (tree << 1) | 1) <= num_nodes )
				search_nearest(query_point, points, n, rad, subtree, acc, heap);
			//	search the left subspace only when splitting hyperplane is within the radius
			if( rad[0] < 0 || dist * dist < rad[0] )
				search_nearest(query_point, points, n, rad, tree << 1, acc, heap);
		}
		else {
			//	search the left subspace first
			search_nearest(query_point, points, n, rad, tree << 1, acc, heap);
			if( ( rad[0] < 0 || dist * dist < rad[0]) && (subtree = (tree << 1) | 1) <= num_nodes )
				search_nearest(query_point, points, n, rad, subtree, acc, heap);
		}
	}

	//	current node
	coord_type r = 0;
	for(int d = 0; d < dim; d++) {
		dist = query_point[d] - _component(points[_nodes[tree]], d);
		r += dist * dist;
	}
	//	rad[0] is the initial radius request (when there are < n points found) or the current n-th smallest distance (when already found n points)
	int bubble, current;
	//	if rad[0] < 0, doing global search
	if( rad[0] < 0 || r < rad[0] ) {
		//	need to add a point to the max heap
		if( acc < n ) {
			bubble = ++acc;
			current = (bubble >> 1);
			//	search up the heap to find a bigger node
			while( current > 0 && rad[current] < r ) {
				rad[bubble] = rad[current];
				heap[bubble] = heap[current];
				bubble = current;
				current >>= 1;
			}
			rad[bubble] = r;
			heap[bubble] = _nodes[tree];
			//	if heap is full, update the sentinel
			if( acc == n )
				rad[0] = rad[1];
		}
		else {
			//	heap is already full, replace the top item, reorganize the heap
			bubble = 1;
			current = 2;
			while( current <= n ) {
				//	find greater child
				if( current < n && rad[current + 1] > rad[current] )
					current++;
				if( r >= rad[current] ) break;
				rad[bubble] = rad[current];
				heap[bubble] = heap[current];
				bubble = current;
				current <<= 1;
			}
			rad[bubble] = r;
			heap[bubble] = _nodes[tree];
			rad[0] = rad[1];
		}
	}
}

}	//	namespace ds
}	//	namespace xntools

#endif	//	XNTOOLS_DS_KDTREE_H
