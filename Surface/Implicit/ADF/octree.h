#ifndef SPACEOCTREE_INCLUDED // -*- C++ -*-
#define SPACEOCTREE_INCLUDED

#ifdef _MSC_VER
#pragma once
#endif

/************************************************************************
   Re-usable octree data structure and related support classes.
  
   Child nodes are ordered as follows; the first node is in the all-negative
   quadrant; other nodes are indexed by their relation to it, looking at the
   binary representation of their child number; ie: if the first bit is set,
   it's the plus x, otherwise negative x; if the second bit is set, then
   it's the plus y, otherwise negative y; and so forth.
  
   Points exactly on the dividing plane are counted as part of the negative
   side.

   $Id: octree.h,v 1.1.1.1 2004/02/04 03:19:51 wensu Exp $
 ************************************************************************/

// #include <gfx/gfx.h>
// #include <gfx/gmVector3.h>
// #include <slim/heap.h>
#include "heap.h"
#include "pool.h"
#include "libgm/gm.h"

// Forward declarations.

class OctreeNode;
class OctreeTraversalNode;
class OctreeTraversalNodeFactory;
class OctreeIterator;
class OctreeFilteredIterator;
class OctreeOrderedIterator;


/**
 * The base octree class.
 */
class Octree
{
public:
    Octree(const gmVector3 &c, double l) : center(c)
    {
        length = l;
        max_actual_depth = 0;
    }

    // Note: subclasses are responsible for allocating and destroying
    // the root.

    virtual ~Octree() {}


    int get_child_for_point(const gmVector3 &p, const gmVector3 &center);
    bool get_leaf_for_point(const gmVector3 &p, OctreeTraversalNode *nodep, OctreeTraversalNode *scratchp);

    bool traverse(OctreeIterator *ip);
    double traverse_in_order(OctreeOrderedIterator *ip);
    bool traverse_with_filter(OctreeFilteredIterator *ip);


protected:
    // These are required since traversal nodes may carry data that the
    // octree doesn't know about, which needs to be propogated or copied
    // during traversal.

    virtual void copy_traversal_node(OctreeTraversalNode *dp, OctreeTraversalNode *sp);

    virtual void init_traversal_node(OctreeTraversalNode *np, OctreeTraversalNode *pp, int child) 
    {	return;}


    bool traverse_node(OctreeTraversalNode *np, OctreeIterator *ip);
    bool traverse_node_with_filter(OctreeTraversalNode *np, OctreeFilteredIterator *ip );

    // For octrees which initialize the octree data members themselves.
    Octree() { max_actual_depth = 0; }

public:
    gmVector3 center;
    double length;
    int max_actual_depth;                      // Root is at depth 0.
    OctreeNode *rootp;
    OctreeTraversalNodeFactory *node_factoryp;
};


/**
 * Offsets of the child centers from the center of the parent node.
 * eg, child 0 is offset -1, -1, -1 times half the size of child's node.
 */
const gmVector3 octree_child_center_offset[] =
{
    gmVector3(-1, -1, -1),
    gmVector3( 1, -1, -1),
    gmVector3(-1,  1, -1),
    gmVector3( 1,  1, -1),
    gmVector3(-1, -1,  1),
    gmVector3( 1, -1,  1),
    gmVector3(-1,  1,  1),
    gmVector3( 1,  1,  1)
};


/**
 * In certain applications, integer versions are more useful.
 * (the cost of casting can be high).
 */
const int octree_child_center_offset_ints[][3] =
{
    {-1, -1, -1},
    { 1, -1, -1},
    {-1,  1, -1},
    { 1,  1, -1},
    {-1, -1,  1},
    { 1, -1,  1},
    {-1,  1,  1},
    { 1,  1,  1}
};


/**
 * Base class for octree nodes.  All subclasses must use nodes inheriting
 * this class.
 */
class OctreeNode
{
public:
    OctreeNode() { childrenpp = NULL; parentp = NULL; }
    ~OctreeNode()
    {
        parentp = NULL;
        if ( childrenpp != NULL )
            for( int i=0; i<8; i++ )
                delete childrenpp[i];
        delete[] childrenpp;
    }


public:
    OctreeNode *  parentp;
    OctreeNode **  childrenpp;
};


/**
 * Used for temporary storage of octree nodes during traversal or insertion.
 * Heapable to allow ordered traversals.
 */
class OctreeTraversalNode : public MxHeapable
{
public:
    gmVector3 center;
    double length;
    int depth;
    OctreeNode * nodep;
};


/**
 * A factory for traversal nodes.  This is necessary since different Octree
 * subclasses may require different OctreeTraversalNode subclasses.
 */
class OctreeTraversalNodeFactory
{
public:
    virtual OctreeTraversalNode *  new_traversal_node() = NULL;
    virtual void free_traversal_node(OctreeTraversalNode *np) = NULL;
    virtual int get_allocated_count() = NULL;
};


/**
 * An object used for iterating or traversing over an octree.
 * The octree traversal method calls the iterate method of the iterator
 * for every node in the octree.  If the iterator ever returns false,
 * traversal is halted, and false is returned by the traversal method.
 * If wants_internal_nodes is false, the iterator will only receive
 * leaf nodes from the traversal method.
 */
class OctreeIterator
{
public:
    virtual bool iterate(OctreeTraversalNode *np) { return false; }

public:
    bool wants_internal_nodes;
};


/**
 * An object used for traversing octree nodes, while pruning out certain
 * subtrees of the octree.  The pruning is controlled by means of the
 * return value of the filtered_iterate method.
 * Note that even if the filter marks certain children as candidates for 
 * traversal, they may not exist in the octree, and therefore would not 
 * be traversed.  However, the iterator is free to create such children 
 * before returning.  The traversal is guaranteed to call the iterator 
 * for a parent before looking at the children of the parent, so this 
 * can be safely used for building the octree.
 * Note that internal nodes are always delivered to filtering iterators.
 * Since this traversal is depth-first, it is guaranteed that only
 * max_actual_depth + 2 traversal nodes will be required.
 * abortp may be set to true by the iterator to prematurely abort iteration.
 */
class OctreeFilteredIterator
{
public:
    virtual unsigned char filtered_iterate(OctreeTraversalNode *np, bool *abortp) = NULL;
};


/**
 * An object used for iterating or traversing over an octree in a set order.
 * The particular iterator defines the order in which traversal will occur,
 * by assigning weights to each node before it is traversed.
 * The traversal keeps a notion of 'maximum current weight', and ends 
 * traversal when the weight of the next node exceeds this maximum.  
 * It is possible that all leaves of the octree may be enqueued at the
 * same time; no better bound on the number of traversal nodes required
 * is possible in the general case.
 * The traversal is guaranteed to call the iterator for a parent before
 * looking at the children of the parent, so this can be safely used for
 * building an octree.
 */
class OctreeOrderedIterator
{
public:
    virtual double ordered_iterate(OctreeTraversalNode *np) = NULL;
    virtual double calculate_weight(OctreeTraversalNode *np) = NULL;
    
public:
    MxHeap nodepq;
    bool wants_internal_nodes;
};


/**
 * An iterator used for pruning unused leaves from an octree.
 * Any child which is unused is deleted.  Subclasses should define
 * is_used appropriately.
 */
class OctreePruner : public OctreeIterator
{
public:
    OctreePruner() { wants_internal_nodes = false; }

    virtual bool is_used(OctreeTraversalNode *np) { return false; }
    virtual void prune(Octree *  op) { octreep = op; op->traverse(this); }

    virtual bool iterate(
        OctreeTraversalNode *  np )
    {
        if ( !this->is_used(np) )
        {
            for( int i=0; i<8; i++ )
                if ( np->nodep->parentp->childrenpp[i] == np->nodep )
                {
                    np->nodep->parentp->childrenpp[i] = NULL;
                    break;
                }

            delete np->nodep;
        }

        return true;
    }

public:
    Octree * octreep;
};


/**
 * A statistics gatherer for octrees.
 */
class OctreeAnalyzer : public OctreeIterator
{
public:
    OctreeAnalyzer()
    {
        octreep = NULL;
        leaf_countp = NULL;
        internal_node_countp = NULL;
        wants_internal_nodes = true;
    }

    ~OctreeAnalyzer()
    {
        if ( leaf_countp != NULL )
            delete[] leaf_countp;
        if ( internal_node_countp != NULL )
            delete[] internal_node_countp;
    }

    void collect_statistics(Octree *  op); 
    virtual void report_statistics();
    virtual void report_statistics_at_depth(int  d);
    virtual void reset();

    virtual bool iterate(OctreeTraversalNode *np);

public:
    Octree *  octreep;
    int total_nodes;

    // These are collected on a per-level basis.

    int * leaf_countp;
    int * internal_node_countp;
};


/**
 * A base factory for traversal nodes.  This is an implementation, supplied
 * as a convenience for subclasses.  (not necessary that it be used).
 * As the name implies, it uses a pool to manage traversal node allocations.
 * The class T must be an OctreeTraversalNode subclass.
 */
template<class T>
class OctreeTraversalNodePool : public OctreeTraversalNodeFactory
{
public:
    OctreeTraversalNodePool(int chunk_size, int num_chunks, bool is_growable = true)
        : node_pool(chunk_size, num_chunks, is_growable) {}


    virtual OctreeTraversalNode *  new_traversal_node() 
    {
        return node_pool.alloc();
    }

    virtual void free_traversal_node(OctreeTraversalNode *np)
    {
        node_pool.dealloc( (T *) np );
    }

    virtual int get_allocated_count()
    {
        return node_pool.get_allocated_count();
    }


public:
    pool<T> node_pool;
};

#endif
