/***********************************************************************
   Implementation of the octree class.
 ***********************************************************************/

#include "octree.h"

int Octree::get_child_for_point(
    const gmVector3& p,
    const gmVector3& c )
{
    int child = 0;
    gmVector3 tmp = p - c;
    if ( tmp[0] > 0 )
        child += 1;
    if ( tmp[1] > 0 )
        child += 2;
    if ( tmp[2] > 0 )
        child += 4;

    return child;
}


void Octree::copy_traversal_node(
    OctreeTraversalNode *  dp,
    OctreeTraversalNode *  sp )
{
    dp->center = sp->center;
    dp->length = sp->length;
    dp->depth = sp->depth;
    dp->nodep = sp->nodep;
}


bool Octree::get_leaf_for_point(
    const gmVector3&            p,
    OctreeTraversalNode *  np,
    OctreeTraversalNode *  n2p )
{
    // Start at the top, and descend til we find a leaf.

    np->center = center;
    np->length = length;
    np->depth = 0;
    np->nodep = rootp;
    this->init_traversal_node(np, np, -1);

    // Determine if the point is actually in the octree.

    gmVector3 tmp = p - center;
    tmp[0] = fabs(tmp[0]);
    tmp[1] = fabs(tmp[1]);
    tmp[2] = fabs(tmp[2]);
    bool ret = true;
    if ( (tmp[0] > length)
         || (tmp[1] > length)
         || (tmp[2] > length) )
        ret = false;

    OctreeTraversalNode *  tnpp[2];
    tnpp[0] = np;
    tnpp[1] = n2p;
    int curr = 0;
    int next = 1;

    while( tnpp[curr]->nodep->childrenpp != NULL )
    {
        // Determine which octant the point is in.

        int child = this->get_child_for_point(p, tnpp[curr]->center);

        // If the child doesn't exist in the tree, return the parent.

        if ( tnpp[curr]->nodep->childrenpp[child] == NULL )
        {
            ret = false;
            break;
        }

        // Adjust the current node appropriately for the next child.

        tnpp[next]->length = 0.5 * tnpp[curr]->length;
        tnpp[next]->center = 
            tnpp[curr]->center + 
            (tnpp[next]->length * octree_child_center_offset[child]);
        tnpp[next]->depth = tnpp[curr]->depth + 1;
        tnpp[next]->nodep = tnpp[curr]->nodep->childrenpp[child];
        this->init_traversal_node(tnpp[next], tnpp[curr], child);

        curr = (curr + 1) & 1;
        next = (next + 1) & 1;
    }

    if ( curr != 0 )
        this->copy_traversal_node(tnpp[0], tnpp[1]);

    return ret;
}


bool Octree::traverse_node(
    OctreeTraversalNode *  np,
    OctreeIterator *       ip )
{
    // If it's a leaf node, or the iterator wants internal nodes, iterate it.

    if ( np->nodep->childrenpp == NULL || ip->wants_internal_nodes )
    {
        if ( !ip->iterate(np) )
            return false;
    }

    // If it's (still) a leaf, there's nothing left to do.

    if ( np->nodep->childrenpp == NULL )
        return true;

    // Recursively traverse the children.

    OctreeTraversalNode *  nextp = node_factoryp->new_traversal_node();
    nextp->depth = np->depth + 1;
    nextp->length = np->length * 0.5;
    for( int j=0; j<8; j++ )
    {
        nextp->nodep = np->nodep->childrenpp[j];
        if ( nextp->nodep == NULL )
            continue;

        nextp->center = np->center + 
                        nextp->length*octree_child_center_offset[j];
        this->init_traversal_node(nextp, np, j);
        if ( !this->traverse_node(nextp, ip) )
            return false;
    }
    node_factoryp->free_traversal_node(nextp);

    return true;
}

    
bool Octree::traverse(
    OctreeIterator *  ip )
{
    OctreeTraversalNode *  np = node_factoryp->new_traversal_node();
    np->center = center;
    np->length = length;
    np->depth = 0;
    np->nodep = rootp;
    this->init_traversal_node(np, np, -1);
    
    bool ret = this->traverse_node(np, ip);
    node_factoryp->free_traversal_node(np);
    return ret;
}


bool Octree::traverse_node_with_filter(
    OctreeTraversalNode *     np,
    OctreeFilteredIterator *  ip )
{
    bool abort = false;
    unsigned char children = ip->filtered_iterate(np, &abort);
    if ( abort )
        return false;

    if ( (children == 0) || (np->nodep->childrenpp == NULL) )
        return true;

    OctreeTraversalNode *  nextp = node_factoryp->new_traversal_node();
    nextp->depth = np->depth + 1;
    nextp->length = np->length * 0.5;
    for( int j=0; j<8; j++ )
    {
        nextp->nodep = np->nodep->childrenpp[j];
        if ( (nextp->nodep == NULL) || ((1 << j) & children) == 0 )
            continue;

        nextp->center = np->center + 
                        nextp->length*octree_child_center_offset[j];
        this->init_traversal_node(nextp, np, j);
        if ( !this->traverse_node_with_filter(nextp, ip) )
        {
            node_factoryp->free_traversal_node(nextp);
            return false;
        }
    }
    node_factoryp->free_traversal_node(nextp);

    return true;
}

    
bool Octree::traverse_with_filter(
    OctreeFilteredIterator *  ip )
{
    OctreeTraversalNode *  np = node_factoryp->new_traversal_node();
    np->center = center;
    np->length = length;
    np->depth = 0;
    np->nodep = rootp;
    this->init_traversal_node(np, np, -1);
    
    bool ret = this->traverse_node_with_filter(np, ip);
    node_factoryp->free_traversal_node(np);
    return ret;
}


double Octree::traverse_in_order(
    OctreeOrderedIterator *  ip )
{
    double ret = gmGOOGOL;
    bool have_result = false;

    // Add the root node.

    OctreeTraversalNode *  nodep = node_factoryp->new_traversal_node();
    nodep->center = center;
    nodep->length = length;
    nodep->nodep = rootp;
    nodep->depth = 0;
    this->init_traversal_node(nodep, nodep, -1);
    nodep->not_in_heap();
    nodep->heap_key(0);
    ip->nodepq.insert(nodep);

    // Iterate while there's nodes remaining.

    nodep = node_factoryp->new_traversal_node();
    while( ip->nodepq.size() != 0 )
    {
        // Free the last traversal node, and get the next one.

        node_factoryp->free_traversal_node(nodep);
        nodep = (OctreeTraversalNode *) ip->nodepq.extract();

        // If the weight of the next node is greater than the current
        // estimate, we're done.

        if ( ret < -nodep->heap_key() )
            break;

        // Call the iterator for the node.

        if ( ip->wants_internal_nodes || (nodep->nodep->childrenpp == NULL) )
        {
            double tmp = ip->ordered_iterate(nodep);
            if ( tmp < ret )
                ret = tmp;

            if ( nodep->nodep->childrenpp == NULL )
                continue;
        }

        // Add the children to the queue.
        // This can be optimized through coherence, but this is simple enough.

        double newlen = nodep->length * 0.5;
        int newdepth = nodep->depth + 1;
        for( int j=0; j<8; j++ )
        {
            if ( nodep->nodep->childrenpp[j] == NULL )
                continue; 

            OctreeTraversalNode *  childp = 
                node_factoryp->new_traversal_node();
            childp->center = nodep->center + 
                             newlen*octree_child_center_offset[j];
            childp->length = newlen;
            childp->nodep = nodep->nodep->childrenpp[j];
            childp->depth = newdepth;
            this->init_traversal_node(childp, nodep, j);

            childp->not_in_heap(); 
            double weight = ip->calculate_weight(childp);
            if ( weight != 0 )           
                childp->heap_key( -weight );
            else
                childp->heap_key( weight );

            ip->nodepq.insert(childp);
        }
    }

    // Free the remaining heap elements, the last seen element,
    // and reset the heap.

    node_factoryp->free_traversal_node(nodep);
    for(unsigned int j=0; j<ip->nodepq.size(); j++ )
        node_factoryp->free_traversal_node(
                          (OctreeTraversalNode *) ip->nodepq.item(j));
    ip->nodepq.reset();

    return ret;
}


// OctreeAnalyzer implementation --------------------------------------------

void OctreeAnalyzer::report_statistics_at_depth(
    int  d )
{
    double len = 2*octreep->length * pow(0.5, d);
		std::cout << "Depth " << d << " (len = " << len << "):\n"
         << "  " << internal_node_countp[d] << " internal / " 
         << leaf_countp[d] << " leaf\n";
}


void OctreeAnalyzer::report_statistics()
{
    int total_leaf = 0;
    int total_internal = 0;
    for( int i=0; i<octreep->max_actual_depth+1; i++ )
    {
        report_statistics_at_depth(i);
        std::cout << "\n";
        total_leaf += leaf_countp[i];
        total_internal += internal_node_countp[i];
    }

	std::cout << "Maximum depth: " << octreep->max_actual_depth << "\n"
         << "Total nodes: " << total_nodes << " ("
         << total_internal << " internal / "
         << total_leaf << " leaf)\n";
}

  
void OctreeAnalyzer::collect_statistics(
    Octree *  op )
{
    this->reset();
    octreep = op;
    leaf_countp = new int[octreep->max_actual_depth + 1];
    internal_node_countp = new int[octreep->max_actual_depth + 1];
    for( int i=0; i<octreep->max_actual_depth+1; i++ )
    {
        leaf_countp[i] = 0;
        internal_node_countp[i] = 0;
    }
    total_nodes = 0;

    octreep->traverse(this);
}
    

bool OctreeAnalyzer::iterate(
    OctreeTraversalNode *  np )
{
    total_nodes++;
    if ( np->nodep->childrenpp == NULL )
        leaf_countp[np->depth]++;
    else
        internal_node_countp[np->depth]++;
	return true;	// not correct, added by Wen, just to make it compile.
}


void OctreeAnalyzer::reset()
{
    if ( leaf_countp != NULL )
    {
        delete[] leaf_countp;
        leaf_countp = NULL;
    }
    if ( internal_node_countp != NULL )
    {
        delete[] internal_node_countp;
        internal_node_countp = NULL;
    }
}
