/***********************************************************************
   Implementation of the ADF class.
 ***********************************************************************/

#include "adf.h"
#include "Surface/gmTNTconvert.h"

REGISTER_IMPLICIT(ADF,"UnaryOp:ADF");

// nan function for windows
double nan() {  return std::numeric_limits<double>::quiet_NaN();};
// ADF Implementation -------------------------------------------------------

// Protected Methods --------------------------------------------------------

void ADF::init_traversal_node(
    OctreeTraversalNode *  np,
    OctreeTraversalNode *  pp,
    int                    child )
{
    ADFTraversalNode *  atp = (ADFTraversalNode *) np;
    ADFTraversalNode *  app = (ADFTraversalNode *) pp;

    // Handle the root.

    if ( child == -1 )
    {
        atp->address.set_to_root();
        return ;
    }

    // Just copy the parent address and adjust it properly.

    atp->address = app->address;
    atp->address.adjust_for_child(app->depth, child);
}


void ADF::copy_traversal_node(
    OctreeTraversalNode *  dp,
    OctreeTraversalNode *  sp )
{
    Octree::copy_traversal_node(dp, sp);

    ADFTraversalNode *  adp = (ADFTraversalNode *) dp;
    ADFTraversalNode *  asp = (ADFTraversalNode *) sp;

    adp->address = asp->address;
    adp->samples[0] = asp->samples[0];
    adp->samples[1] = asp->samples[1];
    adp->samples[2] = asp->samples[2];
    adp->samples[3] = asp->samples[3];
    adp->samples[4] = asp->samples[4];
    adp->samples[5] = asp->samples[5];
    adp->samples[6] = asp->samples[6];
    adp->samples[7] = asp->samples[7];
}

/**
 * constructor using Implicit
 * \param imp a pointer to the implicit that is used to computer the values.
 * \param b a Box defining the region of space to evaluate proc() over.
 * \param error should be the tolerance that you want. ADF guarantees the approximated value is within this error tolerance.
 * \note This constructor just call createADF helper function, because serChildren also calls createADF.
 */
ADF::ADF(Implicit *imp, Box<double> *b, double error)
{
  createADF(imp, b, error);
}

ADF::~ADF()
{
    delete rootp;
    delete node_factoryp;
}

/**
 * create ADF from a implicit, box, error
 * \param imp a pointer to the implicit that is used to computer the values.
 * \param b a Box defining the region of space to evaluate proc() over.
 * \param err should be the tolerance that you want. ADF guarantees the approximated value is within this error tolerance.
 * \note This constructor accuately calls refine and prune to create the data stucture, so it may take a long time.
 */
void ADF::createADF(Implicit *imp, Box<double> *b, double err)
{
  m_f = imp;
  box=b;
  epsilon = gmGOOGOL;
  center=gmVector3(box->center()[0],box->center()[1],box->center()[2]);
  // this should use the maximum length, now we assumes it is a cube.
  // this length is not the width of the cell, but it is half of it.
  length=box->width()*0.5;

  // precompute length table at each depth level
  // length at level 0 is the root leve and it is twice the current length
  length_table[0] = length*2;
  for(int i=1; i<16; i++)
    length_table[i] = length_table[i-1]*0.5;

  // Create the top-level node and its samples.
  ADFNode * ap = new ADFNode();
  ap->parentp = NULL;
  ap->childrenpp = NULL;
  ap->node_type = ADF_NODE_LEAF;
  rootp = ap;

  ADFTraversalNode atn;
  atn.center = center;
  atn.length = length;
  atn.depth = 0;
  atn.nodep = rootp;
  atn.address.set_to_root();

  atn.extract_samples_from(this);

  // Create the node factory.
  node_factoryp = new OctreeTraversalNodePool<ADFTraversalNode>(32, 1);

  // create the data using refiner and pruner
  ADFRefiner refiner;
  refiner.refine(this, err);

  ADFPruner pruner;
  pruner.prune(this);

}

// Public methods -----------------------------------------------------------

/**
 * Overwrite the setChildren for UnaryOp. It then call createADF.
 * @param   children A vector of children (there should only be 1)
 * @returns False if there is not exactly 1 child in the vector.
 */

bool ADF::setChildren(std::vector<Implicit*> children)
{
  if (children.size() != 1)
    return false;

  // set implicit, box, epsilon
  // this is a default. need to be delete once implicit has a bounding box.
  Box<double> bbox(3);
  bbox[0].setInterval(-10,10);
  bbox[1].setInterval(-10,10);
  bbox[2].setInterval(-10,10);
  // default error is 1, makes the calculation faster.
  createADF(children[0],&bbox,1);

  return true;
}

/**
 * Retreives parameters.
 * @param q Parameter array.
 */
void ADF::getq(double* q) 
{
  q[0] = 0;

  // Check for child
  if (m_f!=NULL)
    m_f->getq(&q[1]);
} // end getq

/**
 * Assigns parameters.
 * @param q Parameter array.
 */
void ADF::_setq(double* q) {

  //m_r = q[0];

  // Check for child
  if (m_f!=NULL)
    m_f->_setq(&q[1]);

} // end _setq


/**
 * Evaluation of dFdq.
 * @param x    Point to evaluate at.
 * @param dfdq Array representing dfdq.
 */
void ADF::procq(gmVector3 x, double* dfdq) {

  if (m_f != NULL)//assert(m_f!=NULL);
  {
    dfdq[0] = -1;
    m_f->procq(x, &dfdq[1]);
  }

} // end procq


/**
 * Returns the # of parameters.
 */
int ADF::qlen() 
{
  int retval = 1;
  if (m_f!=NULL) 
    retval += m_f->qlen();
  return retval;
} // end qlen


/**
 * Retreives parameter names.
 * @param qn Paramater name array.
 */
void ADF::getqname(char** qn) {

  qn[0] = "ADF";

  UnaryOp::getqname(qn);

} // end getqname


/**
 * Create a new sample for the given address.
 * First, convert the sample address to a point in space,
 * and then extract the new sample from the distance function, and store it.
 * @note changed after using new hashset and hashtable.
 */

double ADF::get_sample(const ADFSampleAddress& addr, bool c)
{
  #if !defined(hashflag)
  std::map<ADFSampleAddress, ADFSample, ADFSampleMap>::iterator i = sample_map.find(addr);
  #endif
  #if defined(hashflag)
  hashmap<ADFSampleAddress, ADFSample, ADFSampleHasher, ADFSampleMap>::iterator i = sample_map.find(addr);
  #endif

    if ( i != sample_map.end() ) 
    return ((*i).second.distance);
    if ( !c )
        return nan();

    gmVector3 p;
    addr.to_point(this, p);
  double ret;
  if (m_f!=NULL)
    ret = m_f->proc(p);
  ADFSample s(ret);
    sample_map.insert(std::pair<const ADFSampleAddress, ADFSample>(addr, s));

    return ret;
}

/**
 * return the distance and the normal of a point in the ADF
 * @param p a point in ADF
 * @param normal a pointer to return for the normal
 * @param threshold default to 0
 * @return distance from the point to the surface
 */

double ADF::get_distance(const gmVector3& p, gmVector3 *normal, double threshold )
{
    ADFTraversalNode n, n2;
    this->get_leaf_for_point(p, &n, &n2);

    // If it's a leaf, proceed to do a distance calculation and estimate the normal.

    ADFNode *  ap = (ADFNode *) n.nodep;
    if ( ap->node_type == ADF_NODE_LEAF )
    {
        n.extract_samples_from(this);
        double d = n.get_distance(p);
        if ( (normal != NULL) && (d <= threshold) )
            n.estimate_normal(p, normal);
    
        return d;
    }

    // Outside the accurate region.  Just return HUGE.
    // TBD: we could put a better bound on it.

    return gmGOOGOL;
}

/**
 * ADF virtual functions inherented from Geometric, just calls get_distance
 * in the orignal ADF
 * @param x a point in ADF
 * @return evaluation of function
 */
#ifndef INTERVAL_EVAL_ONLY
double ADF::proc(gmVector3 x)
{
  if (m_f == NULL)
    return 0;
  else
    return get_distance(x);
}
#endif

Intervald ADF::proc(Box<double> x)
{
  Intervald retval;
  if (m_f != NULL)
    retval = Intervald(get_distance(convert(x.center())));
  return retval;
}

// ADFSampleAddress implementation ------------------------------------------

/**
 * Offsets of the samples, in terms of addresses.  Here, we're offsetting
 * from the center of the cube.  The order is such that it is most likely
 * to find a violating sample first (assumes more interior == more likely to
 * fail).
 */
const int ADF_sample_address_offset[][3] =
{
    // Middle of cube

    { 0,  0,  0},

    // Middle of faces

    { 0,  0, -1},
    { 0, -1,  0},
    {-1,  0,  0},
    { 0,  0,  1},
    { 0,  1,  0},
    { 1,  0,  0},

    // Middle of edges

    { 0, -1, -1},
    {-1,  0, -1},
    {-1, -1,  0},
    { 0, -1,  1},
    {-1,  0,  1},
    {-1,  1,  0},
    { 0,  1, -1},
    { 1,  0, -1},
    { 1, -1,  0},
    { 0,  1,  1},
    { 1,  0,  1},
    { 1,  1,  0},
};


/**
 * Sample offsets, used for trilinear reconstructions.  Here, in contrast
 * to the above, we use the -,-,- corner as 0,0,0, as is used for the
 * trilinear reconstruction.  The order must match ADF_sample_address_offset.
 */
const gmVector3 ADF_sample_offset[] =
{
    // Middle of cube

    gmVector3(0.5, 0.5, 0.5),

    // Middle of faces

    gmVector3(0.5, 0.5,   0),
    gmVector3(0.5,   0, 0.5),
    gmVector3(  0, 0.5, 0.5),
    gmVector3(0.5, 0.5,   1),
    gmVector3(0.5,   1, 0.5),
    gmVector3(  1, 0.5, 0.5),

    // Middle of edges

    gmVector3(0.5,   0,   0),
    gmVector3(  0, 0.5,   0),
    gmVector3(  0,   0, 0.5),
    gmVector3(0.5,   0,   1),
    gmVector3(  0, 0.5,   1),
    gmVector3(  0,   1, 0.5),
    gmVector3(0.5,   1,   0),
    gmVector3(  1, 0.5,   0),
    gmVector3(  1,   0, 0.5),
    gmVector3(0.5,   1,   1),
    gmVector3(  1, 0.5,   1),
    gmVector3(  1,   1, 0.5)
};


void ADFSampleAddress::adjust_coordinate(
    int  coord,
    int  delta,
    int  depth )
{
    // We start with the highest bit, and shift it right with added depth,
    // so that the depth may fluctuate (up to 15, of course), and so that the
    // address modification is simply an addition or subtraction.

    switch( delta )
    {
        case 1:
            address[coord] += (0x8000 >> depth);
            return ;
        case -1:
            address[coord] -= (0x8000 >> depth);
            return ;
        default:
            return ;
    }
}

/**
 * Root is at half the distance across the octree. So the root of the
 * octree is at 0x4000.
 *
 *   |--------|  0x8000
 *   |----|    0x4000
 *   |--|    0x2000
 *   |-|    0x1000 ...
 */
void ADFSampleAddress::set_to_root()
{
    // Note, the MSBit is for the entire length of the ADF; the root is
    // at the half-length.

    address[0] = 0x4000;
    address[1] = 0x4000;
    address[2] = 0x4000;
}

/**
 * Given the address for a cell, adjust the address so it
 * points to the given sample of the cell.
 * Depth is depth of current cell.
 * Sample is a number from 0 to 19 (same labeling as samples).
 * Sampling occurs at center of each cell, faces, edges (total 19).
 */
void ADFSampleAddress::adjust_for_sample(
    int  depth,
    int  sample )
{
    const int *  delta = ADF_sample_address_offset[sample];
    depth += 1;
    this->adjust_coordinate(0, delta[0], depth);
    this->adjust_coordinate(1, delta[1], depth);
    this->adjust_coordinate(2, delta[2], depth);
}

/**
 * Given the address for a cell, adjust the address so it
 * points to the given sample of the cell.
 * Depth is depth of current cell.
 * Sample is a number from 0 to 8 (same labeling as children).
 * Sampling occurs at (center of) vertices.
 */
void ADFSampleAddress::adjust_for_corner_sample(
    int  depth,
    int  corner )
{
    const int *  delta = octree_child_center_offset_ints[corner];
    depth += 1;
    this->adjust_coordinate(0, delta[0], depth);
    this->adjust_coordinate(1, delta[1], depth);
    this->adjust_coordinate(2, delta[2], depth);
}


/**
 * Samples in the center of one of the eight child cells.
 */
void ADFSampleAddress::adjust_for_child(
    int  depth,
    int  child )
{
    const int *  delta = octree_child_center_offset_ints[child];
    depth += 2;
    this->adjust_coordinate(0, delta[0], depth);
    this->adjust_coordinate(1, delta[1], depth);
    this->adjust_coordinate(2, delta[2], depth);
}

/**
 * Converts sample address to its spatial coordinats.
 * \param adfp Pointer to an ADF object
 * \param coord Coordinate to be converted (0=x, 1=y, 2=z)
 * Helper for to_point.
 */
double ADFSampleAddress::coord_to_point(
    ADF *  adfp,
    int    coord ) const
{
    double ret = adfp->center[coord] - adfp->length;
    for( unsigned short i=0x8000, j=0; i != 0; i >>= 1, j++ )
        if ( (i & address[coord]) != 0 )
            ret += adfp->length_table[j];

    return ret;
}

/**
 * Converts address to coordinates.
 */
void ADFSampleAddress::to_point(
    ADF *  adfp,
    gmVector3&  ret ) const
{
    ret[0] = this->coord_to_point(adfp, 0);
    ret[1] = this->coord_to_point(adfp, 1);
    ret[2] = this->coord_to_point(adfp, 2);
}


// ADFSampleHasher implementation -------------------------------------------
// wen this may not be used because we are using STL set
// This is based on code from a graphics gem for hashing 3D lattice addresses.

#define NBITS4 8
#define RBITS4 8
#define MASK4 0x0360

#define NBITS 8
#define RBITS 8
#define MASK 0x000F3C00
#define HASH(a,b,c) ((((a&MASK)<<NBITS|b&MASK)<<NBITS|c&MASK)>>RBITS)

int ADFSampleHasher::operator()(const ADFSampleAddress& addr) const
{
#if !defined (newhasher)
  return HASH(addr.address[0], addr.address[1], addr.address[2]);
//  this is a test to show that the address are not used efficiently, because last 8 bits are always 0
//  if (addr.address[0] & 0x00FF!=0)
//    i=i;
#endif

#if defined (newhasher)
// this is a better hasher
//  if (addr.address[0]!=0)
//    int i=0;
  int d=0x00FFFFFF;
  d=d&addr.address[0];
  d=d<<8;
  d=d|0x0000FFFF;
  d=d&(0x00FF0000|addr.address[1]);
  d=d|0x000000FF;
  int c=(0xFFFF0000|addr.address[2]);
  c=c >> 8;
  d=d&c;
  return d;
#endif
}

/**
 * ADFSampleSet implementation
 * this is used to replace ADFSampleHasher
 * hash_set need an operator equal
 */
bool ADFSampleSet::operator()(const ADFSampleAddress& a1, const ADFSampleAddress& a2) const
{
  // compare x, y, z component
  if (a1.address[0]==a2.address[0])
    if (a1.address[1]==a2.address[1])
      if (a1.address[2]==a2.address[2])
        return false;
      else
        return (a1.address[2]<a2.address[2]);
    else
      return (a1.address[1]<a2.address[1]);
  else
    return (a1.address[0]<a2.address[0]);
}

/**
 * ADFSampleMap implementation
 * this is used to replace ADFSampleHasher
 * hash_set need an operator less than
 */
bool ADFSampleMap::operator()(const ADFSampleAddress& a1, const ADFSampleAddress& a2) const
{
  // compare x, y, z component
  if (a1.address[0]==a2.address[0])
    if (a1.address[1]==a2.address[1])
      if (a1.address[2]==a2.address[2])
        return false;
      else
        return (a1.address[2]<a2.address[2]);
    else
      return (a1.address[1]<a2.address[1]);
  else
    return (a1.address[0]<a2.address[0]);
}


// ADFRefiner implementation ------------------------------------------------

/**
 * First, test that the node contains an iso-surface less than the accuracy_isosurface. 
 * If not, then there's no need to split. Of course, this test shouldn't be done on
 * the root node (and by not doing so, we're guaranteed to create the child samples
 * for the root node). 
 */
bool ADFRefiner::should_split(ADFTraversalNode * atp)
{
  // Validate the node according to the actual distances.  If
    // validation succeeds, mark the node as a leaf and return.

    ADFNode *  ap = (ADFNode *) atp->nodep;
    atp->extract_samples_from(adfp);

    if ( atp->is_reconstruction_valid(adfp) )
    {
        ap->node_type = ADF_NODE_LEAF;
        return false;
    }

    // Otherwise, we have to split.  Mark the node internal.

    ap->node_type = ADF_NODE_INTERNAL;
    return true;
}


/** Traverses to child nodes. Creates children if necessary.
 */
bool ADFRefiner::iterate(OctreeTraversalNode *np )
{
    ADFTraversalNode *  atp = (ADFTraversalNode *) np;
    ADFNode * ap = (ADFNode *) np->nodep;

    // Update the depth appropriately.

    if ( atp->depth > adfp->max_actual_depth )
        adfp->max_actual_depth = atp->depth;

    // Check to see if the node needs to be split.
    // should_split also sets the node_type appropriately for atp.

    if ( this->should_split(atp) )
    {
        // Create new nodes for the children.
 
        ap->childrenpp = (OctreeNode **) new ADFNode *[8];
        for( int i=0; i<8; i++ )
        {
            ap->childrenpp[i] = new ADFNode();
            ap->childrenpp[i]->parentp = ap;
            ((ADFNode *) ap->childrenpp[i])->node_type = ADF_NODE_LEAF;
        }
    }

    return true;
}
           
/** Checks tolerance and if not satisfied, continues to traverse.
 * \param ap Pointer to an ADF object
 * \param e Global tolerance
 */
void ADFRefiner::refine(ADF * ap, double e)
{
    // If the new tolerance is worse, do nothing.  Otherwise, record the
    // new tolerance, and traverse the octree.

    if ( ap->epsilon <= e )
        return ;

    ap->epsilon = e;
    adfp = ap;
    adfp->traverse(this);
}


// ADFTraversalNode implementation -------------------------------------

/**
 * Collect the weights for bilinear interpolation within the unit square.
 * The order of weights is quadrant order, of course.  x and y should be
 * transformed to between 0 and 1, as well.
 *
 * @param x The x coordinate of the point to be weighted.
 * @param y The y coordinate of the point to be weighted.
 * @param weights On return, the weights.  Must have room for 4 items.
 */
static void collect_bilinear_weights(
    double    x,
    double    y,
    double *  weights )
{
    double xp = 1 - x;
    double yp = 1 - y;

    weights[0] = x  * y ;
    weights[1] = xp * y ;
    weights[2] = x  * yp;
    weights[3] = xp * yp;
}


/**
 * Collect the weights for trilinear interpolation within a unit cube.
 * The order of weights corresponds to the sample order, of course.
 *
 * @param p The point to be weighted.  Should be in the unit cube.
 * @param weights On return, the weights.  Must have room for 8 items.
 */
static void collect_trilinear_weights(
    const gmVector3& p,
    double *    weights )
{
    gmVector3 q(1 - p[0], 1 - p[1], 1 - p[2]);
    weights[0] = q[0] * q[1] * q[2];
    weights[1] = p[0] * q[1] * q[2];
    weights[2] = q[0] * p[1] * q[2];
    weights[3] = p[0] * p[1] * q[2];
    weights[4] = q[0] * q[1] * p[2];
    weights[5] = p[0] * q[1] * p[2];
    weights[6] = q[0] * p[1] * p[2];
    weights[7] = p[0] * p[1] * p[2];
}

/** Grabs samples from ADF and stores locally.
 * TraversalNode has its own sample array.
 */
void ADFTraversalNode::extract_samples_from(ADF *  adfp )
{
    ADFSampleAddress orig = address;
    for( int i=0; i<8; i++ )
    {
        // We adjust at depth - 1, since we want to add the half-length
        // of this depth, not the child depths.

        address.adjust_for_corner_sample(this->depth, i);
        samples[i] = adfp->get_sample(address, true);
        address = orig;
    }
}

/** Checks ADF node reconstruction against samples.
 * \todo this is place to change if we want to change the sampling testing with the new interval testing.
 */
bool ADFTraversalNode::is_reconstruction_valid(ADF* adfp)
{
    ADFSampleAddress orig = address;
    for ( int i=0; i<19; i++ )
    {
        // Get the reconstructed distance.

        double d1 = this->get_local_coord_distance(ADF_sample_offset[i]);

        // Get the actual distance.  We need to save and restore this node's
        // address while doing so.

        address.adjust_for_sample(this->depth, i);
        double d2 = adfp->get_sample(address, true);
        address = orig;

        // If the difference is great, the reconstruction is invalid.

        if ( fabs(d1 - d2) > adfp->epsilon )
            return false;
    }

    return true;
}

/** Returns the reconstructed distance at world-coordinate point p.
 */
double ADFTraversalNode::get_distance(const gmVector3&  p )
{
    // Transform the input into the coordinate system of the cell.
    // This results in a vector of values from 0,0,0 in the nnn corner
    // to 1,1,1 in the ppp corner.

    gmVector3 pp = p - this->center;
    pp /= (this->length*2);
    pp[0] += 0.5;
    pp[1] += 0.5;
    pp[2] += 0.5;

    return this->get_local_coord_distance(pp);
}


/** Returns the reconstructed distance at cell-coordinate point p.
 */
double ADFTraversalNode::get_local_coord_distance(const gmVector3&  pp )
{
    // Trilinearly interpolate.

    double weights[8];
    collect_trilinear_weights(pp, weights);
    double ret = this->samples[0]*weights[0] +
                 this->samples[1]*weights[1] +
                 this->samples[2]*weights[2] +
                 this->samples[3]*weights[3] +
                 this->samples[4]*weights[4] +
                 this->samples[5]*weights[5] +
                 this->samples[6]*weights[6] +
                 this->samples[7]*weights[7];

    return ret;
}

/** Returns (in normal) the central-differenced gradient at
 * world-coordinate point p.
 */
void ADFTraversalNode::estimate_normal(
    const gmVector3&             p,
    gmVector3 *                  normal )
{
    // Transform into the cell coordinate system.

    gmVector3 pp = p - this->center;
    pp /= (this->length*2);
    pp[0] += 0.5;
    pp[1] += 0.5;
    pp[2] += 0.5;

    // Essentially use central differences across the cell.

    // X difference.

    double weights[4];
    collect_bilinear_weights(pp[1], pp[2], weights);
    (*normal)[0]  =   weights[0]*this->samples[0] 
                    + weights[1]*this->samples[2] 
                    + weights[2]*this->samples[4] 
                    + weights[3]*this->samples[6] 
                    - weights[0]*this->samples[1] 
                    - weights[1]*this->samples[3] 
                    - weights[2]*this->samples[5] 
                    - weights[3]*this->samples[7];

    // Y difference.

    collect_bilinear_weights(pp[0], pp[2], weights);
    (*normal)[1]  =   weights[0]*this->samples[0] 
                    + weights[1]*this->samples[1] 
                    + weights[2]*this->samples[4] 
                    + weights[3]*this->samples[5] 
                    - weights[0]*this->samples[2] 
                    - weights[1]*this->samples[3] 
                    - weights[2]*this->samples[6] 
                    - weights[3]*this->samples[7];
    
    // Z difference.

    collect_bilinear_weights(pp[0], pp[1], weights);
    (*normal)[1]  =   weights[0]*this->samples[0] 
                    + weights[1]*this->samples[1] 
                    + weights[2]*this->samples[2] 
                    + weights[3]*this->samples[3] 
                    - weights[0]*this->samples[4] 
                    - weights[1]*this->samples[5] 
                    - weights[2]*this->samples[6] 
                    - weights[3]*this->samples[7];
}


// ADFPruner implementation ------------------------------------------------

/** Pruner removes unused samples (and maybe cells).
 * Samples are unused if the cell isn't subdivided.
 * @note changed after new hashset and hashmap
 */
void ADFPruner::prune(Octree * op)
{
  #if !defined(hashflag)
  used_samplesp = new std::map<ADFSampleAddress, bool, ADFSampleMap>();
  #endif
  #if defined(hashflag)
  used_samplesp = new hashmap<ADFSampleAddress, bool, ADFSampleHasher, ADFSampleMap>();
  #endif

  OctreePruner::prune(op);
  // Delete all of the samples which weren't used.
  ADF *  adfp = (ADF *) op;

  /* old hasher
  */
  #if !defined(hashflag)
  std::map<ADFSampleAddress, ADFSample, ADFSampleMap>::iterator i;
  std::set<ADFSampleAddress, ADFSampleSet> to_remove;
  #endif
  #if defined(hashflag)
  hashmap<ADFSampleAddress, ADFSample, ADFSampleHasher, ADFSampleMap>::iterator i;
  hashset<ADFSampleAddress, ADFSampleHasher, ADFSampleSet> to_remove;
  #endif
  for( i = adfp->sample_map.begin(); i != adfp->sample_map.end(); ++i )
  {
    if ( used_samplesp->find(i->first) == used_samplesp->end() )
      to_remove.insert(i->first);
  }

  #if !defined(hashflag)
  std::set<ADFSampleAddress, ADFSampleSet>::iterator j;
  #endif
  #if defined(hashflag)
  hashset<ADFSampleAddress, ADFSampleHasher, ADFSampleSet>::iterator j;
  #endif
  for( j = to_remove.begin(); j != to_remove.end(); ++j )
    adfp->sample_map.erase(*j);
  // find duplicates and remove them.
  delete used_samplesp;
}

/** Indicates whether node np is used. Always returns true.
 * (Never want to prune cells, just samples.)
 * Side effect of marking the corner samples of leaf cells
 * as used.
 */
bool ADFPruner::is_used(OctreeTraversalNode * np)
{
    ADFNode *  ap = (ADFNode *) np->nodep;

    if ( ap->node_type == ADF_NODE_LEAF )
    {
        // Mark each of the corners of the node as used.

        ADFTraversalNode *  atp = (ADFTraversalNode *) np;
        ADFSampleAddress a;
        for( int i=0; i<8; i++ )
        {
            a = atp->address;
            a.adjust_for_corner_sample(atp->depth, i);
            (*used_samplesp)[a] = true;
        }
    }

    return true;
}


// ADFAnalyzer implementation ----------------------------------------------

bool ADFAnalyzer::iterate(
    OctreeTraversalNode *  np )
{
    OctreeAnalyzer::iterate(np);

    ADFNode *  ap = (ADFNode *) np->nodep;

    if ( node_type_countsp[0] == NULL )
    {
        node_type_countsp[0] = new int[octreep->max_actual_depth+1];
        node_type_countsp[1] = new int[octreep->max_actual_depth+1];
        node_type_countsp[2] = new int[octreep->max_actual_depth+1];
        node_type_countsp[3] = new int[octreep->max_actual_depth+1];

        for( int i=0; i<4; i++ )
            for( int j=0; j<octreep->max_actual_depth+1; j++ )
                node_type_countsp[i][j] = 0;
    }

    node_type_countsp[ap->node_type][np->depth]++ ;

    return true;
}


void ADFAnalyzer::report_statistics_at_depth(int d)
{
    OctreeAnalyzer::report_statistics_at_depth(d);
    std::cout << "  " << node_type_countsp[ADF_NODE_INTERIOR][d] << " interior / "
         << node_type_countsp[ADF_NODE_EXTERIOR][d] << " exterior / "
         << node_type_countsp[ADF_NODE_LEAF][d] << " true leaf\n";
}


void ADFAnalyzer::report_statistics()
{
    OctreeAnalyzer::report_statistics();

    int totals[4];
    totals[0] = totals[1] = totals[2] = totals[3] = 0;
    for( int i=0; i<octreep->max_actual_depth + 1; i++ )
    {
        totals[0] += node_type_countsp[0][i];
        totals[1] += node_type_countsp[1][i];
        totals[2] += node_type_countsp[2][i];
        totals[3] += node_type_countsp[3][i];
    }

    ADF *  adfp = (ADF *) octreep;
    std::cout << "  Total interior: " << totals[ADF_NODE_INTERIOR] << "\n"
         << "  Total exterior: " << totals[ADF_NODE_EXTERIOR] << "\n"
         << "  Total true leaves: " << totals[ADF_NODE_LEAF] << "\n"
         << "  Total samples: " << adfp->sample_map.size() << "\n";
}


void ADFAnalyzer::reset()
{
    OctreeAnalyzer::reset();

    if ( node_type_countsp[0] != NULL )
    {
        delete[] node_type_countsp[0];
        delete[] node_type_countsp[1];
        delete[] node_type_countsp[2];
        delete[] node_type_countsp[3];
        node_type_countsp[0] = NULL;
        node_type_countsp[1] = NULL;
        node_type_countsp[2] = NULL;
        node_type_countsp[3] = NULL;
    }
}
