#ifndef SPACEADF_INCLUDED // -*- C++ -*-
#define SPACEADF_INCLUDED

#ifdef _MSC_VER
#pragma once
#pragma warning(disable : 4786)
#endif


/************************************************************************
   Adaptively-sampled distance field.
   Third-generation:
    - re-usable octree-based implementation
    - hash table for samples
 ************************************************************************/
// controls if use the hashtable or the regular map
#define hashflag 1
// controls if to use the new hasher.
#define newhasher 1

// replace hash_map by map in STL #include <STLPort/hash_map>
#include "Surface/Implicit/Operator/UnaryOp.h"
#include "octree.h"
#include "libgm/gm.h"
#include "hashset.h"
#include "hashmap.h"
#include <map>
#include <set>
#include <limits>

/**
 * nan function for windows
 * this was changed after move the code from cygwin to VC++ since the STL is different
 */
double nan();

// Forward declarations.

class ADF;
class ADFNode;
class ADFTraversalNode;

/**
 * An ADF Sample Address.  Generally this is only used for manipulation
 * and maintenance of the address as an ADF is traversed.
 */
class ADFSampleAddress
{
public:
    void adjust_for_child(int depth, int child);
    void adjust_for_corner_sample(int depth, int corner);
    void adjust_for_sample(int depth, int sample);
    void adjust_coordinate(int coord, int delta, int depth);
   
    void set_to_root();

    double coord_to_point(ADF *adfp, int coord) const;
    void to_point(ADF *adfp, gmVector3 &ret ) const;

    const ADFSampleAddress& operator=(const ADFSampleAddress&  a)
    {
        address[0] = a.address[0];
        address[1] = a.address[1];
        address[2] = a.address[2];
        return a;
    }

    bool operator==(const ADFSampleAddress&  a) const
    {
        return ( (a.address[0] == address[0])
                 && (a.address[1] == address[1])
                 && (a.address[2] == address[2]) );
    }

public:
    /**
     * The x, y, and z components of the address.
     */
    unsigned short address[3];
};


/**
 * A hasher for an ADFSampleAddress.
 */
class ADFSampleHasher
{
public:
    int operator()(const ADFSampleAddress &addr) const;
};

/**
 * A compare operation for an ADFSampleAddress
 * used with STL set
 */
class ADFSampleSet
{
public:
  bool operator()(const ADFSampleAddress& a1, const ADFSampleAddress& a2) const;
};

/**
 * A compare operation for an ADFSampleAddress
 * used with STL set
 */
class ADFSampleMap
{
public:
  bool operator()(const ADFSampleAddress& a1, const ADFSampleAddress& a2) const;
};

/**
 * An ADF Sample.  If the value is NaN, then the sample is invalid.
 * Basically, we need a simple wrapper for a double, which allows an
 * invalid value, to be produced for addresses which aren't in the
 * ADF's sample hash map.
 */
class ADFSample
{
public:
  ADFSample(double  d) { distance = d; }
  ADFSample(const ADFSample&  s) { distance = s.distance; }
  ADFSample() { distance = nan(); }

  bool is_valid() const { return (distance != nan()); }

public:
    double distance;
};


/**
 * An ADF based on an octree.
 * The accuracy_isosurface denotes the isosurface below which the tolerance is guaranteed.  Allows only having accurate reconstruction near the surface.
 * In this new verion, accuracy_isosurface is not used because the ADF is storing a field, it needs to be accurate every where, including places where the values are not close to zero.
 */
class ADF : public Octree, public UnaryOp
{
protected:
    virtual void copy_traversal_node(OctreeTraversalNode *dp, OctreeTraversalNode *sp );
    virtual void init_traversal_node(OctreeTraversalNode *np, OctreeTraversalNode *pp, int child );

public:
  ADF() { m_f = NULL;}  /// default constructor needed by surface
//  ADF(DistanceFunction *df, double iso = gmGOOGOL);
  ADF(Implicit *m_f, Box<double> *box, double error);  /// constructor you normally use.
  ~ADF();  /// default destructor

    double get_sample(const ADFSampleAddress &addr, bool create = true);
    double get_distance(const gmVector3 &p, gmVector3 *normal = NULL, double calc_normal_threshold = gmGOOGOL);  /// old verion's proc

#if !defined(hashflag)
  std::map<ADFSampleAddress, ADFSample, ADFSampleMap> sample_map;  /// this new std::map replaced the old hashmap used in STLPort
#endif
#if defined(hashflag)
  hashmap<ADFSampleAddress, ADFSample, ADFSampleHasher, ADFSampleMap> sample_map;  /// this new std::map replaced the old hashmap used in STLPort
#endif

  virtual bool setChildren(std::vector<Implicit*> children);     /// Sets the operand of this operation.
  /* This should be removed. Everywhere the ADF accesses "implicit"
   * should be changed to an access of m_f. The member m_f is a pointer
   * to the implicit object the operator object actually operates on. 
   * The member m_f is already defined in UnaryOp, so it is inherited
   * and does not need to be defined here.
   * the implicit is called i_m, inhereited from UnaryOp
   * DistanceFunction *distance_function;
   *  Implicit *implicit; // used to create the ADF.
   */

#ifndef INTERVAL_EVAL_ONLY
    virtual double proc(gmVector3 x);    /// Evaluation of function.
#endif
    virtual Intervald proc(Box<double> x);

    virtual void procq(gmVector3, double*); /// Evaluation of dFdq.
    virtual void getq(double*);          /// Retreives parameters.
    virtual void _setq(double*);         /// Assigns parameters.
    virtual int qlen();         /// Returns the # of parameters.

  virtual void getqname(char** qn);      /// Retreives parameter names.
  MAKE_NAME();  /// factory call

  void createADF(Implicit *imp, Box<double> *b, double err=0.1);  /// calculate the ADF
  double length_table[16];  /// precomputered lengths for the data stucture
  double epsilon;  /// the error tolerance
  Box<double> *box;  /// an array of two gmVector3 that stores the two corner points
};


/**
 * Node types.  Mostly for denoting the reason a leaf was made a leaf.
 * Interior and exterior nodes are outside the accuracy iso-surface and
 * entirely interior or exterior to the surface.
 */
#define ADF_NODE_INTERNAL (0)
#define ADF_NODE_LEAF     (1)
#define ADF_NODE_INTERIOR (2)
#define ADF_NODE_EXTERIOR (3)


/**
 * An ADF node.
 */
class ADFNode : public OctreeNode
{
public:
    unsigned char node_type;
};


/**
 * Used for temporary storage of ADF nodes during traversal or construction.
 * Note that samples are not automatically extracted during traversal, for
 * efficiency; extract_samples_from must be called before performing any
 * real operations with the node.
 * get_local_coord_distance assumes the point is in the local coordinate
 * frame of the node.
 */
class ADFTraversalNode : public OctreeTraversalNode
{
public:
    void extract_samples_from(ADF *adfp);
    bool is_reconstruction_valid(ADF *adfp);
    double get_distance(const gmVector3 &p);
    double get_local_coord_distance(const gmVector3 &p);
    void estimate_normal(const gmVector3 &p, gmVector3 *normal);


public:
    ADFSampleAddress address;
    double samples[8];
};


/**
 * An iterator for use with mesh octrees, used to insert triangles into the octree.
 */
class ADFRefiner : public OctreeIterator
{
public:
    ADFRefiner() { wants_internal_nodes = false; }

    void refine(ADF *ap, double epsilon);

    virtual bool iterate(OctreeTraversalNode *np);

    
protected:
    virtual bool should_split(ADFTraversalNode *np);


public:
    ADF *adfp;
};


/**
 * Pruner for an ADF.
 * Prunes all unused samples.
 */
class ADFPruner : public OctreePruner
{
public:
    virtual void prune(Octree *op);
    virtual bool is_used(OctreeTraversalNode *np);


public:

#if !defined(hashflag)
  std::map<ADFSampleAddress, bool, ADFSampleMap> * used_samplesp;
#endif
#if defined(hashflag)
  hashmap<ADFSampleAddress, bool, ADFSampleHasher, ADFSampleMap> * used_samplesp;
#endif

};


/**
 * Statistics-gathering for an ADF.
 */
class ADFAnalyzer : public OctreeAnalyzer
{
public:
    ADFAnalyzer()
    {
        node_type_countsp[0] = NULL;
        node_type_countsp[1] = NULL;
        node_type_countsp[2] = NULL;
        node_type_countsp[3] = NULL;
    }

    ~ADFAnalyzer()
    {
        if ( node_type_countsp[0] != NULL )
        {
            delete[] node_type_countsp[0];
            delete[] node_type_countsp[1];
            delete[] node_type_countsp[2];
            delete[] node_type_countsp[3];
        }
    }


    virtual bool iterate(OctreeTraversalNode *np);
    virtual void report_statistics();
    virtual void report_statistics_at_depth(int d);
    virtual void reset();


public:
    int *  node_type_countsp[4];
};

#endif
