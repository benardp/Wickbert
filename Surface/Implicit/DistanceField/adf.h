#ifndef ADF_H
#define ADF_H

#include "mdsimhash.h"
#include "scalarlattice.h"
#include "triangle.h"

#include "renderable.h"

#include <vector>

class OctreeNode{
 public:
  bool is_leaf;
  int parent_index;

  OctreeNode * parent;
  Cell cell;

  Vector3d corner;
  scalar edge_length;

  void * data;

  scalar center_distance;
  union{
    scalar distances[8];
    OctreeNode * children[8];
  };

};

#ifdef WIN32
typedef hash_map<Cell, OctreeNode * ,cell_hash_compare> _ADF;
#else
typedef HASH_VERSION::hash_map<Cell, OctreeNode * , hash_cell, equal_cell > _ADF;
#endif 


class ADF : public _ADF, public Renderable {
 protected:
  scalar voxelsize;
  scalar inv_voxelsize;

 public:
  ADF(const scalar& vsize) { assert(vsize > 0); voxelsize = vsize; inv_voxelsize = 1.0 / vsize; }
  virtual ~ADF(){}

  scalar VoxelSize() const { return voxelsize; }

  static const Cell TraversalOrder[8];

  virtual void Render();
};

class ADFNode {
 public:
  std::vector<bool> types;
  union{
    scalar value;
    ADFNode * ptr;
  } children[8];

  ADFNode() : types(8) {}
};

class ADF2 {
 public:
  scalar edge_length;

  ADFNode root;

  // ADF2(const scalar edge_length) { assert(edge_length > 0); half_edge_length = edge_length / 2.0; }
 

};

#endif
