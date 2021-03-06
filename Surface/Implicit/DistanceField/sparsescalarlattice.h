#ifndef SPARSESCALARLATTICE_H
#define SPARSESCALARLATTICE_H

#include "mdsimhash.h"
#include "scalarlattice.h"
#include "triangle.h"

#ifdef WIN32
typedef std::map<Cell,scalar,cell_hash_compare> _SparseScalarLattice;
#else
typedef HASH_VERSION::hash_map<Cell, scalar, hash_cell, equal_cell > _SparseScalarLattice;
#endif 


class SparseScalarLattice : public _SparseScalarLattice, public ScalarLattice{
 public:
  SparseScalarLattice(scalar vsize = 1.0);

  bool HasValue(const Cell& c) const;
  bool ContainsPoint(const Vector3d& v) const;
  virtual AABoundingBox GetAABBox() const;
  virtual scalar Value(const Cell& c) const;  
};


void OutputToDF3(const SparseScalarLattice& lattice, std::string filename, scalar lower_bound = SCALAR_MAX, scalar upper_bound = SCALAR_MAX);

#endif
