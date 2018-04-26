#ifndef CACHINGSPARSESCALARLATTICE_H
#define CACHINGSPARSESCALARLATTICE_H


#include "sparsescalarlattice.h"




class CachingSparseScalarLattice : public SparseScalarLattice{
 protected:

  struct V8{ Vector3d vals[8]; };
  struct S8{ scalar   vals[8]; };

#ifdef WIN32
typedef HASH_VERSION::hash_map< Cell, V8, cell_hash_compare> VectorCache;
typedef HASH_VERSION::hash_map< Cell, S8, cell_hash_compare> ScalarCache;
#else
typedef HASH_VERSION::hash_map< Cell, V8, hash_cell, equal_cell > VectorCache;
typedef HASH_VERSION::hash_map< Cell, S8, hash_cell, equal_cell > ScalarCache;
#endif


  VectorCache n_cache;
  ScalarCache s_cache;

 public:
  CachingSparseScalarLattice(const SparseScalarLattice& ssl) : SparseScalarLattice(ssl) {}

  scalar   InterpValue(const  Vector3d& v);
  Vector3d InterpNormal(const Vector3d& v);
  scalar   InterpMeanCurvature(const Vector3d& v);

};

#endif
