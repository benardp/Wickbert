#include "cachingsparsescalarlattice.h"



scalar  CachingSparseScalarLattice::InterpValue(const  Vector3d& v){
  Vector3d vscaled = InverseVoxelSize() * v;

  Cell base((int)FLOOR(vscaled.X()),
	    (int)FLOOR(vscaled.Y()),
	    (int)FLOOR(vscaled.Z()));



  ScalarCache::const_iterator iter = s_cache.find(base);

  //If the values aren't already in the cache....
  if(iter == s_cache.end()){
    scalar * temp = s_cache[base].vals;
    for(int x = 0; x <= 1; x++){
      for(int y = 0; y <= 1; y++){
	for(int z = 0; z <= 1; z++){	
	  temp[4*x+2*y+z] = Value(base + Cell(x,y,z));
	}
      }
    }    
    iter = s_cache.find(base);
  } 

  const scalar * vals = iter->second.vals;
  
  scalar r = vscaled.X() - base.X();
  scalar s = vscaled.Y() - base.Y();
  scalar t = vscaled.Z() - base.Z();

  scalar result = 0;

  for(int x = 0; x <= 1; x++){
    scalar coef_r = ((x == 0) ? (1-r) : r);
    for(int y = 0; y <= 1; y++){
      scalar coef_s = ((y == 0) ? (1-s) : s);
      for(int z = 0; z <= 1; z++){	
	scalar coef_t = ((z == 0) ? (1-t) : t);
	scalar coef = coef_r * coef_s * coef_t;
	result += (scalar) coef * (scalar) vals[4*x+2*y+z];//Value(base + Cell(x,y,z));
      }
    }
  }

//   if(FABS(result - SparseScalarLattice::InterpValue(v)) >= 1e-5){
//     std::cout << "expected " << SparseScalarLattice::InterpValue(v) << " result " << result << "\n";
//     assert(FABS(result - SparseScalarLattice::InterpValue(v)) < 1e-5);
//   }
  
  return result;
}


Vector3d CachingSparseScalarLattice::InterpNormal(const Vector3d& v){
  Vector3d vscaled = InverseVoxelSize() * v;
  
  Cell base((int)FLOOR(vscaled.X()),
	    (int)FLOOR(vscaled.Y()),
	    (int)FLOOR(vscaled.Z()));

  VectorCache::const_iterator iter = n_cache.find(base);

  //If the values aren't already in the cache....
  if(iter == n_cache.end()){
    Vector3d * temp = n_cache[base].vals;
    for(int x = 0; x <= 1; x++){
      for(int y = 0; y <= 1; y++){
	for(int z = 0; z <= 1; z++){	
	  temp[4*x+2*y+z] = Normal(base + Cell(x,y,z));
	}
      }
    }    
    iter = n_cache.find(base);
  }


  const Vector3d * vals = iter->second.vals;
  
  
  scalar r = vscaled.X() - base.X();
  scalar s = vscaled.Y() - base.Y();
  scalar t = vscaled.Z() - base.Z();
  
  Vector3d result(0.0);
  
  for(int x = 0; x <= 1; x++){
    scalar coef_r = ((x == 0) ? (1-r) : r);
    for(int y = 0; y <= 1; y++){
      scalar coef_s = ((y == 0) ? (1-s) : s);
      for(int z = 0; z <= 1; z++){	
	scalar coef_t = ((z == 0) ? (1-t) : t);
	scalar coef = coef_r * coef_s * coef_t;
	result += coef * vals[4*x+2*y+z];//Value(base + Cell(x,y,z));
      }
    }
  }

  //  assert((result - SparseScalarLattice::InterpNormal(v)).norm() < 1e-5);

  return result;
}
 
scalar   CachingSparseScalarLattice::InterpMeanCurvature(const Vector3d& v) {

  return 0;
}
