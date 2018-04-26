#include "sparsescalarlattice.h"

#include "opengl_utils.h"

#include <iostream>
#include <fstream>


SparseScalarLattice::SparseScalarLattice(scalar vsize) : ScalarLattice(vsize) {;}

bool SparseScalarLattice::HasValue(const Cell& c) const{
  return count(c) != 0;
}
  
bool SparseScalarLattice::ContainsPoint(const Vector3d& v) const{
  Vector3d vscaled = InverseVoxelSize() * v;
  
  Cell base((int)FLOOR(vscaled.X()),(int)FLOOR(vscaled.Y()), (int)FLOOR(vscaled.Z()));    
  
  for(int x = 0; x <= 1; x++){
    for(int y = 0; y <= 1; y++){
      for(int z = 0; z <= 1; z++){	
	if(!HasValue(base + Cell(x,y,z)))
	  return false;
      }
    }
  }
  return true;
}

AABoundingBox SparseScalarLattice::GetAABBox() const{
  AABoundingBox aabb;
  
  for(const_iterator i = begin(); i != end(); i++){
    const Cell c = (*i).first;
    aabb.include(VoxelSize() * c.toVector3d());
  }
  
  return aabb;
}


scalar SparseScalarLattice::Value(const Cell& c) const{
  SparseScalarLattice::const_iterator i = find(c);
  if(i == end()){
    throw InvalidCellException(c);
  } else {
    return i->second;
  }
}



void OutputToDF3(const SparseScalarLattice& lattice, std::string filename, scalar lower_bound, scalar upper_bound){
  std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary);

  if(!file){ std::cerr << "Unable to open file: " << filename << std::endl; return; }
  
  scalar min_value =  SCALAR_MAX;
  scalar max_value = -SCALAR_MAX;

  Cell bottom( INT_MAX, INT_MAX, INT_MAX);
  Cell top   (-INT_MAX,-INT_MAX,-INT_MAX);


  if(lattice.size() == 0){
    std::cerr << "Lattice contains no data.  Not writing file: " << filename << std::endl;
    return;
  }


  
  for(SparseScalarLattice::const_iterator i = lattice.begin(); i != lattice.end(); i++){
    const Cell& c       = (*i).first;
    const scalar& value = (*i).second;
    
    min_value = std::min(min_value,value);
    max_value = std::max(max_value,value);
    
    for(int n = 0; n < 3; n++){
      bottom[n] = std::min(bottom[n],c[n]);
      top[n] = std::max(top[n],c[n]);
    }
  }

  std::cout << "\n\nmax value " << max_value << " min value " << min_value << "\n";
  std::cout << "bottom " << lattice.VoxelSize()*bottom.toVector3d() << " top " << lattice.VoxelSize()*top.toVector3d() << "\n";

  if(!(lower_bound == SCALAR_MAX && upper_bound == SCALAR_MAX)){ //the user supplied some range
    assert(lower_bound < upper_bound);
    min_value = lower_bound;
    max_value = upper_bound;
  }


  int nx = top.X() - bottom.X() + 1;
  int ny = top.Y() - bottom.Y() + 1;
  int nz = top.Z() - bottom.Z() + 1;
  
  
  if(min_value >= max_value){ //lattice is everywhere constant
    max_value = min_value + 1;
  }
  
  //Assumes Little Endian!
  file.put(nx >> 8);
  file.put(nx & 0xff);
  file.put(ny >> 8);
  file.put(ny & 0xff);
  file.put(nz >> 8);
  file.put(nz & 0xff);

  for(int z = bottom.Z(); z <= top.Z(); z++){
    for(int y = bottom.Y(); y <= top.Y(); y++){
      for(int x = bottom.X(); x <= top.X(); x++){
	int intvalue = 0;
	Cell c(x,y,z);
	if(lattice.HasValue(c)){
	  scalar value = lattice.Value(c);
	  value = std::min(upper_bound,value); //Clamp to proper range
	  value = std::max(lower_bound,value);

	  intvalue += (int) (0xff * ((value - min_value)/(max_value - min_value)));;
	}
	//	file.put(value >> 8);
	file.put(intvalue & 0xff);
      }
    }
  }
  


  std::cout << "Lattice dimensions (" << nx << "," << ny << "," << nz << ")\n";

  file.close();
}


