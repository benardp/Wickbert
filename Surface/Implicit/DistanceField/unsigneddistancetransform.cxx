#include "distancetransform.h"
#include "stopwatch.h"


scalar TriangleSquaredDistance::SquaredDistanceTo(const Vector3d& p) const{
  return (p - NearestPointTo(p)).norm2();
}



UnsignedDistanceTransform::UnsignedDistanceTransform(const TriangleMesh2& mesh){

  const std::vector< Vector3d >& vertices = mesh.vertices;
  const std::vector< IndexTriple >& faces = mesh.faces;

  for(unsigned int i = 0; i < faces.size(); i++){
    const IndexTriple& corners = faces[i];
    triangles.push_back(Triangle(vertices[corners.A()],vertices[corners.B()],vertices[corners.C()]));
    tri_distances.push_back(TriangleSquaredDistance(triangles.back()));
  }
}


void UnsignedDistanceTransform::Voxelize(SparseScalarLattice& lattice, scalar max_dist){
  assert(max_dist > 0);

  Stopwatch sw;

  Vector3d diag = 1.0001 * Vector3d(max_dist,max_dist,max_dist);
  
  unsigned int total_calculations = 0;
  unsigned int total_successful = 0;

  const scalar voxelsize = lattice.VoxelSize();

  int last_print = 0;
  for(unsigned int i = 0; i < triangles.size(); i++){
    if((int) (100*((scalar) i / (scalar) triangles.size())) >= last_print + 5){
      last_print = (int) (100*((scalar) i / (scalar) triangles.size()));
      std::cout << "  " << last_print << "% finished" << std::endl;
    }

    const Triangle& tri = triangles[i];
    const TriangleSquaredDistance& TSD = tri_distances[i];

    AABoundingBox bbox = tri.BoundingBox();

    Vector3d v0 = bbox.v0 - diag;
    Vector3d v1 = bbox.v1 + diag;

    v0 /= voxelsize;
    v1 /= voxelsize;

    // only do the lattice points *contained* within the bbox
    Cell start(CEIL(v0.X()),CEIL(v0.Y()),CEIL(v0.Z()));
    Cell end(FLOOR(v1.X()),FLOOR(v1.Y()),FLOOR(v1.Z()));

    Vector3d current;

    for(int x = start.X(); x <= end.X(); x++){
      for(int y = start.Y(); y <= end.Y(); y++){

	scalar last_distance = SCALAR_MAX;

	for(int z = start.Z(); z <= end.Z(); z++){
	  Cell c(x,y,z);

	  current = Vector3d(voxelsize * x, voxelsize * y, voxelsize * z);

	  scalar distance = SQRT(TSD.SquaredDistanceTo(current));

	  total_calculations++;

	  //Don't try doing points father way than this (errors occur otherwise)
 	  if(distance > max_dist) { 
  	    if(distance >= last_distance){ break; }  //if distance does not decrease we're done with this line
 	    z += CEIL((distance - max_dist)/voxelsize) - 1;
  	    last_distance = distance;
 	    continue;
 	  }

	  last_distance = distance;

	  total_successful++;

	  SparseScalarLattice::iterator iter = lattice.find(c);
	  if(iter == lattice.end()){
	    lattice[c] = distance;
	  }
	  else if(distance < iter->second){
	    iter->second = distance;
	  } 
	}
      }
    }
  }

  scalar efficiency = 100*(scalar) total_successful /(scalar) total_calculations;
	std::cout << "Finished Voxelizing Mesh:  " << lattice.size() << " lattice points.  Efficiency: " << efficiency << "%   time:" << sw.time() << " seconds" << std::endl;
}

