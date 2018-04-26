#include "signeddistancetransform.h"

#include "stopwatch.h"
#include "mdsimhash.h"
#include "Surface/OpenMesh/SurfaceMesh.h"

//MEMORY LEAK it seems as if the Octree is not cleaned up... but it's only the destructor and not sooo bad...

bool ContainsAllCorners(const SparseScalarLattice& lattice, const Cell& cell);
void BuildAdaptiveOctree(OctreeNode * current, const scalar max_dist, const scalar abs_error, 
			 int max_subdivisions,
			 const std::vector<unsigned int>& nearby_triangles,
			 unsigned int& cell_count,
			 const std::vector< TriangleSignedDistance >& tri_distances);


// Methods used to Voxelize surface
class Index2Set{
 private:
  unsigned int pair[2];
 public:
  Index2Set() { pair[0] = 0; pair[1] = 0; }
  Index2Set(unsigned int a,unsigned int b) { pair[0] = std::min(a,b); pair[1] = std::max(a,b); }
  
  const unsigned int& A() const {return pair[0];}
  const unsigned int& B() const {return pair[1];}
};


#ifdef WIN32

struct index2set_hash_compare{
  enum
    {
      bucket_size = 4,
      min_buckets = 8
    };
  size_t operator()(const Index2Set& p) const{
    return p.A() ^ (p.B() << 16);
  }
  bool operator()(const Index2Set& p, const Index2Set& q) const{
	 if(p.A() < q.A())
		 return true;
	 if(p.A() > q.A())
		 return false;
	 if(p.B() < q.B())
		 return true;
	 return false;    
  }
};

typedef hash_map< Index2Set,Vector3d,index2set_hash_compare> EdgeMap;


#else
struct hash_index2set{
  size_t operator()(const Index2Set& p) const{
    return p.A() ^ (p.B() << 16);
  }
};

struct equal_index2set{
  bool operator()(const Index2Set& p, const Index2Set& q) const{
    return (p.A() == q.A()) && (p.B() == q.B());
  }
};
typedef HASH_VERSION::hash_map< Index2Set, Vector3d, hash_index2set, equal_index2set > EdgeMap;
#endif




SignedDistanceTransform::SignedDistanceTransform(const TriangleMesh2& mesh){

  const std::vector< Vector3d >& vertices = mesh.vertices;
  const std::vector< IndexTriple >& faces = mesh.faces;

  std::vector< Vector3d > weighted_facenormals(faces.size());

  //Compute face normals
  for(unsigned int i = 0; i < faces.size(); i++){
    const IndexTriple& tri = faces[i];
    Vector3d normal = cross(vertices[tri.B()] - vertices[tri.A()],vertices[tri.C()] - vertices[tri.A()]);
    assert(normal.norm() > 0.0);
    normalize(normal);
    weighted_facenormals[i] = normal;
  }

  
  //Compute angle-weighted edge normals
  EdgeMap edgemap;
  for(unsigned int i = 0; i < faces.size(); i++){
    const IndexTriple& tri = faces[i];
    
    Index2Set edges[3] = { Index2Set(tri.A(),tri.B()),
			   Index2Set(tri.B(),tri.C()),
			   Index2Set(tri.C(),tri.A()) };
    for(int j = 0; j < 3; j++){
      if(edgemap.count(edges[j]) == 0){
	edgemap[edges[j]] = weighted_facenormals[i];
      } else {
	edgemap[edges[j]] += weighted_facenormals[i];
      }
    }
  }

  std::vector< Vector3d > weighted_vertex_normals(vertices.size());
  //Compute angle-weighted vertex normals
  for(unsigned int i = 0; i < faces.size(); i++){
    const IndexTriple& tri = faces[i];
    
    const Vector3d& normal = weighted_facenormals[i];
    const Vector3d& v0 = vertices[tri.A()];
    const Vector3d& v1 = vertices[tri.B()];
    const Vector3d& v2 = vertices[tri.C()];

    Vector3d edge01 = (v1 - v0).normalized();
    Vector3d edge12 = (v2 - v1).normalized();
    Vector3d edge20 = (v0 - v2).normalized();

    //Must clamp ACOS parameter to [-1,1] (could lie outside this range due to inexact arithmetic)
    scalar angleA = ACOS(std::max((scalar)-1.0,std::min((scalar)1.0,-(edge01*edge20))));
    scalar angleB = ACOS(std::max((scalar)-1.0,std::min((scalar)1.0,-(edge01*edge12))));
    scalar angleC = ACOS(std::max((scalar)-1.0,std::min((scalar)1.0,-(edge12*edge20))));

    weighted_vertex_normals[tri.A()] += angleA * normal;
    weighted_vertex_normals[tri.B()] += angleB * normal;
    weighted_vertex_normals[tri.C()] += angleC * normal;
  }

 
  for(unsigned int i = 0; i < faces.size(); i++){
    const IndexTriple& corners = faces[i];

    triangles.push_back(Triangle(vertices[corners.A()],vertices[corners.B()],vertices[corners.C()]));

    TriangleSignedDistance TSD(triangles.back());

    TSD.angle_weighted_normals[TriangleNearestPoint::CORNER_A_REGION] = weighted_vertex_normals[corners.A()];
    TSD.angle_weighted_normals[TriangleNearestPoint::CORNER_B_REGION] = weighted_vertex_normals[corners.B()];
    TSD.angle_weighted_normals[TriangleNearestPoint::CORNER_C_REGION] = weighted_vertex_normals[corners.C()];
    TSD.angle_weighted_normals[TriangleNearestPoint::EDGE_AB_REGION] = edgemap[Index2Set(corners.A(),corners.B())];
    TSD.angle_weighted_normals[TriangleNearestPoint::EDGE_AC_REGION] = edgemap[Index2Set(corners.A(),corners.C())];
    TSD.angle_weighted_normals[TriangleNearestPoint::EDGE_BC_REGION] = edgemap[Index2Set(corners.B(),corners.C())];
    TSD.angle_weighted_normals[TriangleNearestPoint::FACE_REGION] = weighted_facenormals[i];

    tri_distances.push_back(TSD);
  }
}

SignedDistanceTransform::SignedDistanceTransform(const SurfaceMesh& mesh){
	typedef std::pair<SurfaceMesh::VertexHandle,SurfaceMesh::VertexHandle> VertexHandlePair;
	std::map<VertexHandlePair,Vector3d> edgeNormals;
	
	for (SurfaceMesh::ConstEdgeIter e_it = mesh.edges_begin(); e_it!=mesh.edges_end(); ++e_it) 
	{
		SurfaceMesh::HalfedgeHandle h1 (mesh.halfedge_handle(e_it.handle(),0));
		SurfaceMesh::HalfedgeHandle h2 (mesh.halfedge_handle(e_it.handle(),1));

		SurfaceMesh::Point n1 (mesh.normal(mesh.face_handle(h1)));
		SurfaceMesh::Point n2 (mesh.normal(mesh.face_handle(h1)));
		SurfaceMesh::Point edgeN((n1+n2).normalize());
		Vector3d edgeNormal(edgeN[0],edgeN[1],edgeN[2]);
		edgeNormals[VertexHandlePair(mesh.from_vertex_handle(h1), mesh.to_vertex_handle(h1))]=edgeNormal;
		edgeNormals[VertexHandlePair(mesh.to_vertex_handle(h1), mesh.from_vertex_handle(h1))]=edgeNormal;
	}

	for (SurfaceMesh::ConstFaceIter f_it = mesh.faces_begin(); f_it!=mesh.faces_end(); ++f_it) 
	{
		std::vector<gmVector3> vertices;
		std::vector<SurfaceMesh::VertexHandle> vertexHandles;

		for (SurfaceMesh::ConstFaceVertexIter fv_it = mesh.cfv_iter(f_it.handle()); fv_it; ++fv_it) 
		{
			SurfaceMesh::Point p=mesh.point(fv_it.handle());
			vertexHandles.push_back(fv_it.handle());
			gmVector3 pV(p[0],p[1],p[2]);
			vertices.push_back(pV);
		}
		assert(vertices.size()==3);
		
		TriangleSignedDistance TSD(vertices[0],vertices[1],vertices[2]);
		SurfaceMesh::Point n0(mesh.normal(vertexHandles[0]));
		TSD.angle_weighted_normals[TriangleNearestPoint::CORNER_A_REGION] = Vector3d(n0[0],n0[1],n0[2]);
		SurfaceMesh::Point n1(mesh.normal(vertexHandles[1]));
		TSD.angle_weighted_normals[TriangleNearestPoint::CORNER_B_REGION] = Vector3d(n1[0],n1[1],n1[2]);
		SurfaceMesh::Point n2(mesh.normal(vertexHandles[2]));		
		TSD.angle_weighted_normals[TriangleNearestPoint::CORNER_C_REGION] = Vector3d(n2[0],n2[1],n2[2]);

		TSD.angle_weighted_normals[TriangleNearestPoint::EDGE_AB_REGION] = edgeNormals[VertexHandlePair(vertexHandles[0],vertexHandles[1])];
		TSD.angle_weighted_normals[TriangleNearestPoint::EDGE_AC_REGION] = edgeNormals[VertexHandlePair(vertexHandles[0],vertexHandles[2])];
		TSD.angle_weighted_normals[TriangleNearestPoint::EDGE_BC_REGION] = edgeNormals[VertexHandlePair(vertexHandles[1],vertexHandles[2])];
		SurfaceMesh::Point fn=mesh.normal(f_it.handle());
		TSD.angle_weighted_normals[TriangleNearestPoint::FACE_REGION] = Vector3d( fn[0], fn[1], fn[2]);

		tri_distances.push_back(TSD);
	}
}

double SignedDistanceTransform::signedDistanceTo(const gmVector3 & pos) const
{
	double dist=DBL_MAX;
	
	Vector3d p(pos[0],pos[1],pos[2]);
	for (std::vector< TriangleSignedDistance >::const_iterator iter=tri_distances.begin();iter!=tri_distances.end();++iter)
	{
		double currDist=iter->SignedDistanceTo(p);
		if (fabs(dist)>fabs(currDist))
			dist=currDist;
	}
	return dist;
}



void SignedDistanceTransform::Voxelize(SparseScalarLattice& lattice, scalar max_dist){
  assert(max_dist > 0);

  Stopwatch sw;


  const scalar voxelsize = lattice.VoxelSize();

  Vector3d diag = 1.0001 * Vector3d(max_dist,max_dist,max_dist);
  
  unsigned int total_calculations = 0;
  unsigned int total_successful = 0;

  int last_print = 0;
  for(unsigned int i = 0; i < triangles.size(); i++){
    if((int) (100*((scalar) i / (scalar) triangles.size())) >= last_print + 5){
      last_print = (int) (100*((scalar) i / (scalar) triangles.size()));
      std::cout << "  " << last_print << "% finished" << std::endl;
    }

    const Triangle& tri = triangles[i];
    const TriangleSignedDistance& TSD = tri_distances[i];

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

	  scalar signed_distance = TSD.SignedDistanceTo(current);
	  scalar distance = FABS(signed_distance);

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
	    lattice[c] = signed_distance;
	  }
	  else if(distance < FABS(iter->second)){
	    iter->second = signed_distance;
	  } 
	}
      }
    }
  }

  scalar efficiency = 100*(scalar) total_successful /(scalar) total_calculations;
	std::cout << "Finished Voxelizing Mesh:  " << lattice.size() << " lattice points.  Efficiency: " << efficiency << "%   time:" << sw.time() << " seconds" << std::endl;
}


#ifdef WIN32
typedef hash_map<Cell, std::vector<unsigned int>, cell_hash_compare> IndexLattice;
#else
typedef HASH_VERSION::hash_map<Cell, std::vector<unsigned int>, hash_cell, equal_cell > IndexLattice;
#endif 

void SignedDistanceTransform::Voxelize(ADF& adf, scalar max_dist, int max_subdivisions, scalar abs_error){
  SparseScalarLattice lattice(adf.VoxelSize());
  Voxelize(adf, lattice, max_dist, max_subdivisions, abs_error);
}

void SignedDistanceTransform::Voxelize(ADF& adf, SparseScalarLattice& lattice, scalar max_dist, int max_subdivisions, scalar abs_error){

  assert(max_dist > 0);
  assert(adf.VoxelSize() == lattice.VoxelSize());

  Stopwatch sw;

  IndexLattice index_lattice;

  //Must include all corners of every cell that has points max_dist away from surface
  max_dist += SQRT((scalar)3)*adf.VoxelSize();

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
    const TriangleSignedDistance& TSD = tri_distances[i];

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

	  scalar signed_distance = TSD.SignedDistanceTo(current);
	  scalar distance = FABS(signed_distance);

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

	  index_lattice[c].push_back(i);

	  SparseScalarLattice::iterator iter = lattice.find(c);
	  if(iter == lattice.end()){
	    lattice[c] = signed_distance;
	  }
	  else if(distance < FABS(iter->second)){
	    iter->second = signed_distance;
	  } 
	}
      }
    }
  }

  scalar efficiency = 100*(scalar) total_successful /(scalar) total_calculations;
  std::cout << "Finished Voxelizing Mesh:  " << lattice.size() << " lattice points.  Efficiency: " << efficiency << "%   time:" << sw.time() << " seconds" << std::endl;


  std::vector<unsigned int> indices_last_referenced(tri_distances.size(), 0);

  unsigned int cell_count = 0;

  for(SparseScalarLattice::const_iterator i = lattice.begin(); i != lattice.end(); i++){
    const Cell& corner = i->first;

    if(ContainsAllCorners(lattice,corner)){  
      cell_count++;

      //Make a new ADF Node for this cell
      OctreeNode * node = new OctreeNode();
      
      //Fill in the node's data
      node->corner = voxelsize * corner.toVector3d();
      node->edge_length = voxelsize;
      node->is_leaf = true;

      for(unsigned int n = 0; n < 8; n++){
	node->distances[n] = lattice[corner + ADF::TraversalOrder[n]];
      }
      
      //Make a record of which triangles were near this cell
      std::vector<unsigned int> nearby_triangle_indices;
      for(unsigned int n = 0; n < 8; n++){
	const std::vector<unsigned int>& neighbors = index_lattice[corner + ADF::TraversalOrder[n]];
	for(unsigned int m = 0; m < neighbors.size(); m++){
	  const unsigned int current_index = neighbors[m];
	  if(indices_last_referenced[current_index] != cell_count){
	    indices_last_referenced[current_index] = cell_count;
	    nearby_triangle_indices.push_back(current_index);
	  }
	}
      }
      
      BuildAdaptiveOctree(node, max_dist, abs_error, max_subdivisions, nearby_triangle_indices, cell_count, tri_distances);

      //Insert this octree into the adf
      adf[corner] = node;
    }
  }

  std::cout << "Total Octree Cells " << cell_count << std::endl;
}

void LinearInterpolate19(scalar known[], scalar interpolated_values[]){
  //Edges  
  interpolated_values[0]  = 0.50*(known[0] + known[1]);
  interpolated_values[1]  = 0.50*(known[0] + known[2]);
  interpolated_values[3]  = 0.50*(known[1] + known[3]);
  interpolated_values[4]  = 0.50*(known[2] + known[3]);

  interpolated_values[5]  = 0.50*(known[0] + known[4]);
  interpolated_values[7]  = 0.50*(known[1] + known[5]);
  interpolated_values[11] = 0.50*(known[2] + known[6]);
  interpolated_values[13] = 0.50*(known[3] + known[7]);

  interpolated_values[14] = 0.50*(known[4] + known[5]);
  interpolated_values[15] = 0.50*(known[4] + known[6]);
  interpolated_values[17] = 0.50*(known[5] + known[7]);
  interpolated_values[18] = 0.50*(known[6] + known[7]);

  //Face Centers
  interpolated_values[2] = 0.50*(interpolated_values[0] + interpolated_values[4]);
  
  interpolated_values[6]  = 0.50*(interpolated_values[0] + interpolated_values[14]);
  interpolated_values[8]  = 0.50*(interpolated_values[1] + interpolated_values[15]);
  interpolated_values[10] = 0.50*(interpolated_values[3] + interpolated_values[17]);
  interpolated_values[12] = 0.50*(interpolated_values[4] + interpolated_values[18]);

  interpolated_values[16] = 0.50*(interpolated_values[14] + interpolated_values[18]);

  //Middle Point
  interpolated_values[9] = 0.50*(interpolated_values[2] + interpolated_values[16]);
}

void LinearInterpolate19(const Vector3d& corner, const scalar edge_length, Vector3d samples[]){
  scalar half_edge_length = edge_length / (scalar) 2.0;

  for(unsigned int i = 0; i < 19; i++){
    samples[i] = corner;
  }

  //Increment Zs
  samples[ 0][2] += half_edge_length;
  samples[ 6][2] += half_edge_length;
  samples[14][2] += half_edge_length;
  samples[ 2][2] += half_edge_length;
  samples[ 9][2] += half_edge_length;
  samples[16][2] += half_edge_length;
  samples[ 4][2] += half_edge_length;
  samples[12][2] += half_edge_length;
  samples[18][2] += half_edge_length;


  samples[ 3][2] += edge_length;
  samples[ 7][2] += edge_length;
  samples[10][2] += edge_length;
  samples[13][2] += edge_length;
  samples[17][2] += edge_length;
  

  //Increment Ys
  samples[ 1][1] += half_edge_length;
  samples[ 2][1] += half_edge_length;
  samples[ 3][1] += half_edge_length;
  samples[ 8][1] += half_edge_length;
  samples[ 9][1] += half_edge_length;
  samples[10][1] += half_edge_length;
  samples[15][1] += half_edge_length;
  samples[16][1] += half_edge_length;
  samples[17][1] += half_edge_length;


  samples[ 4][1] += edge_length;
  samples[11][1] += edge_length;
  samples[12][1] += edge_length;
  samples[13][1] += edge_length;
  samples[18][1] += edge_length;


  //Increment Xs  
  samples[ 5][0] += half_edge_length;
  samples[ 6][0] += half_edge_length;
  samples[ 7][0] += half_edge_length;
  samples[ 8][0] += half_edge_length;
  samples[ 9][0] += half_edge_length;
  samples[10][0] += half_edge_length;
  samples[11][0] += half_edge_length;
  samples[12][0] += half_edge_length;
  samples[13][0] += half_edge_length;


  samples[14][0] += edge_length;
  samples[15][0] += edge_length;
  samples[16][0] += edge_length;
  samples[17][0] += edge_length;
  samples[18][0] += edge_length;
}


void TrueDistances19(const std::vector<unsigned int>& nearby_triangles,
		     const std::vector< TriangleSignedDistance >& tri_distances,
		     Vector3d samples[],
		     scalar true_distances[]){

  for(unsigned int i = 0; i < 19; i++){
    true_distances[i] = SCALAR_MAX;
  }


  for(unsigned int i = 0; i < nearby_triangles.size(); i++){
    const TriangleSignedDistance& TSD = tri_distances[nearby_triangles[i]];
    
    for(unsigned int j = 0; j < 19; j++){
      scalar signed_distance = TSD.SignedDistanceTo(samples[j]);
      if(FABS(signed_distance) < FABS(true_distances[j])){
	true_distances[j] = signed_distance;
      }
    }
  }
}
		     

void InitializeChildren(OctreeNode * current, const scalar true_distances[]){
  //Change the node type
  current->is_leaf = false;
  
  //Store the current distances
  scalar corner_distances[8];
  for(unsigned int i = 0; i < 8; i++){
    corner_distances[i] = current->distances[i];
  }

  //Add and initialize children
  const Vector3d& corner = current->corner;
  scalar half_edge_length = current->edge_length / (scalar) 2.0;
  for(unsigned int i = 0; i < 8; i++){
    current->children[i] = new OctreeNode();
    
    current->children[i]->is_leaf     = true;
    current->children[i]->corner      = corner;
    current->children[i]->edge_length = half_edge_length;
  }
  

  OctreeNode * node;

  node = current->children[0];
  //  node->corner = corner;
  node->distances[0] = corner_distances[0];
  node->distances[1] = true_distances[0];
  node->distances[2] = true_distances[1];
  node->distances[3] = true_distances[2];
  node->distances[4] = true_distances[5];
  node->distances[5] = true_distances[6];
  node->distances[6] = true_distances[8];
  node->distances[7] = true_distances[9];

  node = current->children[1];
  node->corner[2] += half_edge_length;
  node->distances[0] = true_distances[0];
  node->distances[1] = corner_distances[1];
  node->distances[2] = true_distances[2];
  node->distances[3] = true_distances[3];
  node->distances[4] = true_distances[6];
  node->distances[5] = true_distances[7];
  node->distances[6] = true_distances[9];
  node->distances[7] = true_distances[10];

  node = current->children[2];
  node->corner[1] += half_edge_length;
  node->distances[0] = true_distances[1];
  node->distances[1] = true_distances[2];
  node->distances[2] = corner_distances[2];
  node->distances[3] = true_distances[4];
  node->distances[4] = true_distances[8];
  node->distances[5] = true_distances[9];
  node->distances[6] = true_distances[11];
  node->distances[7] = true_distances[12];

  node = current->children[3];
  node->corner[2] += half_edge_length;
  node->corner[1] += half_edge_length;
  node->distances[0] = true_distances[2];
  node->distances[1] = true_distances[3];
  node->distances[2] = true_distances[4];
  node->distances[3] = corner_distances[3];
  node->distances[4] = true_distances[9];
  node->distances[5] = true_distances[10];
  node->distances[6] = true_distances[12];
  node->distances[7] = true_distances[13];

  node = current->children[4];
  node->corner[0] += half_edge_length;
  node->distances[0] = true_distances[5];
  node->distances[1] = true_distances[6];
  node->distances[2] = true_distances[8];
  node->distances[3] = true_distances[9];
  node->distances[4] = corner_distances[4];
  node->distances[5] = true_distances[14];
  node->distances[6] = true_distances[15];
  node->distances[7] = true_distances[16];

  node = current->children[5];
  node->corner[2] += half_edge_length;
  node->corner[0] += half_edge_length;
  node->distances[0] = true_distances[6];
  node->distances[1] = true_distances[7];
  node->distances[2] = true_distances[9];
  node->distances[3] = true_distances[10];
  node->distances[4] = true_distances[14];
  node->distances[5] = corner_distances[5];
  node->distances[6] = true_distances[16];
  node->distances[7] = true_distances[17];

  node = current->children[6];
  node->corner[1] += half_edge_length;
  node->corner[0] += half_edge_length;
  node->distances[0] = true_distances[8];
  node->distances[1] = true_distances[9];
  node->distances[2] = true_distances[11];
  node->distances[3] = true_distances[12];
  node->distances[4] = true_distances[15];
  node->distances[5] = true_distances[16];
  node->distances[6] = corner_distances[6];
  node->distances[7] = true_distances[18];

  node = current->children[7];
  node->corner[2] += half_edge_length;
  node->corner[1] += half_edge_length;
  node->corner[0] += half_edge_length;
  node->distances[0] = true_distances[9];
  node->distances[1] = true_distances[10];
  node->distances[2] = true_distances[12];
  node->distances[3] = true_distances[13];
  node->distances[4] = true_distances[16];
  node->distances[5] = true_distances[17];
  node->distances[6] = true_distances[18];
  node->distances[7] = corner_distances[7];
}


void FilterTriangles(const std::vector<unsigned int>& nearby_triangles,
		     const std::vector< TriangleSignedDistance >& tri_distances,
		     std::vector<unsigned int>& filtered_triangles,
		     const Vector3d& center_sample,
		     const scalar edge_length)
{
  std::vector< scalar > distances(nearby_triangles.size());
  
  
  scalar min_distance = SCALAR_MAX;
  for(unsigned int i = 0; i < nearby_triangles.size(); i++){
    const TriangleSignedDistance& TSD = tri_distances[nearby_triangles[i]];  
    const scalar distance = FABS(TSD.SignedDistanceTo(center_sample));

    distances[i] = distance;
    min_distance = std::min(min_distance, distances[i]);			
  }

  //Triangles farther away than this can't possibly be near this octree cell
  const scalar upper_bound = (SQRT(3.0) / 2.0) * edge_length + min_distance;
  
  for(unsigned int i = 0; i < distances.size(); i++){
    if(distances[i] <= upper_bound){
      filtered_triangles.push_back(nearby_triangles[i]);
    }
  }
  
  //  std::cout << "was " << nearby_triangles.size() << " now " << filtered_triangles.size() << "\n"; 

}
		     
void BuildAdaptiveOctree(OctreeNode * current, 
			 const scalar max_dist, const scalar abs_error, 
			 int max_subdivisions,
			 const std::vector<unsigned int>& nearby_triangles,
			 unsigned int& cell_count,
			 const std::vector< TriangleSignedDistance >& tri_distances){

  if(max_subdivisions == 0){
    return;
  }

  cell_count++;


  //First compute values at interpolated samples
  scalar interpolated_values[19];
  LinearInterpolate19(current->distances, interpolated_values);


  //Compute the actual sample points needed
  Vector3d samples[19];
  LinearInterpolate19(current->corner, current->edge_length, samples);

  //Filter out distant triangles
  std::vector<unsigned int> filtered_triangles;
  FilterTriangles(nearby_triangles, tri_distances, filtered_triangles, samples[9], current->edge_length);

  //Now compute the true distances to the same points
  scalar true_distances[19];
  TrueDistances19(filtered_triangles, tri_distances, samples, true_distances);


  //See if we need to subdivide
  bool within_tolerance = true;
  for(unsigned int i = 0; i < 19; i++){
    if(FABS(interpolated_values[i] - true_distances[i]) > abs_error){
      within_tolerance = false;
      break;
    }
  }
  if(within_tolerance){ return;  } 
  
  
  // Subdividing //
  
  InitializeChildren(current, true_distances);

  
  for(unsigned int i = 0; i < 8; i++){
    BuildAdaptiveOctree(current->children[i], 
			max_dist, abs_error, 
			max_subdivisions - 1,
			filtered_triangles,
			cell_count,
			tri_distances);
  }

}



bool ContainsAllCorners(const SparseScalarLattice& lattice, const Cell& cell){
  for(int x = 0; x <= 1; x++){
    for(int y = 0; y <= 1; y++){
      for(int z = 0; z <= 1; z++){
	if(!lattice.HasValue(cell + Cell(x,y,z))){
	  return false;
	}
      }
    }
  }

  return true;
}

scalar TriangleSignedDistance::SignedDistanceTo(const Vector3d& p) const{
  VoronoiRegion region;
  
  Vector3d surface_to_p = p - NearestPointToVerbose(p, region);
  
  scalar distance = surface_to_p.norm();

  if(angle_weighted_normals[region] * surface_to_p < 0){
    return -distance;
  } else {
    return distance;
  }
}





void SignedDistanceTransform::Test(){
  //Make sure we allways have the same sequence of random numbers
  srand(0);

  //The 19 intepolated samples in a cube with edge length = 2.0
  Vector3d correct_samples[19] = {
    Vector3d(0,0,1),
    Vector3d(0,1,0),
    Vector3d(0,1,1),
    Vector3d(0,1,2),
    Vector3d(0,2,1),      
    
    Vector3d(1,0,0),
    Vector3d(1,0,1),
    Vector3d(1,0,2),
    Vector3d(1,1,0),
    Vector3d(1,1,1),
    Vector3d(1,1,2),
    Vector3d(1,2,0),
    Vector3d(1,2,1),
    Vector3d(1,2,2),
    
    Vector3d(2,0,1),
    Vector3d(2,1,0),
    Vector3d(2,1,1),
    Vector3d(2,1,2),
    Vector3d(2,2,1)
  };



  //LinearInterpolate19(Vector3d)
  {
    std::cerr << "SignedDistanceTransform::LinearInterpolate19(vector):";
    bool test_passed = true;

    Vector3d samples[19];
    LinearInterpolate19(Vector3d(0,0,0), 2, samples);
    
    for(int i = 0; i < 19; i++){
      if(samples[i] != correct_samples[i]){
	std::cerr << " Error: sample[" << i << "] = " << samples[i] << " != " << correct_samples[i] << std::endl;
      }
    }
    if(test_passed){
      std::cerr << " passed\n";
    }
  }

  //LinearInterpolate19(scalar)
  {
    std::cerr << "SignedDistanceTransform::LinearInterpolate19(scalar):";
    bool test_passed = true;

    scalar known[8];
    for(int i = 0; i < 8; i++){
      known[i] = rand_scalar(1.0);
    }

    scalar correct_values[19];
    for(int i = 0; i < 19; i++){
      scalar r = correct_samples[i].X() / 2.0;
      scalar s = correct_samples[i].Y() / 2.0;
      scalar t = correct_samples[i].Z() / 2.0;

      correct_values[i] 
	= (1 - r) * (1 - s) * (1 - t) * known[0]
	+ (1 - r) * (1 - s) * (    t) * known[1]
	+ (1 - r) * (    s) * (1 - t) * known[2]
	+ (1 - r) * (    s) * (    t) * known[3]
	+ (    r) * (1 - s) * (1 - t) * known[4]
	+ (    r) * (1 - s) * (    t) * known[5]
	+ (    r) * (    s) * (1 - t) * known[6]
	+ (    r) * (    s) * (    t) * known[7];
    }

    scalar interpolated_values[19];
    LinearInterpolate19(known, interpolated_values);

    for(int i = 0; i < 19; i++){
      if(interpolated_values[i] != correct_values[i]){
	std::cerr << " Error: interpolated_values[" << i << "] = " << interpolated_values[i] << " != " << correct_values[i] << std::endl;
	test_passed = false;
      }
    }
    if(test_passed){
      std::cerr << " passed\n";
    }


  }
  
  //InitializeChildren()
  {
    std::cerr << "SignedDistanceTransform::InitializeChildren():";
    bool test_passed = true;

    //Initialize the test node
    OctreeNode * current = new OctreeNode();
    current->is_leaf = true;
    current->corner = Vector3d(0,0,0);
    current->edge_length = 2.0;

    //Create some random data for the distance values
    scalar true_distances[19];
    scalar corner_distances[8];
    for(int i = 0; i < 19; i++){ true_distances[i] = rand_scalar(1.0); }
    for(int i = 0; i < 8; i++){ current->distances[i] = corner_distances[i] = rand_scalar(1.0); }
    
    InitializeChildren(current, true_distances);


    /* 
     * Basic tests
     */
    assert(current->is_leaf == false);
    for(int i = 0; i < 8; i++){ assert(current->children[i]->is_leaf == true); }
    for(int i = 0; i < 8; i++){ assert(current->children[i]->edge_length == current->edge_length / 2.0 ); }


    Vector3d correct_corners[8] = {  Vector3d(0,0,0),  Vector3d(0,0,1),  Vector3d(0,1,0), Vector3d(0,1,1),
				     Vector3d(1,0,0),  Vector3d(1,0,1),  Vector3d(1,1,0), Vector3d(1,1,1)};
    for(int i = 0; i < 8; i++){
      if(current->children[i]->corner != correct_corners[i]){
	std::cerr << " Error: children[" << i << "]->corner = " << current->children[i]->corner << " != " << correct_corners[i] << std::endl;
      }
    }


    /*
     * Advanced tests
     */ 
    int correct_distances_index27[8][8] = {
      { 0, 1, 3, 4, 9, 10, 12, 13 },
      { 1, 2, 4, 5, 10, 11, 13, 14 },
      { 3, 4, 6, 7, 12, 13, 15, 16 },
      { 4, 5, 7, 8, 13, 14, 16, 17 },
      { 9, 10, 12, 13, 18, 19, 21, 22},
      { 10, 11, 13, 14, 19, 20, 22, 23},
      { 12, 13, 15, 16, 21, 22, 24, 25},
      { 13, 14, 16, 17, 22, 23, 25, 26}};

    scalar correct_distances27[27];
    correct_distances27[0] = corner_distances[0];
    correct_distances27[2] = corner_distances[1];
    correct_distances27[6] = corner_distances[2];
    correct_distances27[8] = corner_distances[3];
    correct_distances27[18] = corner_distances[4];
    correct_distances27[20] = corner_distances[5];
    correct_distances27[24] = corner_distances[6];
    correct_distances27[26] = corner_distances[7];

    correct_distances27[1] = true_distances[0];
    correct_distances27[3] = true_distances[1];
    correct_distances27[4] = true_distances[2];
    correct_distances27[5] = true_distances[3];
    correct_distances27[7] = true_distances[4];

    correct_distances27[9] = true_distances[5];
    correct_distances27[10] = true_distances[6];
    correct_distances27[11] = true_distances[7];
    correct_distances27[12] = true_distances[8];
    correct_distances27[13] = true_distances[9];
    correct_distances27[14] = true_distances[10];
    correct_distances27[15] = true_distances[11];
    correct_distances27[16] = true_distances[12];
    correct_distances27[17] = true_distances[13];

    correct_distances27[19] = true_distances[14];
    correct_distances27[21] = true_distances[15];
    correct_distances27[22] = true_distances[16];
    correct_distances27[23] = true_distances[17];
    correct_distances27[25] = true_distances[18];


    for(int i = 0; i < 8; i++){
      for(int j = 0; j < 8; j++){
	scalar expected = correct_distances27[correct_distances_index27[i][j]];
	scalar found = current->children[i]->distances[j];
	if(found != expected){
	  std::cerr << " Error: children[" << i << "]->distances[" << j << "] = " << found << " != " << expected << std::endl;
	  test_passed = false;
	}
      }
    }

    if(test_passed){  std::cerr << " passed\n"; }
  }

    

}
