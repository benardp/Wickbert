/**  @file DistanceField.cpp
 *  @author Matei Stroila
 *  @date 10/7/05
 */

#include "DistanceField.h"

REGISTER_IMPLICIT(DistanceField,"DistanceField");

DistanceField::DistanceField(std::string& _filename)
{		
	filename = _filename;
	setMesh();
}

DistanceField::DistanceField(void)
{
	dist_transform = NULL;
	sparse_lattice = NULL;
	mesh = NULL;
	
	new SurfParamString(this, &filename,"<empty>","Mesh File", "Obj Mesh File");
	new SurfParamDouble(this, &voxel_size, 0.1, "voxel_size","Grid Size");
	new SurfParamDouble(this, &max_distance, 10.0, "max_distance", "Max Dist", "Maximum distance from the zero-isosurface (in grid coordinates)");
	new SurfParamButton(this,new DistanceFieldLoadOBJ(this),"load","Load Mesh", "Load the triangular mesh->");
}

double DistanceField::proc(const gmVector3 & x)
{	
	scalar interpolated_value;
	Vector3d test_point(x[0],x[1],x[2]);
	if(sparse_lattice->ContainsPoint(test_point))
	   interpolated_value = sparse_lattice->InterpValue(test_point);
	else
	   interpolated_value = -max_value;

	return interpolated_value;
}

gmVector3 DistanceField::grad(const gmVector3 & x)
{
	Vector3d interpolated_vector(-0.1,-0.1,-0.1);
	Vector3d test_point(x[0],x[1],x[2]);
	if(sparse_lattice->ContainsPoint(test_point))
		interpolated_vector = sparse_lattice->InterpNormal(test_point);
	
	return gmVector3(interpolated_vector[0],interpolated_vector[1],interpolated_vector[2]);
}


void DistanceField::setMesh(void)
{
	//load an OBJ mesh->  TRIANGLE MESHES ONLY!
	if(mesh) delete mesh;
	mesh = new TriangleMesh2(filename.c_str());
	if (!mesh->IsLoaded()) return;
	//scale the mesh
	Vector3d bbox_size = mesh->BoundingBox().v1 - mesh->BoundingBox().v0;
	double scale_factor = bbox_size[0];
	for(int i = 1; i < 3; i++) if(bbox_size[i] > scale_factor) scale_factor = bbox_size[i];
	scale_factor = 2.0/scale_factor;
	mesh->Scale(Vector3d(scale_factor,scale_factor,scale_factor));
	setSignedDistanceTransfrom();
}

void DistanceField::setSignedDistanceTransfrom()
{
	/*
	 * Setup a signed distance transform object.  This precomputes
	 * some angle weighted normals used later to determine sign
	 */
	if(dist_transform) delete dist_transform;
	dist_transform = new SignedDistanceTransform(*mesh);
	
	/*
	 * SparseScalarLattice is a hash data structure that maps integer grid 
	 * points to scalar values.  The only argument to the constructor is the
	 * grid size (voxel edge length).
	 *
	 * Underlying container is an STL HashMap<Cell,scalar>
	 *
	 */
	if(sparse_lattice) delete sparse_lattice;
	sparse_lattice = new SparseScalarLattice((scalar) voxel_size);
	
	/*
	 *  Compute the distance transform on the sparse lattice to a given distance from
	 *  the zero-isosurface.  When using trilinear interpolation for normals or scalar values
	 *  make sure that this distance is sufficiently large to accomodate your queries.
	 *
	 */
	
	dist_transform->Voxelize(*sparse_lattice, (scalar) max_distance * sparse_lattice->VoxelSize() ); 
	
	/*
	 *  Output some info to illustrate how to access the scalar lattice
	 */
	std::cout << "Number of values in the lattice: " << sparse_lattice->size() << std::endl;
	std::cout << std::endl;
	
	max_value = SCALAR_MIN;
	scalar min_value = SCALAR_MAX;
	for(SparseScalarLattice::const_iterator iter = sparse_lattice->begin(); iter != sparse_lattice->end(); iter++){
		max_value = std::max(max_value, iter->second);
		min_value = std::min(min_value, iter->second);
	}
	
	std::cout << "Maximum value: " << max_value << std::endl;
	std::cout << "Minimum value: " << min_value << std::endl;
	std::cout << std::endl;
	
	/*
	 *  Scale all the values by a factor of two
	 */
	//for(SparseScalarLattice::iterator iter = sparse_lattice->begin(); iter != sparse_lattice->end(); iter++){
	//	iter->second *= 2.0;
	//}
	
	
	/*
	 *  Output a few elements of the lattice (remember, hash map is unordered)
	 */
	int i = 0;
	for(SparseScalarLattice::const_iterator iter = sparse_lattice->begin(); 
		i < 5; 
		i++, iter++){
		const Cell&   current_cell  = iter->first;
		const scalar& current_value = iter->second;
		
		std::cout << "Grid point [" <<  current_cell << "]  has value  " << current_value << std::endl;
		//Could also use sparse_lattice.Value(current_cell)    
		
		assert(current_value == sparse_lattice->Value(current_cell));
	}
	std::cout << std::endl;
	
	
	/*
	 *  Use trilinear interpolation on a point in R^3, in this case we'll use a vertex on the pig
	 *  which is guaranteed to lie within the lattice (i.e. all 8 corners of voxel containing
													   *  this point are in the lattice)
	 *
	 *  Since the true distance at this point is 0 (it lies on the mesh surface) we expect the 
	 *  interpolated value to be near zero.
	 * 
	 *  InterpNormal() and InterpMeanCurvature() work similarly
	 *
	 *  For better performance, use CachingScalarLattice().  This is much faster when making lots of 
	 *  Interp queries on a fixed lattice.
	 *   
	 */
	Vector3d test_point = mesh->vertices[0];
	scalar interpolated_value = sparse_lattice->InterpValue(test_point);
	
	std::cout << "Test Point (" <<  test_point << ") interpolated value at test point: " << interpolated_value << std::endl;
	std::cout << std::endl;
	
}


DistanceField::~DistanceField(){
	
	if(dist_transform) delete dist_transform;
	if(sparse_lattice) delete sparse_lattice;
	if(mesh) delete mesh;
} 
