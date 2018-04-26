/**  @file DistanceField.h
 *  @author Matei Stroila (using Nathan Bell's distance field code)
 *  @date 10/7/05
 */

#ifndef DistanceField_h
#define DistanceField_h

#include "Surface/Implicit/Implicit.h"

//Distance transform
#include "signeddistancetransform.h"
#include "cachingsparsescalarlattice.h"

/** Distance field is an Implicit designed to accelerate the distance of a
 *  query point to a surface mesh using a uniform grid.
 */
class DistanceField : public Implicit
{
public:
	DistanceField();
	DistanceField(std::string&); ///< Explicit constructor if we know the file name of the obj file.
	
	double proc(const gmVector3 & x);
	gmVector3 grad(const gmVector3 & x);
	
	void setMesh(void);

	~DistanceField(); 
	MAKE_NAME();
	
private:
		
	void setSignedDistanceTransfrom();
	
	std::string filename;
	double voxel_size;
	double max_distance;
	TriangleMesh2 *mesh;
	SignedDistanceTransform *dist_transform;
	SparseScalarLattice *sparse_lattice;
	scalar max_value;
};

class DistanceFieldLoadOBJ : public SurfParamButton::Callback
{
public:
	DistanceField *df;
	DistanceFieldLoadOBJ(DistanceField *me) {df = me;}
	virtual void onbuttonpress() {df->setMesh();}
};

#endif

