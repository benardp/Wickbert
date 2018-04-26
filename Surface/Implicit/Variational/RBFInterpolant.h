//RBFInterpolant.h
//Loads a mesh and interpolates it using radial basis functions placed at each of the vertices
//Mike Mullan 

#ifndef RBFINTERPOLANT_H
#define RBFINTERPOLANT_H

#include "Surface/Implicit/Variational/RBF.h"
#include "Surface/OpenMesh/TriMesh.h"
#include <stdio.h>

#define epsilon .01

class RBFMesh : public TriMesh 
{
	public:
	bool readFile(const char* filename=NULL);	
	
};

//A class designed to interpolate polygonal meshes with radial basis functions, and support their level set propagation
class RBFInterpolant : public RBF
{
	
	public:
	
	RBFMesh rbfMesh;  //Mesh we are interpolating
	RBFMesh targetMesh;  //Mesh we want to morph to

	std::vector<gmVector3> normals;

	bool dipoles;  //Are we using the dipole model, or just standard on-surface centers?
	bool morph;   //Are we morphing to another mesh?

	RBFInterpolant() 
	{
		centers.clear();
	}
	
	RBFInterpolant(bool theDipoles) 
	{
		centers.clear();
		dipoles = theDipoles;
	}
	
	//Place RBF centers at each of the vertices
	void placeCenters();
	
	void procq(const gmVector3 & x, double* dfdq);
	void getq(double* q);
	void _setq(double* q);
	unsigned int qlen();
	void getqname(char** qn);
	
	gmMatrix3 computeHessian(gmVector3 point);  //Hessian used for normal vector propagation
	// should be moved to particle
	//void interpolate(Particles *ps, std::valarray<bool> flexible);
	double targetShape(gmVector3 point);  //Distance function to target mesh

	bool moving;  //Is the surface currently being propagated?
};

#endif
