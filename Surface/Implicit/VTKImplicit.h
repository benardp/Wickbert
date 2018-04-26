/*
 *  VTKImplicit.h
 *  gps
 *
 *  Created by Matei Stroila on 11/4/05.
 *  
 *
 */

#ifndef VTKImplicit_h
#define VTKImplicit_h

#include "Surface/Implicit/Implicit.h"
#include "Particles.h"
#include "ParticlePosition.h"
#include "ParticleBoundingBox.h"
#include "ParticleSystem.h"

#include "vtkImplicitVolume.h"
#include "vtkImplicitModeller.h"
#include "vtkOBJReader.h"
#include "vtkPolyData.h"

class VTKImplicit : public Implicit
{
public:
	VTKImplicit();	
	double proc(const gmVector3 & x);
	gmVector3 grad(const gmVector3 & x);
	
	void setMesh(void);
	
	~VTKImplicit(){}; 
	MAKE_NAME();
	
private:
	
	vtkImplicitVolume* impVol;
	vtkImplicitModeller* impModeller;
	vtkOBJReader* reader;
	vtkPolyData* mesh;
	std::string filename;
	Particles *pVertices;
	ParticlePosition* pos;
	ParticleBoundingBox* bbox;
	bool attachParticles();
};

class VTKImplicitLoadOBJ : public SurfParamButton::Callback
{
public:
	VTKImplicit *vtkImp;
	VTKImplicitLoadOBJ(VTKImplicit *me) {vtkImp = me;}
	virtual void onbuttonpress() {vtkImp->setMesh();}
};


#endif
