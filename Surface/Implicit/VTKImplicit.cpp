/*
 *  VTKImplicit.cpp
 *  gps
 *
 *  Created by Matei Stroila on 11/4/05.
 *  
 *
 */

#ifdef WB_USE_VTK
#include "VTKImplicit.h"
#include <vtkstd/exception>

REGISTER_IMPLICIT(VTKImplicit,"VTKImplicit");

VTKImplicit::VTKImplicit()
{
	new SurfParamString(this,&filename,"/home/users/stroila2/Data/cow-800.obj","Mesh file",
						"OBJ Mesh file");
	new SurfParamButton(this,new VTKImplicitLoadOBJ(this),"load","Load Mesh", "Load the triangular mesh");

	pVertices = NULL;
	pos = NULL;
	bbox = NULL;
	new SurfParRefParam(this,&pVertices,"floaters","pCenters","Centers", "Collection of particle whose positions provide the center locations for the mesh vertices.");
	reader = vtkOBJReader::New();
	impVol = vtkImplicitVolume::New();
	impModeller = vtkImplicitModeller::New();
}

double VTKImplicit::proc(const gmVector3 & x)
{
	double val = impVol->EvaluateFunction(x[0], x[1], x[2]);	
	return val;
}

gmVector3 VTKImplicit::grad(const gmVector3 & x)
{
	double n[3];
	double px[3];
	for(size_t i = 0; i < 3; i++) px[i] = x[i];
	impVol->EvaluateGradient(px,n);
	return gmVector3(n[0],n[1],n[2]);	
}
	

void VTKImplicit::setMesh(void)
{
	reader->SetFileName(filename.c_str());
	mesh = reader->GetOutput();
	try
	{
		reader->Update();
	}
    catch (vtkstd::exception& e)
	{
		std::cerr << "Unable to read OBJ file " << e.what();
		return;
	}
	
	if(!attachParticles()) return;
	
	impModeller->SetInputConnection(reader->GetOutputPort());
	impModeller->Update();
	
	impModeller->ComputeModelBounds();
	double bounds[6];
	impModeller->GetModelBounds(bounds);	
	
	bbox->setq(bounds);
	
	impModeller->Update();
	impVol->SetVolume(impModeller->GetOutput());
}

bool VTKImplicit::attachParticles()
{
	vtkPoints* vtkpoints = mesh->GetPoints();
	if(!vtkpoints) return false;
	int numPoints = vtkpoints->GetNumberOfPoints();
	
	pos = pVertices->getAttribute<ParticlePosition>(std::string("ParticlePosition"));
	if (!pos){ 	
		return false;
	}

	bbox = pVertices->getAttribute<ParticleBoundingBox>(std::string("ParticleBoundingBox"));
	if (!bbox){ 		
		return false;
	}

	double x[3];
	int newparticle;
	for( int n = 0; n < numPoints; n++)
	{
		vtkpoints->GetPoint(n, x);
		newparticle = pVertices->addParticle();
		gmVector3 posn(x[0],x[1],x[2]);
		(pos->x).push_back(posn);
	}
	return true;
	
}


#endif
