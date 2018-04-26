/**
* @file MeshInterrogator.h
* @author John C. Hart
* @date 4 Jan. 2005
*/

#ifndef MESHINTERROGATOR_H
#define MESHINTERROGATOR_H

#include "SurfaceInterrogator.h"
class SurfaceMesh;

class MeshInterrogator : public SurfaceInterrogator
{	

public:
	
	MAKE_PARTICLESTUFF_NAME();

	// temperary store the index points to the surface for file saving and loading
	//int index;

	MeshInterrogator(Particles* ps=NULL, SurfaceMesh* ms = NULL, const std::string& name = std::string("MeshInterrogator"));

	virtual void attachAttributes();

	virtual void setSurface(Surface *s);
	virtual Surface* getSurface() { return surface; }

	virtual void setMesh(SurfaceMesh *ms);
	virtual SurfaceMesh* getMesh();
};

#endif

