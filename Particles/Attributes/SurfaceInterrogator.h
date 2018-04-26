/**
* Definition of SurfaceInterrogator.
* @file SurfaceInterrogator.h
* @author Wen Su
*/

#ifndef SURFACEINTERROGATOR_H
#define SURFACEINTERROGATOR_H

#include "Particles.h"
#include "ParticleAttribute.h"

class Surface;

class SurfaceInterrogator : public ParticleAttribute
{	

public:
	
	MAKE_PARTICLESTUFF_NAME();

	Surface* surface;
	std::string surfname;

	// temperary store the index points to the surface for file saving and loading
	//int index;

	SurfaceInterrogator(Particles* ps=NULL, Surface* surface = NULL, const std::string& name = std::string("SurfaceInterrogator"));

	virtual void attachAttributes();

	virtual void setSurface(Surface *s) { surface = s; }
	virtual Surface* getSurface() { return surface; }
};

#endif

