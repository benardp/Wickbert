/*
@file ParticleShaderDiskVisibility.h
@author Matei Stroila, but subclassed from ParticleShaderDisk by John Hart
This class draws a disk for a particle object if the disk is visible
*/

#ifndef PARTICLESHADERDISKVISIBILITY_H
#define PARTICLESHADERDISKVISIBILITY_H

#include "ParticleShaderDisk.h"
#include "ParticleVisibility.h"
#include "ViewDependence.h"

class ParticleShaderDiskVisibility: public ParticleShaderDisk
{

protected:
	gmVector3 *myCameraPosition;

public:
	MAKE_PARTICLESTUFF_NAME();

	ParticleShaderDiskVisibility(Particles *ps=NULL);
	virtual void attachAttributes();
   
	virtual void drawParticle(int i);
	virtual void drawPre();

	ParticleVisibility *vis;
	ViewDependence *view;
	bool useVisibility;
	bool cameraPosChanged;
	gmVector3 *oldCameraPosition;
};

#endif
