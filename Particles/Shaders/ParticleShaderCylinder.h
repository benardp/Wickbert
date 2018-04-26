/*
@file ParticleShaderCylinder.h
@author Wen Su
This class draws a disk for a particle object.
*/

#ifndef PARTICLESHADERCYLINDER_H
#define PARTICLESHADERCYLINDER_H

#include "ParticleShader.h"

class ParticleShaderCylinder:public ParticleShader
{
protected:

	GLUquadricObj* quad;
	float base;
	float top;
	float height;
	int sides;
	int capped;

public:

	MAKE_PARTICLESTUFF_NAME();
	
	ParticleShaderCylinder(Particles *ps=NULL);

	virtual ~ParticleShaderCylinder();

	virtual void drawParticle(int i);

};

#endif
