/*
@file ParticleShaderDisk.h
@author Wen Su
This class draws a disk for a particle object.
*/

#ifndef PARTICLESHADERDISK_H
#define PARTICLESHADERDISK_H

#include "ParticleShader.h"
#include "ParticleAttribute.h"

class ParticleShaderDisk: public ParticleShader
{

protected:
	GLUquadricObj* quad;
	double radius;
	double scale;
	int sides;

	std::string radius_source;
	ParticleAttribute *radius_attr;
	DoubleVector *radius_data;

public:

	MAKE_PARTICLESTUFF_NAME();

	ParticleShaderDisk(Particles *ps=NULL);
	virtual void attachAttributes();
   
	virtual ~ParticleShaderDisk();

	virtual void drawParticle(int i);
};

#endif
