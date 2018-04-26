/*
@file Sphere.h
@author Wen Su
This class draws a disk for a particle object.
*/

#ifndef PARTICLESHADERSPHERE_H
#define PARTICLESHADERSPHERE_H

#include "ParticleShader.h"

class ParticleShaderSphere:public ParticleShader
{
public:

	GLUquadricObj* quad;
	double radius;
	int sides;

public:

	MAKE_PARTICLESTUFF_NAME();
	
	ParticleShaderSphere(Particles *ps=NULL);

	virtual ~ParticleShaderSphere();

	virtual void drawParticle(int i);

};

#endif
