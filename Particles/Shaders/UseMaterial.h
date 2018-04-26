/*
@file UseMaterial.h
@author John C. Hart
*/

#ifndef USEMATERIAL_H
#define USEMATERIAL_H

#include "ParticleShader.h"

class ParticleMaterial;

class UseMaterial: public ParticleShader
{

protected:
	ParticleMaterial *material;

public:

	MAKE_PARTICLESTUFF_NAME();

	UseMaterial(Particles *ps=NULL);

	virtual void drawParticle(int i);
};

#endif
