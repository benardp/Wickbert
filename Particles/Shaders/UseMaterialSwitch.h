/*
@file UseMaterialSwitch.h
@author Tony Kaap
*/

#ifndef USEMATERIALSWITCH_H
#define USEMATERIALSWITCH_H

#include "ParticleShader.h"
#include "ParticleScalar.h"
#include "ParticleMaterial.h"

// This class allows for a particle system to be shaded
// by three different materials according to some
// per-particle value stored in a ParticleScalar attribute
class UseMaterialSwitch: public ParticleShader
{

protected:
	ParticleScalar *switchScalar;
	ParticleMaterial *materialNegative;
	ParticleMaterial *materialZero;
	ParticleMaterial *materialPositive;

public:

	MAKE_PARTICLESTUFF_NAME();

	UseMaterialSwitch(Particles *ps=NULL);

	virtual void drawParticle(int i);
};

#endif
