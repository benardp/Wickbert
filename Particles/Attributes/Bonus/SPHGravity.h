/*
* @file SPHGravity.h
* @author Yuan Zhou
*/
#ifndef SPHGRAVITY_H
#define SPHGRAVITY_H

#include "Particles.h"
#include "ParticleBehavior.h"
#include "ParticleAcceleration.h"

class SPHGravity : public ParticleBehavior
{
public:
	ParticleAcceleration *acceleration;

	// create SPH attribute for Particles p, default constructor
	SPHGravity(Particles *ps = NULL);

	MAKE_PARTICLESTUFF_NAME();
	
	virtual void attachAttributes();

	virtual void applyForce();
};

#endif
