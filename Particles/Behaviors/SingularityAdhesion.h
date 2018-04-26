/**
@file SingularityInterrogator.h
@author Wen Su and John C. Hart
Adds a particle force in direction of decreasing gradient magnitude.
*/

#ifndef SINGULARITYADHESION_H
#define SINGULARITYADHESION_H

#include "Particles.h"
#include "ParticleBehavior.h"

class ImplicitInterrogator;
class ParticlePosition;
class ParticleVelocity;

/** Adds a particle force in direction of decreasing gradient magnitude. */
class SingularityAdhesion : public ParticleBehavior
{
public:
	
	MAKE_PARTICLESTUFF_NAME();
	
	/// parameters to control for the user
	double strength;	///< scale applied to grad mag
	double clamp;		///< maximum strength of force


	/// default constructor
	SingularityAdhesion(Particles* ps=NULL, const std::string& name = std::string("SingularityAdhesion"));
	
	ImplicitInterrogator *impInt;
	ParticlePosition *position;
	ParticleVelocity *velocity;

	// move opposite of gradient
	virtual void applyForce();
};

#endif

