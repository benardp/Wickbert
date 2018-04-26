/**
* Declaration of the repulsion force.
* @file ParticleRepulsion.h
* @date August 31, 2001
* @author Ed Bachta
*/

#ifndef PARTICLEREPULSION_H
#define PARTICLEREPULSION_H

#include "ParticleBehavior.h"
class ParticleLocality;
class ParticlePosition;
class ParticleVelocity;
class AdaptiveRepulsionData;
#define REPULSION_NEIGHBORS 6

/**
* This is the repulsion force described by Witkin and Heckbert S96
* which evenly distributes particles across the zero-set of an
* implicit function. The equations are derived assuming a single
* component of finite area.
*
* Required Attributes: AdaptiveRepulsionData, ParticleLocality.
* Holds two parameters. The parameter beta is a fudge added to the denominator
* in case it is zero. The parameter rho is a factor controlling the rate at
* which energy approaches a desired energy level.
*/
class ParticleRepulsion : public ParticleBehavior
{	
public:
	MAKE_PARTICLESTUFF_NAME();
	
	/// Reference to a repulsion data attribute.
	AdaptiveRepulsionData* rep_data;

	/// Reference to a locality object.
	ParticleLocality* p_locality;

	ParticlePosition* position;
	ParticleVelocity* velocity;

	double beta; ///< Constant to prevent divide-by-zero.
	double rho;  ///< Feedback coefficient.
	
	/// Creates a particle repulsion attribute for Particles p.
	ParticleRepulsion(Particles *ps=NULL);

	/// Applies a force that repels nearby particles.
	void applyForce();
	
	// update repulsion data
	void integrate();
	
	// call update on locality
	void cleanup();	
};

#endif
