/**
* @file OneSidedRepulsion.h
* @date Jan 23, 2006
* @author John Hart
*/

#ifndef ONESIDEDREPULSION_H
#define ONESIDEDREPULSION_H

#include "Particles.h"
#include "ParticleBehavior.h"
#include "SurfaceAdhesion.h"
#include "AdaptiveRepulsionData.h"
class ParticleLocality;
class ParticlePosition;
class ParticleVelocity;

#define REPULSION_NEIGHBORS 6

/**
* repels based on position of other particles, but does not apply a force to the other particles.
*
* Required Attributes: AdaptiveRepulsionData, ParticleLocality.
* Holds two parameters. The parameter beta is a fudge added to the denominator
* in case it is zero. The parameter rho is a factor controlling the rate at
* which energy approaches a desired energy level.
*/
class OneSidedRepulsion : public ParticleBehavior
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
