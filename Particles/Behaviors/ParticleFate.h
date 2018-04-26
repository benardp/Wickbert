/**
 * Declaration of ParticleFate.h
 * @file ParticleFate.h
 * @date November 28, 2001
 * @author Ed Bachta
 */

#ifndef PARTICLEFATE_H
#define PARTICLEFATE_H

#include "ParticleBehavior.h"

class AdaptiveRepulsionData;
class ParticleNormal;
class ParticleBoundingBox;
class ParticlePosition;
class ParticleVelocity;
class ImplicitInterrogator;
class ParticleAge;

/**
 * This is the class which governs particle birth and death as described by
 * Witkin and Heckbert. Particles which are large and have low energy
 * birth new particles and particles which are overcroweded are killed.
 */
class ParticleFate : public ParticleBehavior
{
	AdaptiveRepulsionData* rep_data; ///< Witkin-Heckbert repulsion data.
	ParticleNormal* p_orient;   ///< Particle orientations.
	ParticleBoundingBox* p_bounds;   ///< Bounding box.
	ParticlePosition* position;
	ParticleVelocity* velocity;
	ParticleAge *pAge;
public:
	MAKE_PARTICLESTUFF_NAME();

	ParticleFate(Particles *ps=NULL, const std::string &name=std::string("ParticleFate"));

	double gamma; ///< Equilibrium speed (multiple of repulsion radius).
	double nu;    ///< Fraction of E_hat, for fissioning.
	double delta; ///< Fraction of sigma_hat, for death.
	int population;

	ImplicitInterrogator *impint;	///< Pointer to ImplicitInterrogator.
	double fateval;				///< Abs value threshold for fission/death.

	void attachAttributes();

	/// Allows for an awareness of surface diameter.
	double setSurfaceDiameter(double d);

	/// Allow the user to specify the desired radius.
	void setDesiredRadius(double &r); 

	/// We determine particle fate during the cleanup step.
	void cleanup();

};

#endif


