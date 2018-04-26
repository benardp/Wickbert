/**
 * Declaration of ParticleValueDeath.h
 * @file ParticleValueDeath.h
 * @date May 5, 2004
 * @author Mike Flavin
 */

#ifndef __PARTICLEVALUEDEATH_H__
#define __PARTICLEVALUEDEATH_H__

#include "ParticleBehavior.h"
class ImplicitInterrogator;

/**
 * This class is very simple - it checks if the value of the function at a floater location is
 *  above a certain threshold, and if it is, it kills the particle.  This is really only useful
 *  for implicits where we may have particles that are, for some reason, constantly moving AWAY
 *  from the surface.  Right now I intend to use it only to help with Compact RBFs, but I
 *  suppose that it could be useful in other situations as well.
 */
class ParticleValueDeath : public ParticleBehavior
{
public:
	double threshold; /// The threshold function value at which to kill the particles.
	ImplicitInterrogator* imp_int;

	ParticleValueDeath(Particles *ps=NULL, const std::string& name=std::string("ParticleValueDeath"));

	/// We determine particle fate during the cleanup step.
	void cleanup();

	MAKE_PARTICLESTUFF_NAME();
};

#endif