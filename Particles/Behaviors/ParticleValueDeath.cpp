/**
 * Implementation of ParticleValueDeath.
 * @file ParticleValueDeath.cpp 
 * @date May 5, 2004 (updated 22 Dec. 2005)
 * @author Mike Flavin
 */

#include "ParticleValueDeath.h"
#include "ImplicitInterrogator.h"


REGISTER_PARTICLESTUFF(ParticleValueDeath,"Behavior:ParticleValueDeath");

/**
 * Creates a ParticleValueDeath object. 
 */
ParticleValueDeath::ParticleValueDeath(Particles *ps, const std::string& name) 
	: ParticleBehavior(ps, name)
{
	new PSParamDouble(this, &threshold, 5.0, "threshold", "kill threshold",
		"Particles whose absolute value exceeds this threshold are killed.");

	new Attached<ImplicitInterrogator>(this, &imp_int);
}

/**
* We kill particles during the cleanup step.
*/
void ParticleValueDeath::cleanup()
{
	// Iterate backward to preserve indexing while removing particles
	for (int i=ps->size()-1; i>=0; i--) 
	{
		// If function value at current location is above the threshold, kill the particle.
		if (fabs(imp_int->proc(i)) > threshold)
			ps->removeParticle(i);
	}
}
