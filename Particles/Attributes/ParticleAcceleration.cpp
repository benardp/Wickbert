/**
 * @file ParticleAcceleration.cpp
 * @author Yuan Zhou
 */

#include "ParticleAcceleration.h"

REGISTER_PARTICLESTUFF(ParticleAcceleration,"Attribute:ParticleAcceleration");

/**
 * Add particle age to a system of particles.
 * @param ps   The owning particle system.
 * @param name The name of this object.
 */


ParticleAcceleration::ParticleAcceleration(Particles *ps, const std::string& name)
	: ParticleAttribute(ps, name) 
{
	if (ps)
		acc.resize(ps->size());

	new PSParamgmVector3PerParticle(this,&acc,
		"accel","acceleration","The acceleration of an individual particle");
}

/**
 * Add a normal corresponding to the new particle.
 * @param i Index of the new particle.
 */
void ParticleAcceleration::particleAdded() 
{
	acc.push_back(gmVector3());
}

/**
 * Callback for particle removal.
 * @param i Index of particle to be removed.
 * @see Particles::particleRemoved
 */
void ParticleAcceleration::particleRemoved(unsigned int i) 
{
	acc[i] = acc.back();
	acc.pop_back();
}

/// this gives a calling convertion for children to get the acceleration of index i
gmVector3 ParticleAcceleration::getAcceleration(unsigned int i)
{
	return acc[i];
}

/// this gives a calling convertion for children to set the acceleration of index i
void ParticleAcceleration::setAcceleration(unsigned int i, gmVector3 p)
{
	acc[i]=p;
}
