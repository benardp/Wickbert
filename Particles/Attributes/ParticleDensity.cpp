/**
 * @file ParticleDensity.cpp
 * @author Yuan Zhou
 */

#include "ParticleDensity.h"

REGISTER_PARTICLESTUFF(ParticleDensity,"Attribute:ParticleDensity");

/**
 * Add particle age to a system of particles.
 * @param ps   The owning particle system.
 * @param name The name of this object.
 */
ParticleDensity::ParticleDensity(Particles *ps, const std::string& name)
	: ParticleAttribute(ps, name) 
{
	if (ps)
		den.resize(ps->size());

	new PSParamDoublePerParticle(this,&den,
		"density","particle density","Density of material represented by particle.");
}

/**
 * Add a normal corresponding to the new particle.
 * @param i Index of the new particle.
 */
void ParticleDensity::particleAdded() 
{
	den.push_back(double());
}

/**
 * Callback for particle removal.
 * @param i Index of particle to be removed.
 * @see Particles::particleRemoved
 */
void ParticleDensity::particleRemoved(unsigned int i) 
{
	den[i] = den.back();
	den.pop_back();
}




/// this gives a calling convertion for children to get the density of index i
double ParticleDensity::getDensity(unsigned int i)
{
	return den[i];
}

/// this gives a calling convertion for children to set the density of index i
void ParticleDensity::setDensity(unsigned int i, double d)
{
	den[i]=d;
}

void ParticleDensity::clear()
{
	den.clear();
}
