/**
 * @file ParticlePosition.cpp
 * @author Wen Su
 */

#include "ParticlePosition.h"

REGISTER_PARTICLESTUFF(ParticlePosition,"Attribute:ParticlePosition");

/**
 * Add particle age to a system of particles.
 * @param ps   The owning particle system.
 * @param name The name of this object.
 */
ParticlePosition::ParticlePosition(Particles *ps, const std::string& name, bool xparam)
	: ParticleAttribute(ps, name) 
{
	if (ps)
		x.resize(ps->size());
	changed = true;

    if (xparam)
	    new PSParamgmVector3PerParticle(this,&x,"pos","position","Particle positions");
}

extern double psRandom();

/**
 * Add a normal corresponding to the new particle.
 * @param i Index of the new particle.
 */
void ParticlePosition::particleAdded() 
{
	x.push_back(gmVector3(psRandom(),psRandom(),psRandom()));
	changed = true;
}

/**
 * Callback for particle removal.
 * @param i Index of particle to be removed.
 * @see Particles::particleRemoved
 */
void ParticlePosition::particleRemoved(unsigned int i) 
{
	x[i] = x.back();
	x.pop_back();
	changed = true;
}

/**
 * Returns the position of the ith particle.  Note: cannot be inlined
 * or we lose v-table dependent code.
 */
gmVector3 ParticlePosition::getPosition(unsigned int i)
{
	return x[i];
}


/// this gives a calling conversion for children to set the position of index i
void ParticlePosition::setPosition(unsigned int i, gmVector3 p)
{
	x[i]=p;
	changed = true;
}

void ParticlePosition::clear()
{
	x.clear();
	changed = true;
}
