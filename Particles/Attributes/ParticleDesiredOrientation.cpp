/**
 * Implementation of ParticleDesiredOrientation.
 * @file ParticleDesiredOrientation.cpp
 * @author Mike Flavin (based off of work by Ed Bachta)
 */

#include "ParticleDesiredOrientation.h"

REGISTER_PARTICLESTUFF(ParticleDesiredOrientation,"Attribute:ParticleDesiredOrientation");

/**
 * Add particle orientation to a system of particles.
 * If particles exist in the system, default their normals
 * to (0.0,0.0,1.0).
 * @param ps   The owning particle system.
 * @param name The name of this object.
 */
ParticleDesiredOrientation::ParticleDesiredOrientation(Particles *ps,const std::string& name) 
  : ParticleAttribute(ps, name)
{
	setParticleSystem(ps);

	new PSParamgmVector3PerParticle(this,&n,
		"initdir","initial orientation","Initial orientation of each particle");
}

void ParticleDesiredOrientation::setParticleSystem(Particles *new_ps)
{
	ParticleAttribute::setParticleSystem(new_ps);

	if (ps)
	{
		n.resize(ps->size());
		gradMag2.resize(ps->size());
		for (unsigned int i=0;i<ps->size();i++) {
			n[i] = gmVector3(0.0,0.0,0.0);
			gradMag2[i] = 0.0;
		}
	}
}

void ParticleDesiredOrientation::clear()
{
	n.clear();
	gradMag2.clear();
}


/**
 * Add a normal corresponding to the new particle.
 * @param i Index of the new particle.
 */
void ParticleDesiredOrientation::particleAdded()
{
	n.push_back(gmVector3());
	gradMag2.push_back(0.0);
}

/**
 * Callback for particle removal.
 * @param i Index of particle to be removed.
 * @see Particles::particleRemoved
 */
void ParticleDesiredOrientation::particleRemoved(unsigned int i)
{
	n[i] = n.back();
	n.pop_back();
	gradMag2[i] = gradMag2.back();
	gradMag2.pop_back();
}
