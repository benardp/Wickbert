/**
 * Implementation of ParticleAge.
 * @file ParticleAge.cpp
 * @author Doug Nachand
 */

#include "ParticleAge.h"

REGISTER_PARTICLESTUFF(ParticleAge,"Attribute:ParticleAge");

/**
 * Add particle age to a system of particles.
 * @param ps   The owning particle system.
 * @param name The name of this object.
 */
ParticleAge::ParticleAge(Particles *ps, const std::string& name)
	: ParticleAttribute(ps, name) 
{
	new PSParamIntPerParticle(this,&t,"age","particle age","Current lifetime of particle");
	new PSParamInt(this,&maxLifeTime,100,"maxLifeTime","Maximum Life Time", "Maximum Life Time");
	new PSParamBool(this,&useParticleAge,true,"useParticleAge","useParticleAge", "Set true to use the aging");
}

void ParticleAge::reset()
{
	for(unsigned int i=0; i<t.size(); i++){
		t[i] = 0;
		resetFlag[i] = true;
	}
}

void ParticleAge::resetParticleAge(const int i, const bool b)
{
	resetFlag[i] = resetFlag[i] * (int) b;
}


void ParticleAge::clear()
{ 
	t.clear();
	resetFlag.clear();
}

/**
 * Add a normal corresponding to the new particle.
 * @param i Index of the new particle.
 */
void ParticleAge::particleAdded() 
{
	t.push_back(0);
	resetFlag.push_back(true);
}

/**
 * Callback for particle removal.
 * @param i Index of particle to be removed.
 * @see Particles::particleRemoved
 */
void ParticleAge::particleRemoved(unsigned int i) 
{
	t[i] = t.back();
	t.pop_back();

	resetFlag[i] = resetFlag.back();
	resetFlag.pop_back();
}

void ParticleAge::prepare()
{
	if(!useParticleAge) return;

	for(unsigned int i=0; i<t.size(); i++){
		if(resetFlag[i]){
			t[i] = 0;
		}
		t[i]++;
		resetFlag[i] = true;
	}
}