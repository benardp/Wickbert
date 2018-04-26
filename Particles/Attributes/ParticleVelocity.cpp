/**
 * @file ParticleVelocity.cpp
 * @author Wen Su
 */

#include "ParticleVelocity.h"

REGISTER_PARTICLESTUFF(ParticleVelocity,"Attribute:ParticleVelocity");

/**
 * Add particle age to a system of particles.
 * @param ps   The owning particle system.
 * @param name The name of this object.
 */
ParticleVelocity::ParticleVelocity(Particles *ps, const std::string& name)
	: ParticleAttribute(ps, name) 
{
	if (ps)
		v.resize(ps->size());

	new PSParamgmVector3PerParticle(this,&v,"vel","velocity",
		"Change in position per unit time");
}

/**
 * Add a normal corresponding to the new particle.
 * @param i Index of the new particle.
 */
void ParticleVelocity::particleAdded() 
{
	v.push_back(gmVector3());
}

/**
 * Callback for particle removal.
 * @param i Index of particle to be removed.
 * @see Particles::particleRemoved
 */
void ParticleVelocity::particleRemoved(unsigned int i) 
{
	v[i] = v.back();
	v.pop_back();
}

int ParticleVelocity::qlenpp()
{
	return 3; 
}
void ParticleVelocity::getqpp(double *q, int i)
{
	q[0] = v[i][0];
	q[1] = v[i][1];
	q[2] = v[i][2];
}
void ParticleVelocity::setqpp(double *q, int i)
{
	v[i][0] = q[0];
	v[i][1] = q[1];
	v[i][2] = q[2];
}
void ParticleVelocity::qnamepp(char **qn)
{
	qn[0] = "Velocity x";
	qn[1] = "Velocity y";
	qn[2] = "Velocity z";
}
