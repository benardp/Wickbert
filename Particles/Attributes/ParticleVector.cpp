/**
* Implementation of the attribute representing per particle vector data.
* @file ParticleVector.cpp
* @date 27 June. 2005
* @author Matei N. Stroila
*/

#include "ParticleVector.h"


REGISTER_PARTICLESTUFF(ParticleVector,"Attribute:ParticleVector");

/**
* Add the attribute
 * @param ps   The owning particle system.
 * @param name The name of this object.
 */
ParticleVector::ParticleVector(Particles *ps, const std::string& name)
: ParticleAttribute(ps, name) 
{
	if (ps)
		data.resize(ps->size());	
	new PSParamgmVector3PerParticle(this,&data,"data","data","Particle vector data");
}

/**
* Add data corresponding to the new particle.
 * @param i Index of the new particle.
 */
void ParticleVector::particleAdded() 
{
	data.push_back(gmVector3());
}

/**
* Callback for particle removal.
 * @param i Index of particle to be removed.
 * @see Particles::particleRemoved
 */
void ParticleVector::particleRemoved(unsigned int i) 
{
	data[i] = data.back();
	data.pop_back();
}


void ParticleVector::clear()
{
	data.clear();
}

int ParticleVector::qlenpp()
{
	return 3; 
}
void ParticleVector::getqpp(double *q, int i)
{
	q[0] = getVector(i)[0];
	q[1] = getVector(i)[1];
	q[2] = getVector(i)[2];
}
void ParticleVector::setqpp(double *q, int i)
{
	data[i][0] = q[0];
	data[i][1] = q[1];
	data[i][2] = q[2];
}
void ParticleVector::qnamepp(char **qn)
{
	qn[0] = "Vector x";
	qn[1] = "Vector y";
	qn[2] = "Vector z";
}
