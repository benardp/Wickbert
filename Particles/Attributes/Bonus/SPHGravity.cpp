/**
* @file SPHGravity.cpp
* @author Yuan Zhou
*/

#include "SPHGravity.h"
#include "ParticleAcceleration.h"

REGISTER_PARTICLESTUFF(SPHGravity,"Behavior:SPHGravity");

SPHGravity::SPHGravity(Particles *ps) 
: ParticleBehavior(ps, std::string("SPHGravity") )
{
	acceleration = NULL;
}

void SPHGravity::attachAttributes()
{
	ParticleBehavior::attachAttributes();
	attachAttribute(acceleration,std::string("ParticleAcceleration"));
}

void SPHGravity::applyForce()
{
	unsigned int i;
	for(i=0; i<ps->size(); i++)
	{
		acceleration->acc[i] += gmVector3(0.0,0.0, -9.8); 
	}
}