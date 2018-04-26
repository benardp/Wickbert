/**
 * @file VelocityIntegrator.cpp 
 * @author Wen Su
 */

#include "VelocityIntegrator.h"
#include "ParticleMesh.h"

REGISTER_PARTICLESTUFF(VelocityIntegrator,"Behavior:VelocityIntegrator");

/**
 * Creates a ParticleFate object. First it looks for the required 
 * attributes in the particle system, and then it sets default 
 * constants.
 */
VelocityIntegrator::VelocityIntegrator(Particles *ps, const std::string& name) 
	: ParticleBehavior(ps, name)
{
}

void VelocityIntegrator::integrate()
{
	for(unsigned int i=0;i<ps->size();i++)
	{
		position->setPosition(i,position->getPosition(i)+ps->dt*velocity->v[i]);
		velocity->v[i] = gmVector3();
	}
}
