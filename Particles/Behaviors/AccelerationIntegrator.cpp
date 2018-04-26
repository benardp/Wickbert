/**
 * @file AccelerationIntegrator.cpp 
 * @author Yuan Zhou
 */

#include "AccelerationIntegrator.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"
#include "ParticleAcceleration.h"

REGISTER_PARTICLESTUFF(AccelerationIntegrator,"Behavior:AccelerationIntegrator");

/**
 * Creates a ParticleFate object. First it looks for the required 
 * attributes in the particle system, and then it sets default 
 * constants.
 */
AccelerationIntegrator::AccelerationIntegrator(Particles *ps, const std::string& name) 
	: ParticleBehavior(ps, name)
{
	new PSParamDouble(this,&time_step,0.03,"dt","time step",
		"Smaller values provide more accurate integration of"
		"acceleration into velocity but require more computation");
	new Attached<ParticlePosition>(this,&pos);
	new Attached<ParticleVelocity>(this,&vel);
	new Attached<ParticleAcceleration>(this,&accel);
}

void AccelerationIntegrator::integrate()
{
	if (!accel || !vel || !pos) return;

	for(unsigned int i=0;i<ps->size();i++) {
		pos->x[i] += time_step * vel->v[i];
		vel->v[i] += time_step * accel->acc[i];
	}
}
