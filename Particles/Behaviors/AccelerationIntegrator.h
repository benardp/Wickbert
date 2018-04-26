/**
 * @file AccelerationIntegrator.h
 * @author Yuan Zhou
 */

#ifndef  ACCELERATIONINTEGRATOR_H
#define ACCELERATIONINTEGRATOR_H

#include "Particles.h"
#include "ParticleBehavior.h"
class ParticlePosition;
class ParticleVelocity;
class ParticleAcceleration;

class AccelerationIntegrator : public ParticleBehavior
{
public:
	MAKE_PARTICLESTUFF_NAME();

	double time_step;
	ParticleAcceleration *accel;
	ParticleVelocity *vel;
	ParticlePosition *pos;

	AccelerationIntegrator(Particles *ps=NULL, const std::string& name=std::string("AccelerationIntegrator"));

	void integrate();
};

#endif
