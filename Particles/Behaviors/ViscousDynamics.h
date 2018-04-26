/**
 * @file ViscousDynamics.h (originally VelocityIntegrator) 
 * @author by Wen Su
 */

#ifndef VISCOUSDYNAMICS_H
#define VISCOUSDYNAMICS_H

#include "Particles.h"
#include "ParticleBehavior.h"
class ParticlePosition;
class ParticleVelocity;
class ParticleScalar;

class ViscousDynamics : public ParticleBehavior
{
public:
	MAKE_PARTICLESTUFF_NAME();

	double dt;

	ParticlePosition *position;
	ParticleVelocity *velocity, *history;
//	ParticleScalar *speed;

	ViscousDynamics(Particles *ps=NULL, const std::string& name=std::string("ViscousDynamics"));

	void integrate();
};

#endif


