/**
 * @file VelocityIntegrator.h
 * @author Wen Su
 */

#ifndef VELOCITYINTEGRATOR_H
#define VELOCITYINTEGRATOR_H

#include "Particles.h"
#include "ParticleBehavior.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"

class VelocityIntegrator : public ParticleBehavior
{
public:
	MAKE_PARTICLESTUFF_NAME();

	VelocityIntegrator(Particles *ps=NULL, const std::string& name=std::string("VelocityIntegrator"));

	void integrate();
};

#endif


