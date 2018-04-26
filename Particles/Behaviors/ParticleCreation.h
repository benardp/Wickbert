//@file ParticleCreation.h
//@author Wen Su
//This behavior creates new particles within a bounding volume.

#ifndef PARTICLECREATION_H
#define PARTICLECREATION_H

#include "ParticleBehavior.h"

class ParticleBoundingBox;
class ParticlePosition;

class ParticleCreation : public ParticleBehavior
{
private:
	ParticleBoundingBox *bounds;
	ParticlePosition *position;

public:
	bool createNewParticles;
	int limitSize;

	MAKE_PARTICLESTUFF_NAME();

	/// default constructor
	ParticleCreation(Particles *ps=NULL, const std::string& name=std::string("ParticleCreation"));

	/// Applies a force that repels nearby particles.
	void cleanup();
};

#endif
