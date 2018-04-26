//@file KeepInBounds.h
//@date 2004-02-19
//@author Wen Su
//This is a behavior that find the attribute ParticleBoundingBox and keep the particles inside this box.

#ifndef KEEPINBOUNDS_H
#define KEEPINBOUNDS_H

#include "ParticleBehavior.h"
class ParticlePosition;
class ParticleVelocity;
class ParticleBoundingBox;

class KeepInBounds: public ParticleBehavior
{
private:
	ParticlePosition *position;
	ParticleVelocity *velocity;
	ParticleBoundingBox *bbox;

public:
	MAKE_PARTICLESTUFF_NAME();

	/// Create a bounding box.
	KeepInBounds(Particles *ps=NULL, const std::string& name=std::string("KeepInBounds"));

	/// Constrains particles so that they do not leave the box.
	void applyConstraint();
};

#endif
