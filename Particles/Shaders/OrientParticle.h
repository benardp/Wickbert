/*
@file OrientParticle.h
@author John C. Hart
Applies the orientation stored in Attribute:ParticleOrientation to the particles.
Rotates the OpenGL modelview coordinate system to the appropriate orientation
*/

#ifndef ORIENTPARTICLE_H
#define ORIENTPARTICLE_H

#include "ParticleShader.h"

class ParticleNormal;
class ParticlePosition;

class OrientParticle : public ParticleShader
{

protected:
	ParticlePosition *position;
	ParticleNormal *orientation;

public:

	MAKE_PARTICLESTUFF_NAME();

	OrientParticle(Particles *ps=NULL);

	virtual void drawParticle(int i);
};

#endif
