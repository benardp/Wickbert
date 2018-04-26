/*
@file CopyParticle.h
@author John C. Hart
*/

#ifndef COPYPARTICLE_H
#define COPYPARTICLE_H

#include "ParticleShader.h"
#include "ParticlePosition.h"

/** This class is designed to copy the currently selected particle to
 * probably a different collection of Particles in the system.
 */

class CopyParticle : public ParticleShader
{
public:

	std::string targetparticles;
	ParticlePosition *pos;
	std::string dstPosName;
	int eventcode;

	MAKE_PARTICLESTUFF_NAME();

	CopyParticle(Particles *ps=NULL);
   
	virtual void event(int e);
};

#endif
