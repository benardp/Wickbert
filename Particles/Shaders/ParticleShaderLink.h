/*
@file ParticleShaderLink.h
@author Wen Su
This class draws links with neightbors.
*/

#ifndef PARTICLESHADERLINK_H
#define PARTICLESHADERLINK_H

#include "ParticleShader.h"
#include "ParticleLocalityGrid.h"
class ParticleLocality;
class AdaptiveRepulsionData;

class ParticleShaderLink: public ParticleShader
{
public:
	MAKE_PARTICLESTUFF_NAME();

	ParticleShaderLink(Particles *ps=NULL);

	//ParticleLocality *p_locality;
	ParticleLocalityGrid *p_locality;
	AdaptiveRepulsionData *ardata;

	virtual void drawPre();

	/// Draw the links connecting the selected particle
	/// to the other particles that it repels.
	virtual void drawParticle(int i);
};

#endif
