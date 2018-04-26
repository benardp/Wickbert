/*
@file ParticleShaderText.h
@author Tony Kaap
*/

#ifndef PARTICLESHADERTEXT_H
#define PARTICLESHADERTEXT_H

#include "ParticleShader.h"
#include "ParticleScalar.h"

class ParticleShaderText:public ParticleShader
{
public:

	ParticleScalar * scalar;
	bool useScalar;
	//GLUquadricObj* quad;
	//double radius;
	//int sides;

public:

	MAKE_PARTICLESTUFF_NAME();
	
	ParticleShaderText(Particles *ps=NULL);

	virtual ~ParticleShaderText();

	virtual void drawParticle(int i);

};

#endif


