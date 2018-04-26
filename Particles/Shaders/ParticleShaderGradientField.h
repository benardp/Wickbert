/*
@file ParticleShaderGradientField.h
@author Wen Su
This class draws a gradient field.
*/

#ifndef PARTICLESHADERGRADIENTFIELD_H
#define PARTICLESHADERGRADIENTFIELD_H

class ParticleShader;

class ParticleShaderGradientField: public ParticleShader
{
protected:
	// find implicit surface
	ImplicitInterrogator *impInt;
	// bounds
	ParticleBoundingBox *bounds;
	// parameters
	double vectorLength;
	double stepSize;
	bool normalizeVector;

public:

	MAKE_PARTICLESTUFF_NAME();

	ParticleShaderGradientField(Particles *ps=NULL);
   
	virtual void draw(int s);
};

#endif
