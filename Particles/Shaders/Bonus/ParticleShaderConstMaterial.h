/*
@file ParticleShaderConstMaterial.h
@author John C. Hart
Sets all particles to the same material properties
*/

#ifndef PARTICLESHADERCONSTMATERIAL_H
#define PARTICLESHADERCONSTMATERIAL_H

#include "ParticleShader.h"

class ParticleMaterial;

class ParticleShaderConstMaterial: public ParticleShader
{

protected:
	gmVector4 diffuseFront;
	gmVector4 diffuseBack;
	gmVector4 color;

	ParticleMaterial *material;

public:

	MAKE_PARTICLESTUFF_NAME();

	ParticleShaderConstMaterial(Particles *ps=NULL);

	virtual void drawShape(int i);
};

#endif
