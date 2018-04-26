/*
@file Rainbow.h
@author John C. Hart
*/

#ifndef RAINBOW_H
#define RAINBOW_H

#include <vector>

#include "Particles.h"
#include "ParticleShader.h"

class ParticleMaterial;

class Rainbow : public ParticleShader
{
public:
	MAKE_PARTICLESTUFF_NAME();

	std::string source;		///< variable to index into the rainbow map
	double red;				///< value to map low end
	double violet;			///< value to map high end
	bool wrap;				///< whether or not to cycle color map

	ParticleMaterial *material;	///< target material onto which to store the color

	std::vector<double> *source_data;	///< pointer to per-particle source data to index color map

	Rainbow(Particles *ps=NULL);
	virtual void attachAttributes();

	virtual void drawParticle(int i);
};

#endif
