/*
@file ShowCell.h
@author John C. Hart
*/

#ifndef SHOWCELL_H
#define SHOWCELL_H

#include "ParticleShader.h"
class ParticleLocalityGrid;

class ShowCell : public ParticleShader
{
public:
	MAKE_PARTICLESTUFF_NAME();

	ParticleLocalityGrid *plg;

	ShowCell(Particles *ps=NULL);

	virtual void drawParticle(int i);
};

#endif
