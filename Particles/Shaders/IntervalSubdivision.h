/*
@file IntervalSubdivision.h
@author John C. Hart
@date 7 Dec. 2005
*/

#ifndef INTERVALSUBDIVISION_H
#define INTERVALSUBDIVISION_H

#include "ParticleShader.h"
#include "Surface/Box.h"
class ImplicitInterrogator;
class ParticleBoundingBox;
class Implicit;

class IntervalSubdivision: public ParticleShader
{
public:
	MAKE_PARTICLESTUFF_NAME();

	double res;
	gmVector3 color;
	ImplicitInterrogator *imp_int_A, *imp_int_B, *imp_int_C, *imp_int_D;
	ParticleBoundingBox *pbox;

	Implicit *impA, *impB, *impC, *impD;

	IntervalSubdivision(Particles *ps=NULL);
   
	virtual void drawPre();
	virtual void subdivide(Box3d b);
};

#endif
