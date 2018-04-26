/*
@file ParticleShaderOriented.h
@author Jerry O. Talton III
This class draws an oriented disk for a particle object.
*/

#ifndef PARTICLESHADERORIENTED_H
#define PARTICLESHADERORIENTED_H

#include "ParticleShader.h"
#include "ImplicitInterrogator.h"

class ParticleShaderOriented: public ParticleShader
{

protected:
	GLUquadricObj* quad;
	double radius;
	double scale;
	int sides;

	// this gives information for shader to draw the size of objects
	AdaptiveRepulsionData *ARData;

	ImplicitInterrogator *imp_int;

public:

	MAKE_PARTICLESTUFF_NAME();

	ParticleShaderOriented(Particles *ps=NULL);
	virtual void attachAttributes();
   
	virtual ~ParticleShaderOriented();

	virtual void drawShape(int i);


	/// parameters
	int qlen();
	void getq(double *q);
	void setq(double *q);
	void qname(char **qn);

};

#endif
