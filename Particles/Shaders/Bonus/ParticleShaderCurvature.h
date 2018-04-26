/*
@file ParticleShaderCurvature.h
@author John C. Hart
Sets material properties based on curvature
*/

#ifndef PARTICLESHADERCURVATURE_H
#define PARTICLESHADERCURVATURE_H

#include "ParticleShader.h"
#include "ImplicitInterrogator.h"

class ParticleShaderCurvature: public ParticleShader
{

protected:

	ImplicitInterrogator *imp_int;
	int kind;
	double scale;

public:

	MAKE_PARTICLESTUFF_NAME();

	ParticleShaderCurvature(Particles *ps=NULL);
	virtual void attachAttributes();
   
	virtual void drawShape(int i);

	/// parameters
	int qlen();
	void getq(double *q);
	void setq(double *q);
	void qname(char **qn);

	char *tip() {return "Sets material color according to curvature.";}
	char *qtip(int);



};

#endif
