/**
Declaration of ParticleShaderCluster
@file ParticleShaderCluster.h
@author Wen Su
This class treate each particle as a cluster but also draws all the triangles
*/

#ifndef PARTICLESHADERCLUSTER_H
#define PARTICLESHADERCLUSTER_H

#include <vector>
#include "ParticleShader.h"

class Particles;
class SurfaceInterrogator;

class ParticleShaderCluster : public ParticleShader
{

public:
	double lineWidth;
	int useMaterial;
	int drawOutline;

	MAKE_PARTICLESTUFF_NAME();

	// default constructor
	ParticleShaderCluster(Particles *ps=NULL);

	virtual void draw(int s);

	SurfaceInterrogator *surInt;

	/// parameters
	int qlen();
	void getq(double *q);
	void setq(double *q);
	void qname(char **qn);

	void attachAttributes();
};

#endif