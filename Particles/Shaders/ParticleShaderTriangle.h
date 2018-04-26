/**
Declaration of ParticleShaderTriangle
@file ParticleShaderTriangle.h
@author Wen Su
This class find the possible singular curves by connection neighboring singularities
*/

#ifndef PARTICLESHADERTRIANGLE_H
#define PARTICLESHADERTRIANGLE_H

#include <vector>
#include "ParticleShader.h"
#include "ParticleMesh.h"

class Particles;
class ParticleMesh;
class ParticleNormal;

class ParticleShaderTriangle : public ParticleShader
{
private:
	void drawPolygon(ParticleOpenMesh &mesh);
	void drawWireframe(ParticleOpenMesh &mesh);

public:
	double lineWidth;
	bool faceNormal;
	ParticleMesh *pm;

public:
	MAKE_PARTICLESTUFF_NAME();

	// default constructor
	ParticleShaderTriangle(Particles *ps=NULL);

	virtual void drawPre();

};

#endif