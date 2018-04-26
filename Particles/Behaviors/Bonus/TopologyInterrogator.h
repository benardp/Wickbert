/**
* @file TopologyInterrogator.h
* @date 2003-09-07
* @author Wen Su

Implicit Mesh needs
Attributes:
	ParticleMesh
	ImplicitInterrogator
Behaviors:
	SurfaceAdhesion
	ParticleRepulsion
	TopologyInterrogator
Shaders
	ParticleShaderTriangle
*/

#ifndef TOPOLOGYINTERROGATOR_H
#define TOPOLOGYINTERROGATOR_H

#include "Particles.h"
#include "ParticleBehavior.h"
class ParticleMesh;
class ImplicitInterrogator;

class TopologyInterrogator : public ParticleBehavior
{
public:
	MAKE_PARTICLESTUFF_NAME();

	/// Reference to a locality object.
	ImplicitInterrogator *im;

	/// subdivide control
	double subdivideFlag;

	// this controls subdivision
	double reStart;

	// start up size of the sphere
	double startUpSize;

	/// default contructor
	TopologyInterrogator(Particles *ps=NULL);

	/// find ImplicitInterrogator and ParticleMesh
	virtual void attachAttributes();

	/// clean up calls updates the position
	virtual void cleanup();

	/// keep triangles equal edges length
	virtual void applyForce();

	virtual void applyConstraint();

	int qlen();
	void getq(double *q);
	void setq(double *q);
	void qname(char **qn);

};

#endif
