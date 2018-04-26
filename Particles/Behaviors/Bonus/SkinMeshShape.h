/*
* @file SkinMeshShape.h
* @date 2003-09-20
* @author Wen Su
* This implements Lee Markosian's Siggraph 99. Skin: A constructive approach to modeling free-form shapes.
*/

#ifndef SKINMESHSHAPE_H
#define SKINMESHSHAPE_H

#include "Particles.h"
#include "ParticleBehavior.h"
#include "ParticleMesh.h"
class ImplicitInterrogator;
class ParticlePosition;
class ParticleVelocity;

class SkinMeshShape : public ParticleBehavior
{
public:

	double targetLength;
	int allowSplit;
	int allowCollapse;
	bool allowFlip;
	ImplicitInterrogator *impInt;

    ParticlePosition* position;
	ParticleVelocity* velocity;

	MAKE_PARTICLESTUFF_NAME();

	void attachAttributes();

	void edgeFlip(ParticleOpenMesh &mesh);

	void edgeCollapseAndSplit(ParticleOpenMesh &mesh);

	/// default contructor
	SkinMeshShape(Particles *ps=NULL);

	/// keep triangles equal edges length
	virtual void applyForce();

	/// edgeflips and subdivisions
	virtual void cleanup();

	// calculate targetEdgeLength;
	void calculateTargetEdgeLength();

	// find the minimum of the interrior angles of two triangles abc, bcd
	double findMinAngle(const ParticleOpenMesh::Point &a, const ParticleOpenMesh::Point &b,
		const ParticleOpenMesh::Point &c, const ParticleOpenMesh::Point &d);

	// find the minimum of the interrior angles of triangles abc
	double findMinAngle(const ParticleOpenMesh::Point &a, const ParticleOpenMesh::Point &b,
		const ParticleOpenMesh::Point &c);

	// find the angle between to vectors with dot product
	double angleOf(const ParticleOpenMesh::Point &a, const ParticleOpenMesh::Point &b);
};

#endif
