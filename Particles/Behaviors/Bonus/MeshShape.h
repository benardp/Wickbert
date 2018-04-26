/*
* @file MeshShape.h
* @date 2003-09-20
* @author Wen Su
* This will be used to maintain triangle shape, take the edge length and grow or shrink edge lengths.
*/

#ifndef MESHSHAPE_H
#define MESHSHAPE_H

#include "Particles.h"
#include "ParticleBehavior.h"
class ParticleMesh;
class ImplicitInterrogator;

class MeshShape : public ParticleBehavior
{
public:
	double targetTriangleArea;
	double edgeLengthScaler;
	double weightTolerance;
	double shortestWeightTolerance;
	double longestWeightTolerance;
	double splitStepLimit;
	double collapseStepLimit;

	MAKE_PARTICLESTUFF_NAME();

	/// default contructor
	MeshShape(Particles *ps=NULL);

	/// keep triangles equal edges length
	virtual void applyForce();

	/// edgeflips
	virtual void applyConstraint();

	/// subdivisions
	virtual void cleanup();

	int qlen();
	void getq(double *q);
	void setq(double *q);
	void qname(char **qn);

};

#endif
