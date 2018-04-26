/**
* @file NeighborAlignment.h
* Declaration of the surface constraint.
*/

#ifndef NEIGHBORALIGNMENT_H
#define NEIGHBORALIGNMENT_H

#include "Particles.h"
#include "ParticleBehavior.h"
#include "ParticleOrientation.h"
#include "ImplicitInterrogator.h"
#include "ParticlePosition.h"


/**
*
* Aligns a particle with its neighbors
*
*/
class NeighborAlignment : public ParticleBehavior {
	
private:
	
	/// Particle orientations.
	ParticleOrientation* p_orient;

	/// Particle positions.
	ParticlePosition *position;

	/// Implicit
	ImplicitInterrogator *imp;

	/// The axis to align
	int axis;

public:

	MAKE_PARTICLESTUFF_NAME();

	int qlen();
	void getq(double *q);
	void setq(double *q);
	void qname(char **qn);

	void attachAttributes();
	
	/** Creates a surface adhesion attribute for Particles p.
	*/
	NeighborAlignment(Particles *ps=NULL);

	void applyConstraint();
	
};

#endif
