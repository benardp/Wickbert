/**
* @file SurfaceAdhesion.h
* Declaration of the surface constraint.
*/

#ifndef SURFACEADHESION_H
#define SURFACEADHESION_H

#include "ParticleBehavior.h"

class ParticlePosition;
class ParticleVelocity;
class ParticleNormal;
class ImplicitInterrogator;
class ParticleAge;

/**
* This is a surface constraint described by Witkin and Heckbert S94
* but had been around for a while before. The constraint removes any component
* of the accumulated force normal to the surface, leaving only
* tangential motion.
*
* Required Attributes: ImplicitInterrogator, ParticleNormal.
*
* The ImplicitInterrogator attribute provides the surface the particles
* adhere to. The SurfaceAdhesion behavior also orients the particles
* to be normal to the surface by setting the per-particle surface
* normals contained in the ParticleNormal attribute.
*
* SurfaceAdhesion contains a single parameter phi which is a feedback
* constant used for the penalty force that keeps the particle on
* the surface when numerical error in the dynamic constraint causes
* it to drift off the zero-surface.
*/
class SurfaceAdhesion : public ParticleBehavior {
	
private:
	
	/// Interrogation surface.
	ImplicitInterrogator* imp_int;
	
	/// Particle orientations.
	ParticleNormal* p_orient;

	ParticlePosition* position;
	ParticleVelocity* velocity;

	ParticleAge* pAge;
	
	double singTh; //singularity threshold
	double pAgeThreshold;

public:

	MAKE_PARTICLESTUFF_NAME();
	
	/** Feedback constant
	*/
	double phi;

	/** Creates a surface adhesion attribute for Particles p.
	*/
	SurfaceAdhesion(Particles *ps=NULL);
	
	/** Apply the surface adhesion constraint.
	*/
	void applyConstraint();
	
	virtual void particleAdded();
	
};

#endif
