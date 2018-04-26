/**
* Declaration of a curve adhesion behavior
 * @file CurveAdhesion.h
 * @date 9/27/2005
 * @author Matei N. Stroila
 * @remarks
 */


#ifndef CURVEADHESION_H
#define CURVEADHESION_H

#include "ParticleBehavior.h"

class ParticleNormal;
class ImplicitInterrogator;
class ParticlePosition;
class ParticleVelocity;

class CurveAdhesion : public ParticleBehavior {
	
private:
	
	/// Interrogation surface.
	ImplicitInterrogator *imp_int1, *imp_int2;
	
	/// Particle orientations.
	ParticleNormal* p_orient;
	
	ParticlePosition* position;
	ParticleVelocity* velocity;
	
public:
		
		MAKE_PARTICLESTUFF_NAME();
	
	/** Feedback constant
		*/
	double phi;
	
	/** Creates a curve drive attribute for Particles p.
		*/
	CurveAdhesion(Particles *ps=NULL);
	
	/** Apply the curve drive constraint.
		*/
	void applyConstraint();
	
	virtual void particleAdded();
	
	void attachAttributes();
	
	
};

#endif
