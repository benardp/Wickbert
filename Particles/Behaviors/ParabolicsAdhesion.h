/*
 *  ParabolicsAdhesion.h
 *  GPS
 *
 *  Created by Matei Stroila on 12/9/05.
 * 
 *
 */

#ifndef PARABOLICSADHESION_H
#define PARABOLICSADHESION_H

#include "ParticleBehavior.h"

class ParticleNormal;
class ImplicitInterrogator;
class ParticlePosition;
class ParticleVelocity;
class Implicit;

class ParabolicsAdhesion : public ParticleBehavior {
	
private:
	
	/// Interrogation surface.
	ImplicitInterrogator *imp_int;
	
	/// Particle orientations.
	ParticleNormal* p_orient;
	
	ParticlePosition* position;
	ParticleVelocity* velocity;
	
	Implicit *imp1, *imp2;
	
	bool ImplicitIsSet;
	
	bool setImplicit();
	
public:
		
		MAKE_PARTICLESTUFF_NAME();
	
	/** Feedback constant
		*/
	double phi;
	
	/** Creates a prabolic contour drive attribute for Particles p.
		*/
	ParabolicsAdhesion(Particles *ps=NULL);
	
	~ParabolicsAdhesion();
	
	/** Apply the curve drive constraint.
		*/
	void applyConstraint();
	
	virtual void particleAdded();
	
	void attachAttributes();
	
	
};

#endif
