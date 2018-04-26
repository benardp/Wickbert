/**
* Declaration of the shadow constraint.
 * @file ShadowAdhesion.h
 * @date June 06, 2005
 * @author Matei Stroila
 */

#ifndef SHADOWADHESION_H
#define SHADOWADHESION_H

#include "TwoImplicitIntersectionAdhesion.h"

class ImplicitInterrogator;
class ParticleNormal;
class ParticlePosition;
class ParticleVelocity;
class LightPosition;

class ShadowAdhesion : public TwoImplicitIntersectionAdhesion{
	
private:
	void initParameters();
	double onShadowSurfProc(const gmVector3& pos, const gmVector3& grad) const;

	/// Interrogation surface.
	ImplicitInterrogator* imp_int;
	/// Particle orientations.
	ParticleNormal* p_orient;
	ParticlePosition* position;
	ParticleVelocity* velocity;
	LightPosition *light; 
	gmVector3 lightPos3;
	bool directionalLight;
	
public:
		
		MAKE_PARTICLESTUFF_NAME();
	
	/** Feedback constant
		*/
	double phi;
	
	/** Creates a silhouette drive attribute for Particles p.
		*/
	ShadowAdhesion(Particles *ps=NULL);
	
	double definingManifoldProc(const gmVector3& pos) const;
	
	/** Apply the silhouette drive constraint.
		*/
	void applyConstraint();
	
	virtual void particleAdded();
	
	void attachAttributes();
	
	
};

#endif
