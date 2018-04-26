/**
 * Declaration of the silhouette constraint.
 * @file SilhouetteAdhesion.h
 * @author Matei N. Stroila
 * @date 6/22/2005
 * @remarks 
 */


#ifndef SILHOUETTEADHESION_H
#define SILHOUETTEADHESION_H

#include "TwoImplicitIntersectionAdhesion.h"

class ImplicitInterrogator;
class ParticleNormal;
class ParticlePosition;
class ParticleVelocity;
class ViewDependence;

class SilhouetteAdhesion : public TwoImplicitIntersectionAdhesion {
	
private:
	void initParameters();
	double onSilhouetteSurfProc(const gmVector3& pos, const gmVector3& grad) const;

	/// Interrogation surface.
	ImplicitInterrogator* imp_int;
	
	/// Particle orientations.
	ParticleNormal* p_orient;
	
	ParticlePosition* position;
	ParticleVelocity* velocity;
	ViewDependence *view;

	gmVector3 cameraPos;
	
public:
		
		MAKE_PARTICLESTUFF_NAME();
	
	/** Feedback constant
		*/
	double phi;
	
	/** Creates a silhouette drive attribute for Particles p.
		*/
	SilhouetteAdhesion(Particles *ps=NULL);

	virtual double definingManifoldProc(const gmVector3 & position) const;



	/** Apply the silhouette drive constraint.
		*/
	void applyConstraint();
	
	virtual void particleAdded();
	
	void attachAttributes();
	
		
};

#endif
