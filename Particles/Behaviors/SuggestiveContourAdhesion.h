/**
* Declaration of the suggestive contour adhesion behavior
* @file SuggestiveContourAdhesion.h
* @date 23 June. 2005
* @author Matei N. Stroila
*/
#ifndef SUGGESTIVECONTOURDRIVE_H
#define SUGGESTIVECONTOURDRIVE_H

#include "TwoImplicitIntersectionAdhesion.h"


class ImplicitInterrogator;
class ViewDependence;
class ParticleNormal;
class ParticlePosition;
class ParticleVelocity;
class ParticleCreation;

class SuggestiveContourAdhesion : public TwoImplicitIntersectionAdhesion {
	
public:
		
	MAKE_PARTICLESTUFF_NAME();
	
	/** Creates a silhouette drive attribute for Particles p.
		*/
	SuggestiveContourAdhesion(Particles *ps=NULL);
	virtual double definingManifoldProc(const gmVector3 & position) const;
		
	/** Apply the suggestive drive constraint.
		*/
	void applyConstraint();
	
	virtual void particleAdded();
	
	void attachAttributes();

private:

	double onSuggSurfProc(const gmVector3& pos, const gmVector3& Hessv) const;
	
	/** Feedback constant
		*/
	double phi;

	/// Interrogation surface.
	ImplicitInterrogator* imp_int;
	
	/// Particle orientations.
	ParticleNormal* p_orient;
	
	ParticlePosition* position;
	ParticleVelocity* velocity;
	ViewDependence *view;
	ParticleCreation* pCreation;
	
	gmVector3 cameraPos;
	void initParameters();

	//check if the suggestive contour is defined
	bool isDefined(gmVector3 dir, gmVector3 normal) const;
	
};

#endif
