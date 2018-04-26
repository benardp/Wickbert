/**
 * Declaration of the specular highlights adhesion behavior 
 * @file SpecularAdhesion.h
 * @date 06/26/2005
 * @author Matei N. Stroila
 * @remarks
 */

#ifndef SPECULARADHESION_H
#define SPECULARADHESION_H

#include "TwoImplicitIntersectionAdhesion.h"

class ImplicitInterrogator;
class LightPosition;
class ViewDependence;
class ParticleNormal;
class ParticlePosition;
class ParticleVelocity;

class SpecularAdhesion : public TwoImplicitIntersectionAdhesion{
		
public:		
	MAKE_PARTICLESTUFF_NAME();
	
	/** Creates a silhouette drive attribute for Particles p.
		*/
	SpecularAdhesion(Particles *ps=NULL);
	virtual double definingManifoldProc(const gmVector3 & position) const;

	/** Apply the silhouette drive constraint.
		*/
	void applyConstraint();
	virtual void particleAdded();
	void attachAttributes();

private:
	double onSpecSurfProc(const gmVector3 & gradient, double gradLength) const;
		
	/// Particle orientations.
	ParticleNormal* p_orient;

	/** Feedback constant and shineness
		*/
	double phi, shine;
	
	/// Interrogation surface.
	ImplicitInterrogator* imp_int;
	ParticlePosition* position;
	ParticleVelocity* velocity;
	LightPosition *light; 
	ViewDependence *view;
	
	void initParameters();
	gmVector3 halfVector;
	double halfVectorLength;
};
#endif
