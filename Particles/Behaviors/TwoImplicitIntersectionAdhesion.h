/**
 * Declaration of the specialized curve adhesion class
 * @file TwoImplicitIntersectionAdhesion.h
 * @date 08/13/2006
 * @author Elmar Eisemann
 * @remarks The name of this class might not be perfect... 
 * Specular adhesion, silhouette adhesion, shadow adhesion and so on share one important
 * property. The resulting curves can be described as the intersection of two manifolds. An implicit surface and 
 * its (depending on the curve type), associated manifold.
 * Therefore this base class groups these kind of curves together.
 * This allows their usage in the code without explicit knowledge about the types.
 * In particular the base clase will provide a function definingManifoldProc(pos) that can be used to interrogate the manifold
 * describing the type of the surface, the implicit surface was already accessible through the implicit interrogator. 
 * The second manifold should be described in a way that when the curve bounds a special area (for example: shadow, spec)
 * its value is positive.
 * This allows to test for a point on the implicit whether it belongs to the described area.
 * for outside the implicit this simply results in a "proc" value of the second manifold, thus the name.
 * This also hides all "hacks", as some classes only work with directional light (and thus convert), others only for orthogonal view and so on...
 * So far we tested in the SVG and needed to reproduce these "hacks", now at least they are encapsulated in the classes and it is the author's responsability
 * to provide a meaningful result.
 * ATTENTION: An area might belong to several classes! (for example a specularity is also on the surface described by the silhouette)
 * This should not be seen as a disadvantage. The user can test for each point the corresponding types and then decide 
 * on his own about the shading, for example multiple light sources are possible like this.
 * ATTENTION: There might be some curves that do not bound an interesting area (for example for suggestive contours, 
 * negative and positive radial curvature are arguably both of the same interest)
 */
#ifndef TWOIMPLICITINTERSECTIONADHESION_H
#define TWOIMPLICITINTERSECTIONADHESION_H

#include "ParticleBehavior.h"

class ParticleAge;

class TwoImplicitIntersectionAdhesion : public ParticleBehavior {
		
public:
	//This is a pure virual baseclass, we do not want to register it!
	// Sure we do. Not all intersection curves are predefined, like silhouettes. -jch 9 Jan 07
	//MAKE_PARTICLESTUFF_NAME();
	
	/** Creates a silhouette drive attribute for Particles p.
		*/
	TwoImplicitIntersectionAdhesion(Particles *ps, const std::string & name);
	//returns a double that indicates the proc value of the defining second manifold (not the implicit surface of the interrogator).
	//It is in the responsability of the author to assure, that if an area is meaningful, to orient this function to be NEGATIVE
	//(thus inside the second manifold) in this bounded area on the base implicit.
	virtual double definingManifoldProc(const gmVector3 & position) const=0;

	void attachAttributes();

protected:

	ParticleAge *pAge;
	double pAgeThreshold;

};
#endif
