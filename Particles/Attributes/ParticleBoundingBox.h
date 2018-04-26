//@file ParticleBoundingBox.h
//@date 2004-02-19
//@author Wen Su
//This used to be a behavior, but it really should be an attribute, the substitude behavior is KeepInBounds.
// By default the box {(-10,-10,-10), (10,10,10)} is created.


#ifndef PARTICLEBOUNDINGBOX_H
#define PARTICLEBOUNDINGBOX_H

#include "ParticleAttribute.h"

class ParticlePosition;

class ParticleBoundingBox : public ParticleAttribute
{
public:
	MAKE_PARTICLESTUFF_NAME();

	gmVector3 min; ///< Minimum coordinate.
	gmVector3 max; ///< Maximum coordinate.

	/// Create a bounding box.
	ParticleBoundingBox(Particles *ps=NULL, const std::string& name=std::string("ParticleBoundingBox"));

	/// check if a point is inside
	bool inBounds(const gmVector3 &p) const;

	void tighten();

	/// find the current bounding box (find the minimum and maximum particles' coordinates, and store them in _min and _max) 
	void computeBoundingBox(gmVector3& _min, gmVector3& _max) const;

	/**
	* Compute the intersection @param intersectBbox of the view vector with the bounding box
	*/
	void computeIntersection(	const int i, 
								const gmVector3& myCameraPosition,
								const double bboxOffset,
								gmVector3& intersectBbox ) const;
	  /**
	* Compute the two intersections @param v1 and @param v2 
	* of the view vector passing through @param P with the bounding box
	*/
	void computeIntersection(	const gmVector3& P, 
								const gmVector3& myCameraPosition, 
								const double bboxOffset, 
								gmVector3& v1, gmVector3& v2 ) const;
	
private:
		ParticlePosition *position;
};

class ParticleBoundingBoxTighten : public PSParamButton::Callback
{
public:
	ParticleBoundingBox *pbb;
	ParticleBoundingBoxTighten(ParticleBoundingBox *me) {pbb = me;}
	virtual void onbuttonpress() {pbb->tighten();}
};

#endif
