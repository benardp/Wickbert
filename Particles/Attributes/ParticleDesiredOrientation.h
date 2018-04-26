/**
 * Declaration of ParticleDesiredOrientation.
 * @file ParticleDesiredOrientation.h
 * @author Mike Flavin
 */

#ifndef __PARTICLEDESIREDORIENTATION_H__
#define __PARTICLEDESIREDORIENTATION_H__

#include "Particles.h"
#include "ParticleAttribute.h"

/**
 * ParticleDesiredOrientation is a ParticleAttribute. It
 * represents a *desired* normal for the surface at a point.
 * Currently, this attribute is only meaningful for RBF control
 * points, in which normal constraints are possible.  In the
 * standard least-squares Witkin-Heckbert (Wickbert) scheme,
 * only control point positions are considered, and this attribute
 * is ignored.
 */
class ParticleDesiredOrientation : public ParticleAttribute {

 public:

	MAKE_PARTICLESTUFF_NAME();

	/// Desired normals
	std::vector<gmVector3> n;

	/// Magnitude of the gradients squared
	std::vector<double> gradMag2;

	/// Add desired orientation to a system of particles.
	ParticleDesiredOrientation(Particles *ps=NULL, const std::string& name=std::string("ParticleDesiredOrientation"));

	virtual void setParticleSystem(Particles *new_ps);

	void clear();

	/// Callback for particle addition.
	virtual void particleAdded();

	/// Callback for particle removal.
	virtual void particleRemoved(unsigned int i);

};

#endif
