/**
 * Declaration of ParticleNormal.
 * @file ParticleNormal.h
 * @author Ed Bachta, Jerry O. Talton III, John C. Hart
 */

#ifndef __PARTICLENORMAL_H__
#define __PARTICLENORMAL_H__

#include "ParticleAttribute.h"

/**
 * ParticleNormal is a ParticleAttribute which represents
 * the normal direction of particles in a system with a vector of 3-vectors.
 */
class ParticleNormal : public ParticleAttribute {
public:

	MAKE_PARTICLESTUFF_NAME();

	/// Add particle normal to a system of particles.
	ParticleNormal(Particles *ps=NULL, const std::string& name=std::string("ParticleNormal"), bool nonormal=false);

	virtual void setParticleSystem(Particles *new_ps);

	/// Set the normal of a particle, and thus its orientation
	virtual void setNormal(int i, gmVector3 normal);

	/// Get the normal of a particle
	virtual gmVector3 getNormal(int i);

	/// Set the gravity
	virtual void setGravity(gmVector3 grav);

	/// Get the gravity
	virtual gmVector3 getGravity();
	
	/// Set the squared magnitude of the gradient of a particle
	virtual void setGradMag2(int i, double mag);
	
	/// Get the squared magnitude of the gradient of a particle
	virtual double getGradMag2(int i);

	/// Get the orientation of a particle as a rotation matrix
	virtual gmMatrix4 getMatrix(int i);

	virtual void clear();

	/// Callback for particle addition.
	virtual void particleAdded();

	/// Callback for particle removal.
	virtual void particleRemoved(unsigned int i);

protected:
	gmVector3 gravity;

	/// The orientation of the particles as quaternions
	std::vector<gmVector3> n;

	/// Magnitude of the gradients squared
	std::vector<double> gradMag2;
};

#endif
