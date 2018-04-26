/**
 * Declaration of ParticleOrientation.
 * @file ParticleOrientation.h
 * @author Ed Bachta, Jerry O. Talton III
 */

#ifndef __PARTICLEORIENTATION_H__
#define __PARTICLEORIENTATION_H__

#include "ParticleNormal.h"

/**
 * ParticleOrientation is a ParticleAttribute which represents
 * the orientation of particles in a system with a vector of
 * quaternions.
 */
class ParticleOrientation : public ParticleNormal {
public:

	MAKE_PARTICLESTUFF_NAME();

	/// Set the normal of a particle, and thus its orientation
	void setNormal(int i, gmVector3 normal);

	/// Get the normal of a particle
	gmVector3 getNormal(int i);
	
	/// Set the orientation of a particle
	void setOrientation(int i, gmQuat q);

	/// Get the orientation of a particle
	gmQuat getOrientation(int i);

	/// Get the orientation of a particle as a rotation matrix
	inline gmMatrix4 getMatrix(int i)
	{
		return quat_to_matrix(quat[i]);
	}
	
	
	/// Get a coordinate axis of the orientation
	gmVector3 getAxis(int i, int axis);

	/// Add particle orientation to a system of particles.
	ParticleOrientation(Particles *ps=NULL, const std::string& name=std::string("ParticleOrientation"));

	virtual void setParticleSystem(Particles *new_ps);

	void clear();

	/// Callback for particle addition.
	virtual void particleAdded();

	/// Callback for particle removal.
	virtual void particleRemoved(unsigned int i);

protected:
	/// The orientation of the particles as quaternions
	std::vector<gmQuat> quat;
};

#endif
