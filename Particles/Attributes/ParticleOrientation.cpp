/**
 * Implementation of ParticleOrientation.
 * @file ParticleOrientation.cpp
 * @author Ed Bachta, Jerry O. Talton III
 */

#include "ParticleOrientation.h"
#ifdef _WIN32
	#include <float.h>
#else
	#include <math.h>
	extern "C" int isnan(double);
	#define _isnan isnan
#endif

REGISTER_PARTICLESTUFF(ParticleOrientation,"Attribute:ParticleOrientation");

/**
 * Add particle orientation to a system of particles.
 * If particles exist in the system, default their orientation
 * to the identity quaternion.
 * @param ps   The owning particle system.
 * @param name The name of this object.
 */
ParticleOrientation::ParticleOrientation(Particles *ps,const std::string& name) 
  : ParticleNormal(ps, name, true)
{
}

void ParticleOrientation::setParticleSystem(Particles *new_ps)
{
	ParticleAttribute::setParticleSystem(new_ps);

	if (ps)
	{
		gradMag2.resize(ps->size());
		quat.resize(ps->size());
		for (unsigned int i=0;i<ps->size();i++) {
			quat[i] = gmQuat();
			gradMag2[i] = 0.0;
		}
	}
}

/**
 * Set the normal of a particle.
 * 
 * @param i   Index of particle.
 * @param normal The normal to use to establish the orientation.
 */
void ParticleOrientation::setNormal(int i, gmVector3 normal)
{
	gmVector3 nrm = normal.normalize();
	gmVector3 currNrm = getNormal(i);

	// Only change the normal if it needs to be changed 
	// (This should help prevent the dreaded disappearing particle bugs)
	//if (true || !gmIsZero(dot(nrm, currNrm)))
		quat[i] = vec_to_vec_quat(currNrm, nrm) * quat[i];

	// Fix dreaded disappearing particle bug by resetting particles that get NAN orientations
	if (_isnan((quat[i].vector())[0]) || _isnan((quat[i].vector())[1]) || _isnan((quat[i].vector())[2])
		|| _isnan(quat[i].scalar()))
		quat[i] = gmQuat();
}

/**
 * Get the normal of a particle.  The default normal is (0, 0, 1).
 * 
 * @param i   Index of particle.
 */
gmVector3 ParticleOrientation::getNormal(int i)
{
	return gmIsZero(quat[i].norm()) ? gmVector3(0.0, 0.0, 1.0) : rotate_by_quat(gmVector3(0.0, 0.0, 1.0), quat[i]);
}

/**
 * Set the orientation of a particle.
 * 
 * @param i Index of particle.
 * @param q The new orientation.
 */
void ParticleOrientation::setOrientation(int i, gmQuat q)
{
	quat[i] = q;
}

/**
 * Get the orientation of a particle.
 * 
 * @param i Index of particle.
 */
gmQuat ParticleOrientation::getOrientation(int i)
{
	return quat[i];
}

/**
 * Get the axis of a particle's local coordinate frame.
 * 
 * @param i Index of particle.
 * @param axis The axis to return (1 for x, 2 for y).  Use getNormal for z.
 */
gmVector3 ParticleOrientation::getAxis(int i, int axis)
{
	// Make sure the axis parameter is valid
	assert((axis == 1) || (axis == 2));

	if (axis == 1) return rotate_by_quat(gmVector3(1.0, 0.0, 0.0), quat[i]);
	else return rotate_by_quat(gmVector3(0.0, 1.0, 0.0), quat[i]);
}

void ParticleOrientation::clear()
{
	quat.clear();
	gradMag2.clear();
}


/**
 * Add a quaternion corresponding to the new particle.
 * @param i Index of the new particle.
 */
void ParticleOrientation::particleAdded()
{
	quat.push_back(gmQuat());
	gradMag2.push_back(0.0);
}

/**
 * Callback for particle removal.
 * @param i Index of particle to be removed.
 * @see Particles::particleRemoved
 */
void ParticleOrientation::particleRemoved(unsigned int i)
{
	quat[i] = quat.back();
	quat.pop_back();
	gradMag2[i] = gradMag2.back();
	gradMag2.pop_back();
}
