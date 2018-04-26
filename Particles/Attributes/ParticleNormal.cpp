/**
 * Implementation of ParticleNormal.
 * @file ParticleNormal.cpp
 * @author Ed Bachta, Jerry O. Talton III, John C. Hart
 */

#include "ParticleNormal.h"
#ifdef _WIN32
	#include <float.h>
#else
	#include <math.h>
	extern "C" int isnan(double);
	#define _isnan isnan
#endif

REGISTER_PARTICLESTUFF(ParticleNormal,"Attribute:ParticleNormal");

/**
 * Add particle orientation to a system of particles.
 * If particles exist in the system, default their orientation
 * to the identity quaternion.
 * @param ps   The owning particle system.
 * @param name The name of this object.
 */
ParticleNormal::ParticleNormal(Particles *ps,const std::string& name, bool nonormal) 
  : ParticleAttribute(ps, name)
{
	new PSParamgmVector3(this,&gravity,gmVector3(0.0,1.0,0.0),
		"gravity","gravity","Vector used as an up vector to align ambiguous orientation.");

	if (!nonormal) {
		new PSParamgmVector3PerParticle(this,&n,"n","normal","Particle normal.");
	}

	setParticleSystem(ps);
}

void ParticleNormal::setParticleSystem(Particles *new_ps)
{
	ParticleAttribute::setParticleSystem(new_ps);

	if (ps)
	{
		gradMag2.resize(ps->size());
		n.resize(ps->size());
		for (unsigned int i=0;i<ps->size();i++) {
			n[i] = gmVector3();
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
void ParticleNormal::setNormal(int i, gmVector3 normal)
{
	n[i] = normal;
}

/**
 * Get the normal of a particle.  The default normal is (0, 0, 1).
 * 
 * @param i   Index of particle.
 */
gmVector3 ParticleNormal::getNormal(int i)
{
	return n[i];
}

/**
 * Set the gravity.
 * 
 * @param grav   The gravity.
 */
void ParticleNormal::setGravity(gmVector3 grav)
{
	gravity = grav;
}

/**
 * Get the gravity.
 */
gmVector3 ParticleNormal::getGravity()
{
	return gravity;
}

/**
 * Set the squared magnitude of the gradient of a particle.
 * 
 * @param i Index of particle.
 * @param mag   The magnitude.
 */
void ParticleNormal::setGradMag2(int i, double mag)
{
	gradMag2[i] = mag;
}

/**
 * Get the squared magnitude of the gradient of a particle.
 * 
 * @param i Index of particle.
 */
double ParticleNormal::getGradMag2(int i)
{
	return gradMag2[i];
}

/** Returns an of a number of rotation matrices that rotates the z-axis to
 * the normal direction n[i]. (One degree of freedom is the twist about the normal)
 */
gmMatrix4 ParticleNormal::getMatrix(int i)
{
	gmVector3 nhat = n[i];
	nhat.normalize();
	gmVector3 v = cross(gmVector3(0.0,0.0,1.0), nhat);
	// if n is close to z, then no rotation needed, return identity
	// but may be 180 degrees off
	if (gmIsZero(v.lengthSquared())) {
		gmVector3 xaxis(1.0,0.0,0.0);
		return nhat[2] > 0.0 ? gmMatrix4::identity() : gmMatrix4::rotate(180.0,xaxis);
	}
	v.normalize();
	gmVector3 u = cross(nhat,v);
	// u should be normalized since v and n are unit and perp.

	// rotation matrix is [v u n]
	return gmMatrix4(v[0], u[0], nhat[0], 0.0,
					 v[1], u[1], nhat[1], 0.0,
					 v[2], u[2], nhat[2], 0.0,
					  0.0,  0.0,     0.0, 1.0);
}

void ParticleNormal::clear()
{
	n.clear();
	gradMag2.clear();
}

/**
 * Add a quaternion corresponding to the new particle.
 * @param i Index of the new particle.
 */
void ParticleNormal::particleAdded()
{
	n.push_back(gmVector3());
	gradMag2.push_back(0.0);
}

/**
 * Callback for particle removal.
 * @param i Index of particle to be removed.
 * @see Particles::particleRemoved
 */
void ParticleNormal::particleRemoved(unsigned int i)
{
	n[i] = n.back();
	n.pop_back();
	gradMag2[i] = gradMag2.back();
	gradMag2.pop_back();
}
