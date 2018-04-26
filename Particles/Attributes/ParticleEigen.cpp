/**
 * @file ParticleEigen.cpp
 * @author Terry Fleury
 */

#include "ParticleEigen.h"

REGISTER_PARTICLESTUFF(ParticleEigen,"Attribute:ParticleEigen");

/**
 * Add particle eigen vectors/values to a system of particles.
 * @param ps   The owning particle system.
 * @param name The name of this object.
 */
ParticleEigen::ParticleEigen(Particles *ps, const std::string& name)
	: ParticleAttribute(ps, name) 
{
	if (ps)
    {
      int size = ps->size();
		  eigenVectors.resize(size);
		  eigenValues.resize(size);
		  positiveEigenValues.resize(size);
    }

	new PSParamgmVector3PerParticle(this,&eigenValues,
		"vals","eigenvalues","Magnitude of the principal directions");
	new PSParamgmMatrix3PerParticle(this,&eigenVectors,
		"vecs","eigenvectors","Principal directions of each particle");
	new PSParamIntPerParticle(this,&positiveEigenValues,
		"posevs","positive eigenvalues","Number of eigenvalues that are positive");
}

/**
 * Callback for addition of a particle.
 * @param i Index of the new particle.
 */
void ParticleEigen::particleAdded() 
{
	eigenVectors.push_back(gmMatrix3());
	eigenValues.push_back(gmVector3());
  positiveEigenValues.push_back(0);
}

/**
 * Callback for particle removal.
 * @param i Index of particle to be removed.
 */
void ParticleEigen::particleRemoved(unsigned int i) 
{
  eigenVectors[i] = eigenVectors.back();
	eigenVectors.pop_back();
  eigenValues[i] = eigenValues.back();
	eigenValues.pop_back();
  positiveEigenValues[i] = positiveEigenValues.back();
	positiveEigenValues.pop_back();
}

gmMatrix3 ParticleEigen::getEigenVectors(unsigned int i)
{
	return eigenVectors[i];
}

void ParticleEigen::setEigenVectors(unsigned int i, gmMatrix3 ev)
{
	eigenVectors[i] = ev;
}

gmVector3 ParticleEigen::getEigenValues(unsigned int i)
{
	return eigenValues[i];
}

void ParticleEigen::setEigenValues(unsigned int i, gmVector3 ev)
{
	eigenValues[i] = ev;
}

int ParticleEigen::getPositiveEigenValues(unsigned int i)
{
	return positiveEigenValues[i];
}

void ParticleEigen::setPositiveEigenValues(unsigned int i, int pev)
{
	positiveEigenValues[i] = pev;
}

void ParticleEigen::clear()
{
  eigenVectors.clear();
  eigenValues.clear();
  positiveEigenValues.clear();
}

