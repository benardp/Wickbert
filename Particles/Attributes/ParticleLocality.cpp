/**
 * Implementation of ParticleLocality.
 * @file ParticleLocality.cpp
 * @author Ed Bachta
 */

#include "ParticleLocality.h"
#include "ParticlePosition.h"

REGISTER_PARTICLESTUFF(ParticleLocality,"Attribute:ParticleLocality");

/**
 * Constructor.
 * @param ps   The owning particle system.
 * @param name The name of this locality object.
 */ 
ParticleLocality::ParticleLocality(Particles *ps, const std::string& name)
  : ParticleAttribute(ps, name)
{
	new Attached<ParticlePosition>(this,&position);
	new PSParamBool(this,&useCachedNeighbors,false,"useCachedNeighbors","Use Cache Neighbors","Set true to use cached neighbors.");
}

/**
 * Return possible neighbors of a particle.
 * @param i    Index of the particle in question.
 * @param r				Radius of query sphere
 * @param neighbors		A vector of indices of
 *						particles within query sphere
 *
 */
void ParticleLocality::getNeighbors(const unsigned int i, const double queryRadius, std::list<unsigned int> &neighbors)
{
	if(!useCachedNeighbors) findNeighbors(i,queryRadius);
	neighbors = neighborsMap[i];
}
/**
 * Find possible neighbors of a particle. They are stored in neighborsMap
 * @param i    Index of the particle in question.
 * @param r				Radius of query sphere
 */
void ParticleLocality::findNeighbors(const unsigned int i, const double queryRadius)
{
	gmVector3 x = position->getPosition(i);
	gmVector3 d;
	double qr2 = queryRadius * queryRadius;
	NeighborsT ithNeighbors;
	for (unsigned int j = 0; j < ps->size(); j++) {
		if (j == i) continue;
		d = position->getPosition(j) - x;
		if (d.lengthSquared() <= qr2)
			ithNeighbors.push_back(j);
	}
	neighborsMap[i] = ithNeighbors;
}

