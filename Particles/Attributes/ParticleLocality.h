/**
 * Declaration of ParticleLocality.
 * @file ParticleLocality.h
 * @author Ed Bachta
 */

#ifndef __PARTICLELOCALITY_H__
#define __PARTICLELOCALITY_H__

#include "ParticleAttribute.h"
class ParticlePosition;

/**
 * ParticleLocality is the base class for objects that maintain a notion
 * of the spatial relationship between particles. Such classes can optimize
 * the finding of particle neighbors.
 */
class ParticleLocality : public ParticleAttribute {

private:

	typedef std::list<unsigned int> NeighborsT;
	//Cache neighbors in a map
	std::map<unsigned int,NeighborsT> neighborsMap;
	/**
	* Find possible neighbors of a particle. They are stored in neighborsMap
	* @param i    Index of the particle in question.
	* @param r				Radius of query sphere
	*/
	void findNeighbors(const unsigned int i, const double queryRadius);
	//flag to use or not the cached neighbors
	bool useCachedNeighbors;

public:

	MAKE_PARTICLESTUFF_NAME();

	ParticlePosition *position;

	ParticleLocality(Particles* ps=NULL, const std::string& name = std::string("ParticleLocality"));

	virtual void update() { }

	/** Finds possible neighbors of a particle.
	 * Might use the internally stored query radius to limit
	 * the search to nearby particles. Use setQueryRadius to set.
	 */
	virtual void getNeighbors(const unsigned int i, const double queryRadius, std::list<unsigned int> &neighbors);

};



#endif
