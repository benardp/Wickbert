/**
 * Declaration of ParticleAge.
 * @file ParticleAge.h
 * @author Doug Nachand
 */

#ifndef PARTICLEAGE_H
#define PARTICLEAGE_H

#include "Particles.h"
#include "ParticleAttribute.h"

/**
 * ParticleAge is a ParticleAttribute that holds (but does not update)
 * the age of particles.
 */
class ParticleAge : public ParticleAttribute {

public:

	MAKE_PARTICLESTUFF_NAME();

	/// Particle ages
	std::vector<int> t;
	int maxLifeTime;

	/// Add paritcle age to a system of particles.
	ParticleAge(Particles *ps=NULL, const std::string& name=std::string("ParticleAge"));

	void clear();

	void reset();

	void resetParticleAge(const int i, const bool b);

	/// Callback for particle addition.
	virtual void particleAdded();

	/// Callback for particle removal.
	virtual void particleRemoved(unsigned int i);

	virtual void prepare();

private:

	bool useParticleAge;
	std::vector<unsigned int> resetFlag;

};

#endif
