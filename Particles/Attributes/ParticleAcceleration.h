/*
@file ParticleAcceleration.h
@auther Yuan Zhou
*/

#ifndef PARTICLEACCELERATION_H
#define PARTICLEACCELERATION_H

#include "Particles.h"
#include "ParticleAttribute.h"

class ParticleAcceleration : public ParticleAttribute {

public:

	MAKE_PARTICLESTUFF_NAME();

	/// ParticleAcceleration
	std::vector<gmVector3> acc;

	/// Add paritcle orientation to a system of particles.
	ParticleAcceleration(Particles *ps=NULL, const std::string& name=std::string("ParticleAcceleration"));

	/// Callback for particle addition.
	virtual void particleAdded();

	/// Callback for particle removal.
	virtual void particleRemoved(unsigned int i);

	/// this gives a calling convertion for children to get the acceleration of index i
	virtual gmVector3 getAcceleration(unsigned int i);
	
	/// this gives a calling convertion for children to set the acceleration of index i
	virtual void setAcceleration(unsigned int i, gmVector3 p);

	virtual void attachAttributes() {
		if (ps && acc.size() != ps->size())
			acc.resize(ps->size());
	}
};

#endif
