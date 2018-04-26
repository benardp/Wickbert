#ifndef PARTICLEDENSITY_H
#define PARTICLEDENSITY_H

#include "../Particles.h"
#include "../ParticleAttribute.h"

class ParticleDensity : public ParticleAttribute {

public:

	MAKE_PARTICLESTUFF_NAME();

	/// ParticleDensity
	std::vector<double> den; 

	/// Add paritcle orientation to a system of particles.
	ParticleDensity(Particles *ps=NULL, const std::string& name=std::string("ParticleDensity"));

	/// Callback for particle addition.
	virtual void particleAdded();

	/// Callback for particle removal.
	virtual void particleRemoved(unsigned int i);

	/// this gives a calling convertion for children to get the density of index i
	virtual double getDensity(unsigned int i);
	
	/// this gives a calling convertion for children to set the density of index i
	virtual void setDensity(unsigned int i, double den);

	virtual void clear();
};

#endif
