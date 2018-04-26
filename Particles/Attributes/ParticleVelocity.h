/*
@file ParticleVelocity.h
@auther wensu
this attribute stores position of the particles
*/

#ifndef PARTICLEVELOCITY_H
#define PARTICLEVELOCITY_H

#include "../Particles.h"
#include "../ParticleAttribute.h"

class ParticleVelocity : public ParticleAttribute {

public:

	MAKE_PARTICLESTUFF_NAME();

	/// ParticleVelocity
	std::vector<gmVector3> v;

	/// Add paritcle orientation to a system of particles.
	ParticleVelocity(Particles *ps=NULL, const std::string& name=std::string("ParticleVelocity"));

	/// Callback for particle addition.
	virtual void particleAdded();

	/// Callback for particle removal.
	virtual void particleRemoved(unsigned int i);

	int qlenpp();
	void getqpp(double *q, int i);
	void setqpp(double *q, int i);
	void qnamepp(char **qn);

};

#endif
