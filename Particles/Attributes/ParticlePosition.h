/*
@file ParticlePosition.h
@auther wensu
this attribute stores position of the particles
*/

#ifndef PARTICLEPOSITION_H
#define PARTICLEPOSITION_H

#include "Particles.h"
#include "ParticleAttribute.h"

class ParticlePosition : public ParticleAttribute {
public:

	MAKE_PARTICLESTUFF_NAME();
//	static std::string registry_name;
//	virtual const std::string getClass() { return registry_name; }

	/// ParticlePosition
	std::vector<gmVector3> x; 

	std::string objFilename;

	/// Add particle orientation to a system of particles.
	ParticlePosition(Particles *ps=NULL, const std::string& name=std::string("ParticlePosition"), bool xparam = true);

	/// Callback for particle addition.
	void virtual particleAdded(); 

	/// Callback for particle removal.
	void virtual particleRemoved(unsigned int i);

	/// this gives a calling conversion for children to get the position of index i
	virtual gmVector3 getPosition(unsigned int i);
	
	/// this gives a calling convertion for children to set the position of index i
	virtual void setPosition(unsigned int i, gmVector3 p);

	virtual void clear();

	bool changed;
};

#endif
