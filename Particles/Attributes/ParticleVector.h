/**
* Declaration of the attribute representing  per particle vector data.
* @file ParticleVector.h
* @date 27 June. 2005
* @author Matei N. Stroila
*/


#ifndef PARTICLEVECTOR_H
#define PARTICLEVECTOR_H

#include "Particles.h"
#include "ParticleAttribute.h"

class ParticleVector : public ParticleAttribute {
public:
	
	MAKE_PARTICLESTUFF_NAME();
		
	
	/// Add particle vectorial data to a system of particles.
	ParticleVector(Particles *ps=NULL, const std::string& name=std::string("ParticleVector"));
	
	/// Callback for particle addition.
	void particleAdded();
	
	/// Callback for particle removal.
	void particleRemoved(unsigned int i);
	
	/// this gives a calling convertion for children to get the position of index i
	virtual inline const gmVector3 getVector(unsigned int i){return data[i];};
	
	/// this gives a calling convertion for children to set the position of index i
	virtual inline void setVector(unsigned int i, const gmVector3& p){data[i]=p;};
	
	virtual void clear();
	
	int qlenpp();
	void getqpp(double *q, int i);
	void setqpp(double *q, int i);
	void qnamepp(char **qn);
	
private:

	/// Particle Vector Attribute Data
	std::vector<gmVector3> data; 

};

#endif
