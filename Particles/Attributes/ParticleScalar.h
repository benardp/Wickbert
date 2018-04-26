/**
 * Declaration of ParticleScalar.
 * @file ParticleScalar.h
 * @author Tony Kaap
 */

#ifndef PARTICLESCALAR_H
#define PARTICLESCALAR_H

#include "ParticleAttribute.h"


/**
 * ParticleScalar is a ParticleAttribute which represents
 * any scalar-valued attribute of a particle system.  This
 * was added to prevent the multiplication of atributes that 
 * were all just a per-particle scalar value.
 */
class ParticleScalar : public ParticleAttribute {

public:

	MAKE_PARTICLESTUFF_NAME();

private:
	/// Particle scalar values
	std::vector<double> t;
	int numElements;

	std::string *str;
	Particles * par_ref;
public:
	/// Add a particle scalar attribute to a system of particles.
	// print flag set to false until IO routine implementation completed -jh
	ParticleScalar(Particles *ps=NULL, const std::string& name=std::string("ParticleScalar"), const bool& printFlag=false);

	void PSPtest();

	void clear();

	void reset();

		/// Callback for particle addition.
	virtual void particleAdded();

	std::string rbfFilename; // filename used to load rbfConstraints
	void loadFromFile(std::string fileName);

	/// Callback for particle removal.
	virtual void particleRemoved(unsigned int i);

	bool changed;
	
	// Returns the vector of stored values
	std::vector<double> getScalar();
	
	// Returns one element of the vector of values
	double getScalar(int i);
	
	void setScalar(const int& i, const double& value);
	int getNumElements() {return t.size();}
};

#endif
