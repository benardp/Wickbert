#ifndef PARTICLEATTRIBUTE_H
#define PARTICLEATTRIBUTE_H

#include "ParticleStuff.h"

/**
 * ParticleAttribute defines special data that is required
 * by different types of behaviors. This data is also visualized
 * by renderers via the context of behaviors.
 */
class ParticleAttribute : public ParticleStuff {

public:

	/** Default constructor
	* Assign a name and attach to a particle system.
	* @param ps   The owning Particles object.
	* @param name The name of this attribute.
	*/
	ParticleAttribute(Particles *ps, const std::string &name);
	virtual ~ParticleAttribute() {}
	//virtual void setParticleSystem(Particles *new_ps);

	/** Add attribute to a particle system.
	 * This does not set the particle system of the attribute.
	 * Use ParticleAttribute::setParticleSystem to integrate an attribute
	 * into a particle system.
	 */
	virtual void addMeTo(Particles *new_ps);

	virtual void setName(const std::string &new_name);

	/** This member function can be replaced with specialized versions
	 * that load or compute per-particle values.
	 * \default nop.
	 */
	virtual void prepare() { }

	/** Remove the attribute from the current particle system.
	 * Does not destroy the attribute. Just unlinks it from the current
	 * particle system.
	 */
	virtual void removeMe();
 
	virtual void clear() { }

	virtual void attachAttributes() { attachedattributes.attach(); }

	virtual void attachTo(ParticleStuff *stuff);


	/** iterations is a counter used to keep track of how many times member functions are called */
	int iterations;
	/** total_time used to profile how much wall clock time used by each behavior. 
	* measured is the time to prepare the attribute.
	*/

	double total_time;
};

#endif
