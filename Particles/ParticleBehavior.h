#ifndef PARTICLEBEHAVIOR_H
#define PARTICLEBEHAVIOR_H

#include <time.h>
#include "ParticleStuff.h"

/** 
 * ParticleBehavior defines a behavior of a Particles object.
 *
 * Behaviors may use one or more ParticleAttributes in determining the effect
 * on a particle system. A behavior normally checks for existing attributes
 * in initialization, and creates attributes needed if they are not found.
 *
 * When evaluating behaviors, the order is all applyForce(), then all
 * applyConstraint(), then all integrate(), then all cleanup().
 */
class ParticleBehavior : public ParticleStuff {

public:
	// everyone needs position attribute
	//ParticlePosition *position;
	// everyone needs velocity behavior
	//ParticleVelocity *velocity;

	/** Default constructor
	* Attach to a particle system.
	* @param ps   The owning Particles object.
	* @param name The name of this behavior.
	*/
	ParticleBehavior(Particles* ps=NULL, const std::string &name=std::string("ParticleBehavior"));

	/** Add the behavior to a particle system.
	 * Clobbers any behavior in new_ps with the same name.
	 * Does not update the behaviors particle system pointer.
	 * Should use ParticleBehavior::setParticleSystem to integrate
	 * a behavior into a particle system.
	 */
	virtual void addMeTo(Particles *new_ps);

	/** Removes a behavior from these particles.
		\note Particles do not own their behaviors, therefore, removed behaviors are not deleted,
		only removed from the list. You must manually delete any behaviors you remove in order
		to not leak memory.
	*/
	virtual void removeMe();

	virtual bool moveAfter(ParticleStuff *stuff);

	/** @name Behavior Application Members */
	/*@{*/

	/** Applies forces to all particles.
	 */
	virtual void applyForce() {};

	/** Constrains the current forces on all particles.
	 */
	virtual void applyConstraint() {};

	/** Applies the forces to all particles.
	 */
	virtual void integrate() {};

	/** Manages the number of particles.
	 */
	virtual void cleanup() {};

	/*@}*/

	/** Attach all the attributes needed by the behavior.
	 * Defaults to a nop.
	 */
	virtual void attachAttributes();

	/** iterations is a counter used to keep track of how many times member functions are called */
	int iterations;

	/** total_time used to profile how much wall clock time used by each behavior. */
	double total_time;
};
#endif
