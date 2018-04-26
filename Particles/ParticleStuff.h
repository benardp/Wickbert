#ifndef PARTICLESTUFF_H
#define PARTICLESTUFF_H


/** \page particlespage The Particle System Library
 *
 * The Particle System Library is based on a Particles object. A Particles object is a
 * homogeneous collection of particles that share the same behavior and appearance.
 * For example, "floater" particles used to interrogate an implicit surface would
 * be described by one Particles object, and "control particles" used to deform an implicit
 * surface would be a second Particles object.
 * Several Particles objects (for instance "floaters" and "control particles") can be
 * collected together into a ParticleSystem object, whose member ParticleSystem::particles
 * is a vector of (pointers to) Particles objects.
 *
 * Each Particles object is composed of a collection of behaviors, attributes and shaders.
 * The behaviors of a Particles object are stored in the member Particles::behaviors
 * which is of type ::ParticleBehaviors, a vector of ParticleBehavior objects that ensures
 * behaviors are performed in the proper order. The shared and per-particle data attributes are
 * stored in the member Particles::attributes which is of type ::ParticleAttributes, a map of
 * ParticleAttribute objects that makes it easy for a behavior to find attribute data by name.
 * Finally, the appearance properties of the particles are stored in the member shaders
 * which is of type ParticleShaders, a vector of ::ParticleShaders that too ensures shaders are
 * applied in the proper order.
 *
 * The classes ParticleAttribute, ParticleBehavior and ParticleShader contain some shared
 * functionality which each inherits from the ParticleStuff class.
 *
 * \section Behaviors
 *
 * The motion of particles is represented by an ordered collection of ParticleBehavior objects.
 * Each ParticleBehavior object represents a different component of the behavior of particles.
 * The contribution of a ParticleBehavior object to the behavior of a particle is divided into
 * four computational steps, computed by the four member functions:
 * \li ParticleBehavior::applyForce() - add the behavior's force to each particle
 * \li ParticleBehavior::applyConstraint() - remove the "illegal" component of the total force on a particle
 * \li ParticleBehavior::integrate() - update the particles position (state) based on the constrained force
 * \li ParticleBehavior::cleanup() - create and destroy particles based on population dynamics
 *
 * This sequence of computations is applied by Particles::fullUpdate().
 * This computes the behavior of the Particles object by first calling applyForce() on all
 * ParticleBehaviors (in order), then all behaviors' applyConstraint(), then all integrate()
 * members, then all cleanup() members.
 *
 * \section Shaders
 *
 * Like behaviors, the Particles object contains an ordered collection of ParticleShader
 * objects which are responsible for the appearance of (and user interaction with) the
 * particles. Each subclass of ParticleShader implements appearance by redefining the
 * ParticleShader::drawParticle(int s) member. Each Particles object is rendered by calling
 * ParticleShader::draw(int s) for each of its ParticleShader objects, in order.
 *
 * \section Attributes
 *
 * Per-particle and shared data elements of a Particles object is stored in a
 * ParticleAttribute object. These objects are accessed by name, and they typically contain
 * member elements that compartmentalize various quantities used by the behaviors and shaders.
 *
 * Behaviors and shaders can access the data held in a ParticleAttribute via the
 * templated ParticleStuff::attachAttribute() class. The behavior/shader calls this member
 * in its constructor with two parameters. The first is a pointer of the type of the
 * desired ParticleAttribute that will afterward point to the attribute object. The second
 * parameter is the name of the desired attribute, of type std::string. If attachAttribute
 * cannot find the desired attribute in the current Particles object, it will create it.
 *
 * \section Parameters
 *
 * Following the same interface as the Implicit object, the attributes, behaviors and shaders
 * of a Particles object have parameters accessed by the following interface:
 * \li ParticleStuff::qlen() - returns the number of parameters
 * \li ParticleStuff::getq(double *q) - returns the parameters in a preallocated double array q
 * \li ParticleStuff::setq(double *q) - sets the parameters from double array q
 * \li ParticleStuff::qname(char **qn) - returns the names of the parameters in an array preallocated with
 * character pointers
 *
 * The above parameters are uniform across all particles, e.g. the magnitude of
 * a penalty force.
 * There are also per-particle counterparts
 * \li ParticleStuff::qlenpp(),
 * \li ParticleStuff::getqpp(double *q, int i),
 * \li ParticleStuff::setqpp(double *q, int i), and
 * \li ParticleStuff::qnamepp(char **qn)
 * which apply to data associated with each particle, such as for example the
 * particle's orientation.
 *
 * \section Construction
 *
 * When attributes, behaviors and shaders are created, they need to be integrated into a
 * Particles object. Using a constructor that accepts a pointer to the Particles object as a
 * parameter automatically integrates the attribute/behavior/shader into the Particles object.
 * For example
 * 
 * new ParticleAge(p)
 *
 * will add the attribute ParticleAge to the Particles object p.
 * 
 * If you construct an attribute/behavior without a particle system (for
 * example, from the factory via NEW_PARTICLESTUFF), then you need to call
 * ParticleStuff::setParticleSystem to connect it to the particle system.
 * If your attribute or behavior knows which particle system it references,
 * but it is not known by that particle system, then call its
 * ParticleStuff::addMeTo(Particles *) method.
 */

/** Parent class of shared functionality between
 * ParticleAttribute and ParticleBehavior classes.
 *
 * When creating a new ParticleStuff object, you will need to
 * attach it to a particle system. This can be accomplished in
 * one of two ways:
 * \li Use a constructor with a Particles * argument, or
 * \li Construct the object and call setParticleSystem.
 * 
 * \note There is no ParticleStuff(Particles *,...) constructor. This is
 * because constructors with a Particles parameter need to call the specific
 * setParticleSystem for a behavior or attribute, which is not accessible
 * at the ParticleStuff level of the class hierarchy.
 */

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include "ParticleStuffParameters.h"
#include "ParticleSystem.h"	// for AttachAttribute

typedef std::vector<double>         DoubleVector;
typedef std::vector<std::string>	NameVector;

class ParticleStuff
{

public:
	/// Vector of particle stuff parameters (see PSParam)
	PSParams params;

	/// Vector of per-particle parameters (see PSParamPerParticle)
	PSParamsPerParticle perparticle;

	AttachedAttributes attachedattributes;

	/// Pixmap icon for this object
	char **xpm;

	/** Returns a pointer to an attribute variable given the
	 * reference string "<attributename>:<variablename>".
	 */
	void *getAttributeVariable(std::string ref);

	/** Returns a pointer to my variable given the variable's shortname */
	void *getVariable(std::string ref);

// I believe this is automatically added by the factory
// Included registry_name again to implement getClass()
#if 1
	/// The class name of the object according to the factory.
	std::string registry_name;
#endif

	/** Just a placeholder. Didn't make pure virtual in case newer
	 * development classes haven't yet been entered into the factory.
	 * Should be overwritten. Otherwise will return
	 * ParticleStuff::registry_name, which will be NULL.
	 */
	virtual const std::string getClass();
	
	/** Set the name of the ParticleStuff
	 */
	virtual void setName(const std::string &new_name);

	/** Return the name of the ParticleStuff
	 */
	virtual std::string getName();

	/// Reference to the owning particle system.
	Particles *ps;

	/// The name of this object.
	std::string name;
	
	/** Default constructor
	 * Constructor to set particle system and name.
	 */
	ParticleStuff(Particles *p=NULL, const std::string &new_name=std::string("ParticleStuff"));
	virtual ~ParticleStuff() {}

	/** Set the particle system the ParticleStuff points to.
	 * Also adds the ParticleStuff to the particle system.
	 * If this is a new particle system, this deletes the ParticleStuff
	 * from the existing particle system.
	 *
	 * 2002-10-26 Wen
	 * Pure virtual function at the ParticleStuff level. Made concrete
	 * in ParticleAttribute, ParticleBehavior, etc.
	 * Not pure any more, since all the children do the same thing,
	 * it is easier just to inherent this one.
	 * 
	 * 2003-11-05 Wen
	 * Very important, every one MUST call this before using this class.
	 * This is used to initialize the system. Since the factory can create
	 * a ParticleStuff with a ps=NULL, we must handle initialization step
	 * as a separate function, in order to not calling the same code twice
	 * I decided not to call this function in all the constructors.
	 */
	virtual void setParticleSystem(Particles *new_ps);
	
	/** Add the stuff to a particle system. This should be avoided
	 * as it does not set ParticleStuff::ps. Use
	 * setParticleSystem(Particles *)
	 * instead.
	 */
	virtual void addMeTo(Particles *new_ps) = 0;

	// all particle stuff may need attributes, even attributes themselves
	virtual void attachAttributes() = 0;

	/** Remove the stuff from the current particle system.
	 * Does not destroy stuff. Just unlinks it from the current
	 * particle system.
	 */
	virtual void removeMe() = 0;

	/** Move me to after the indicated stuff.
	 * Defined by behaviors and shaders, which are order dependent. Attributes
	 * are sorted alphabetically and do not depend on their order.
	 */
	virtual bool moveAfter(ParticleStuff *stuff) { return false; }

	/** Called after particles are added to the system.
	*/ 
	virtual void particleAdded() {}

	/** Called after particle are removed from the system.
	*/
	virtual void particleRemoved(unsigned int i) {}

	/** 
	* Print attribute info.
	*/
	virtual void printContent();

	/** @name Global Parameter Manipulation */
	/*@{*/

	/** Returns length of particle stuff parameters.
	 * Defaults to no parameters.
	 */
	virtual int qlen();

	/** Get the parameters.
	 * Defaults to no parameters.
	 */
	virtual void getq(double *);

	/** Get the parameters.
	 * Defaults to no parameters.
	 */
	virtual void getq(DoubleVector &q);

	/** Set particle stuff parameters.
	 * Defaults to no parameters.
	 */
	virtual void setq(double *);

	/** Set particle stuff parameters stored in a DoubleVector.
	 */
	virtual void setq(DoubleVector &q);

	/** Get the names of the parameters.
	 * Defaults to no parameters.
	 */
	virtual void qname(char **);

	/** Get the names of the parameters in a name vector.
	 */
	virtual void qname(NameVector &q);

	/** Get the short names of the parameters in a name vector.
	 */
	virtual void qshortname(NameVector &q);

	/** Returns a tool-tip for the parameter i.
	 */
	virtual char *qtip(int i);

	/** Returns a tool-tip for the ParticleStuff
	 */
	virtual char *tip();

	/*@}*/

	/** @name Per-Particle Parameter Manipulation */
	/*@{*/

	/** # of parameters per particle
	 */
	virtual int qlenpp();

	/** Get the per-particle parameters.
	 * Defaults to no parameters.
	 * Second parameter is the particle index
	 */
	virtual void getqpp(double *, int);
	virtual void getqpp(DoubleVector &q, int);

	/** Set per-particle parameters.
	 * Defaults to no parameters.
	 * Second parameter is the particle index
	 */
	virtual void setqpp(double *, int);
	virtual void setqpp(DoubleVector &q, int);

	/** Get the names of the per-particle parameters.
	 * Defaults to no parameters.
	 */
	virtual void qnamepp(char **);
	virtual void qnamepp(NameVector &q);

	/*@}*/

	/** Attach an attribute to the behavior.
	 * Templated by attribute.
	 * @param attr A pointer to a subclass of ParticleAttribute. Set to zero if ps = 0.
	 * @param name The string name of the attribute in case it is different than the classname.
	 * do a linear search for a dynamic_cast<className> that is not NULL, if no matching name
	 * is found. If that none exist, create a new className.
	 *
	 * A side effect is that if an attribute 
	 * of name <name> is not found, then <name> is replaced with the name of the
	 * first attribute of appropriate type that is found. (Matei, this means don't
	 * change &name back to a const! -jch)
	 */
	template<class attrclass>
		void attachAttribute(attrclass* &attr, std::string& name);

	void error(std::string e);
};

std::ostream &operator<<(std::ostream &out, ParticleStuff *stuff);

template<class attrclass>
void ParticleStuff::attachAttribute(attrclass* &attr, std::string &name)
{
	attr=NULL;
	Particles *p = ps;
	std::string attrname = name;

	/* Here is where we should implement the ability to attach an attribute from a 
	   different particles */
	size_t colon = name.find(":");
	if (colon != std::string::npos) {
		std::string pname = name.substr(0,colon);
		attrname = name.substr(colon+1);
		// find particles in current system named pname
		p = p->particleSystem->findParticles(pname);
	}
	
	//if no attribute found, use the current particle system
	if (!p) p = ps;

	// find by name first
	std::map<std::string, ParticleAttribute *>::iterator it = p->attributes.find(attrname);
	if (it != p->attributes.end())
		attr = dynamic_cast<attrclass *>(it->second);
	else
	{
		// not found, linear search for any of the same type
		for(it = p->attributes.begin(); it != p->attributes.end(); ++it)
		{
			attr = dynamic_cast<attrclass *>(it->second);
			if (attr) {
				// rename to found attribute
				name = attr->name;
				break;
			}
		}
		if (attr==NULL)
		{
			// still not found create it
			attr = new attrclass(p);
			// initialize if just created
			attr->setParticleSystem(p);
		}
	}
}


#endif
