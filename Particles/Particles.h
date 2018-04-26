/*
* Declaration of ParticleStuff, ParticleBehavior, ParticleAttribute and
*   finally Particles.
* @file Particles.h
* @author Ed Bachta, Wojciech Jarosz, Wen Su, John Hart
*/

#ifndef PARTICLES_H
#define PARTICLES_H

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "libgm/gm.h"
#include <vector>
#include <map>
#include <list>
#include <string>
#include <valarray>
#include <algorithm>
#include <iostream>

class ParticleStuff;
class ParticleAttribute;
class ParticleBehavior;
class ParticleShader;
class ParticleSystem;

typedef std::valarray<double>       DoubleArray;
typedef std::vector<double>         DoubleVector;
typedef std::vector<std::string>	NameVector;

/// ParticleAttributes are maintained in a map from names to attributes.
typedef std::map<std::string, ParticleAttribute*> ParticleAttributes;

/// ParticleBehaviors are maintained in a vector.
typedef std::vector<ParticleBehavior *> ParticleBehaviors;

/// ParticleShaders are maintained in a vector.
typedef std::vector<ParticleShader *> ParticleShaders;

/** Particles in an object describing a particle system.
    The particle system consists of
	\li a collection of particles, whose positions are stored in the
	    vector Particles::x and whose velocities are stored in the vector
		Particles::v,
	\li a collection of additional particle attributes stored in the
	    map Particles::attributes, and
	\li a collection of particle behaviors stored in the vector
		Particles::behaviors.

	The motion of the particle system is simulated by repeatedly calling
	fullUpdate(). This function calls all of the behaviors'
	ParticleBehavior::applyForce(), then all of the
	ParticleBehavior::applyConstraint(), then all of the\
	ParticleBehavior::integrate() (plus integration of x and v),
	and finally all of the ParticleBehavior::cleanup() members.
*/

class Particles
{
public:
	/** Name of the particle system.
	 * 
	 * For example: "floaters" or "controllers"
	 */
	std::string name;
	std::string getName() { return name; }
	void setName(const std::string& name) { this->name = name; }

	/** particleSize keeps track of the # of particles
	 * in the system. It is returned by the size() method,
	 * incremented by addParticle() and decremented by
	 * removeParticle().
	 */
	unsigned int particleSize;

	int selectedParticle;
	int selectedShader;

	//! gui can move the position or not
	bool draggable;

	//! Execute the behaviors 
	bool play;

	//! Execute the shaders
	bool show;

	//! child particles objects that are spawn by this one
	// std::vector <Particles *> childPS;

	//! Each Particles knows its parent ParticleSystem
	ParticleSystem *particleSystem;

	//! A collection of behaviors.
	ParticleBehaviors behaviors;
	
	//! A collection of attributes.
	ParticleAttributes attributes;
	
	//! A collection of attributes.
	ParticleShaders shaders;
		
	/** Timestep.
	 *  This should probably be held in a ParticleAttribute, eventually.
	 *  Currently initialized to 0.03 in the constructor.
	 */
	double dt;
	
	/** Creates a particles object with no particles. If ps is non-null, automatically
	 *  inserts the particles object into the system.
	 */
	Particles(ParticleSystem *ps=NULL,const std::string& n = "Unnamed Particles");

	/*! Searches for the behavior based on its name.
		\return	The index of the behavior in the vector.
	*/
	int findBehavior(const std::string &name);

	int findShader(const std::string &name);

	/** linear search for anything that can be cast into this type
	 */
	template<class attributeType> bool findAttribute(attributeType *&attr)
	{
		attr=NULL;
		for(ParticleAttributes::iterator it = attributes.begin(); it!=attributes.end(); ++it)
		{
			attr = dynamic_cast<attributeType *>(it->second);
			if (attr)
				return true;
		}
		return false;
	}

	/** Finds an attribute by name. Templatized for convenient casting of the result.
	 */
	template<class attributeType> attributeType *getAttribute(const std::string& name)
	{
		attributeType *found=NULL;
		std::map<std::string, ParticleAttribute *>::iterator it = attributes.find(name);
		if (it != attributes.end())
			found = dynamic_cast<attributeType *>(it->second);
		return found;
	}

	/** Finds an attribute by name. Returns as a ParticleAttribute *.
	 */
	ParticleAttribute *getAttributeGeneric(const std::string& name)
	{
		std::map<std::string, ParticleAttribute *>::iterator it = attributes.find(name);
		if (it != attributes.end())
			return it->second;
		return NULL;
	}

	template<class behaviorType> behaviorType *getBehavior(const std::string& name)
	{
		behaviorType *found=NULL;
		int index=findBehavior(name);
		if (index!=-1)
			found = dynamic_cast<behaviorType *>(this->behaviors[index]);
		return found;
	}

	template<class shaderType> shaderType *getShader(const std::string& name)
	{
		shaderType *found=NULL;
		int index=findShader(name);
		if (index!=-1)
			found = dynamic_cast<shaderType *>(this->shaders[index]);
		return found;
	}

	template <class attributeType> attributeType* getAttributeByType()
	{
		attributeType *result=0;
		std::map<std::string, ParticleAttribute *>::iterator it=attributes.begin();
		for (;it!=attributes.end();++it)
		{
			if ((result=dynamic_cast<attributeType *>(it->second)))
			{
				return result;
			}
		}
		return result;
	}

	template <class attributeType> void getAllAttributesByType(std::vector<attributeType*> * result)
	{
		assert(result);

		std::map<std::string, ParticleAttribute *>::iterator it=attributes.begin();
		for (;it!=attributes.end();++it)
		{
			if (dynamic_cast<attributeType *>(it->second))
			{
				result->push_back(dynamic_cast<attributeType *>(it->second));
			}
		}
	}





	//template<class behaviorType> behaviorType *getBehaviorByType(const std::string& ambiguity="")
	//{
	//	behaviorType *result=0;
	//	
	//	const ParticleBehaviors::iterator end=behaviors.end();
	//	for (ParticleBehaviors::iterator iter =behaviors.begin();iter!=end;++iter)
	//	{ 
	//		if (dynamic_cast<behaviorType *>(*iter))
	//		{
	//			if ((*iter)->name==ambiguity)
	//				return dynamic_cast<behaviorType *>(*iter);

	//			if (result==0)
	//				result=dynamic_cast<behaviorType *>(*iter);
	//		}
	//	}
	//	return result;
	//}

	//to make it Mac compatible and because I did not really like the solution with the name
	template<class behaviorType> behaviorType *getBehaviorByType()
	{
		behaviorType* result;
		const ParticleBehaviors::iterator end=behaviors.end();
		for (ParticleBehaviors::iterator iter =behaviors.begin();iter!=end;++iter)
		{ 
			if ((result=dynamic_cast<behaviorType *>(*iter)))
			{
				return result;
			}
		}
		return 0;
	}



	template<class behaviorType> void getAllBehaviorsByType(std::vector<behaviorType*> * result)
	{
		assert(result);
		
		const ParticleBehaviors::iterator end=behaviors.end();
		for (ParticleBehaviors::iterator iter =behaviors.begin();iter!=end;++iter)
		{ 
			if (dynamic_cast<behaviorType *>(*iter))
			{
				result->push_back(dynamic_cast<behaviorType *>(*iter));
			}
		}
	}


	//see above note...
	//template<class shaderType> shaderType *getShaderByType(const std::string& ambiguity="")
	//{
	//	shaderType *result=0;
	//	
	//	const ParticleShaders::iterator end=shaders.end();
	//	for (ParticleShaders::iterator iter =shaders.begin();iter!=end;++iter)
	//	{ 
	//		if (dynamic_cast<shaderType *>(*iter))
	//		{
	//			if ((*iter)->name==ambiguity)
	//				return dynamic_cast<shaderType *>(*iter);

	//			if (result==0)
	//				result=dynamic_cast<shaderType *>(*iter);
	//		}
	//	}
	//	return result;
	//}
	

	template<class shaderType> shaderType *getShaderByType()
	{
		shaderType *result;
		
		const ParticleShaders::iterator end=shaders.end();
		for (ParticleShaders::iterator iter =shaders.begin();iter!=end;++iter)
		{ 
			if ((result=dynamic_cast<shaderType *>(*iter)))
			{
				return result;
			}
		}
		return 0;
	}

	template<class shaderType> shaderType *getAllShadersByType(std::vector<shaderType*> *result)
	{
		assert(result);
		
		const ParticleShaders::iterator end=shaders.end();
		for (ParticleShaders::iterator iter =shaders.begin();iter!=end;++iter)
		{ 
			if (dynamic_cast<shaderType *>(*iter))
			{
				result->push_back(dynamic_cast<shaderType *>(*iter));
			}
		}
	}



	/// create a new ParticleStuff and initialize it with setParticleSystem
	template<class attrclass>
	attrclass *createParticleStuff()
	{
		attrclass *stuff=new attrclass(this);
		stuff->setParticleSystem(this);
		return stuff;
	}

	/** Removes a particle attribute, behavior or shader. Returns false if not found.
	 */
	bool remove(ParticleStuff *stuff);

	/** Removes all of the attributes.
	 */
	void removeAttributes();
	/** Removes all of the behaviors.
	 */
	void removeBehaviors();
	/** Removes all of the shaders.
	 */
	void removeShaders();

	void removeStuff();
	
	void copyFrom(Particles *p);

	unsigned int addParticle();

	//! Returns the number of particles in the system.
	unsigned int size() {return particleSize;}
	
	//! Remove a particle from the system.
	virtual void removeParticle(unsigned int i);

	//! Remove all particles
	void removeAll() { while (particleSize) removeParticle(particleSize-1); }

	/** Perform all the bookkeeping. Should be run anytime an attribute's, behavior's
	 * or shader's name is changed, or one is deleted, or a connection is changed.
	 */
	void attachAttributes();

	//! Processes one entire time step.
	void fullUpdate();

	//! Prepare attributes before processing
	void prepareAttributes();

	//! Process all forces.
	void applyForces();
	
	//! Process all constraints.
	void applyConstraints();

	//! Do all integration. Zeroes the velocity vector.
	void integrate();

	//! Clean everything up.
	void cleanup();
	
	//! Prints the contents of this particle system.
	void printContent();
	
	//! it calls the shaders to draw itself.
	void draw();

	//! call all the event processors in the shaders
	void event(int e);

	void error(std::string e);

	std::string profile();

	//! Desctructor deletes all behaviors and attribute.
	virtual ~Particles();

	/** Returns a pointer to an attribute variable or a vector of perparticle variables.
	 *  \parameter ref String of the form "<attributename>:<perparticleparametershortname>"
	 */
	void *getVariable(std::string ref);
};

std::ostream &operator<<(std::ostream &out, Particles *p);
std::istream &operator>>(std::istream &in, Particles *p);





/* This stuff below is modified from similar stuff in Implicit.h
 */

#include "tools/factory.h"

/** Shorthand type for a factory that creates particle attributes and behaviors.
 * ParticleStuffFactory is a genericFactory<ParticleStuff>.
 */
typedef genericFactory<ParticleStuff> ParticleStuffFactory;

/** Shorthand type for the ParticleStuffFactory registry.
 * This type is essentially a map of classname strings to
 * functions that create the new class.
 * \sa PARTICLESTUFF_REGISTRY
 */
typedef genericFactory<ParticleStuff>::FN_REGISTRY ParticleStuffRegistry;

/** Macro to register a new Implicit subclass.
 * \param class_name The token name of the class.
 * \param class_string The string name of the class.
 * \note This belongs in the .cpp file, right after the includes.
 */
#define REGISTER_PARTICLESTUFF(class_name,class_string) \
  int grab##class_name; \
  std::string class_name::registry_name (class_string); \
  namespace { \
	registerInFactory<ParticleStuff,class_name> registerMe(class_name::registry_name); \
  }

/*
  bool classname::compatible(ParticleStuff *stuff) { \
    return dynamic_cast<classname *>(stuff) != NULL; \
  }
*/

/** Macro to add name() and objectName functionality to all subclasses of
 * ParticleStuff.
 */
#define MAKE_PARTICLESTUFF_NAME() \
  static std::string registry_name;\
  virtual const std::string getClass() { return registry_name; }
/*  virtual void resetObjectName() {\
	std::ostringstream objectNumberName;\
	objectNumberName << ":" << ++objectNumber;\
	objectName = "~!" + name() + objectNumberName.str();\
  }
  */

/** Macro to get the registry for the ParticleStuffFactory.
 * The registry is a map of classname strings to new class constructors.
 * \return The registry of type ParticleStuffRegistry.
 */
#define PARTICLESTUFF_REGISTRY genericFactory<ParticleStuff>::instance().registry

/** Macro to construct a new subclass of Implicit.
 * \param class_string String name of the class to create.
 * \return A std::auto_ptr<Implicit> to the appropriate subclass.
 * This is pretty much just a Implicit* but the auto_ptr adds some
 * safety mechanisms to avoid memory leaks.
 */
#define NEW_PARTICLESTUFF(class_string) \
  genericFactory<ParticleStuff>::instance().create(class_string)

#endif
