/*! \file Surface.h
 *  \author John C. Hart
 *  Describes the surface class, which is everything this library supports, including implicits and subdivs.
 *
 * \note Mesh.h and Surface method m_mesh have been eliminated until Nate's new half-edge library
 * is used. Current mesh library has a nasty constructor/destructor compile bug.
 */

/** \mainpage
 * The Advanced Surface Library is designed as a core library for
 * research and development of new shape modeling tools. Its main
 * components include a broad library of Surface surface representations
 * in the \ref implicitpage object library, and a rather sophisticated particle
 * system structure for interrogating Surface surfaces, described in the
 * \ref particlespage object library.
 */

#ifndef __SURFACE_H
#define __SURFACE_H

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

/** The library uses the libgm library to represent 2- and 3-vectors.
 */
#include "Surface/Box.h"
#include "Surface/IMatrix.h"

#include "libgm/gm.h"
/** Computes the outer product between to 3-vectors. The result might
 * be the transpose of what it really should be. I haven't checked the
 * ordering.
 * @param x 3-vector factor
 * @param y 3-vector factor
 * @return 3x3 matrix containing pairwise products of x,y elements
 * @note Hate to put this here of all places, but want to keep libgm unchanged
 * and this is the only file I know will be included by everyone.
 */
gmMatrix3 outer( gmVector3 x, gmVector3 y);
/// Interval version of outer product of two Boxes
IMatrix3d outer(Box3d x, Box3d y);







//! The Surface class, which is inherited by all surfaces
/*! So far the only item that all surfaces share is a meshed representation.
*/

#include <string>
#include <vector>
#include <sstream>
#include "SurfParam.h"

class ParticleSystems;
class Surfaces;

class Surface 
{
public:
	Surface();
	virtual ~Surface() {}

	SurfParams params;

	/** A surface can contain references to other surfaces (e.g. a union would
	 * reference two components) or other particles (e.g. an RBF might reference
	 * particles for its centers). These references are by name, and the attach
	 * method looks up the names and sets the appropriate pointers. The attach
	 * method should be called anytime anything is renamed or deleted.
	 */
	virtual void attach(Surfaces *surfs, ParticleSystems *psystems) {
		params.attach(surfs,psystems);
	}

	virtual void output(std::ostream &out) {};

protected:

    /// @name Object Naming
    /// @{
    std::string registry_name;   ///< The name of the Surface class

    /**
     * The name of the Surface object when instantiated. This string should
     * be equal to the CLASS name plus a unique number so that every
     * Surface object has a unique name when first created.
     */
    std::string objectName;

    /**
     * The prefix added to the default object instance name.  This is done
     * so we know when the user has changed the object name.
     */
    static std::string defaultObjectNamePrefix;
    // Initialization of defaultObjectName at beginning of Surface.cpp
 
    /**
     * The number of the instanced object used for the object name.  This
     * number is automatically incremented upon creation of a new Surface
     * object so that every Surface has a unique name.
     */
    // Initialization of objectNumber at beginning of Surface.cpp
    static unsigned int objectNumber;

	void checkObjectName(); ///< See if the objectName has not been set yet
    /** @} end of object naming */

public:
    virtual const std::string name();   ///< Return name of 'class' of Surface
    std::string getObjectName();
    const char* getObjectNameAsCString();
    std::string getFullObjectName();
    const char* getFullObjectNameAsCString();
    void setObjectName(std::string);
    void setObjectName(char*);
    virtual void resetObjectName();
    bool objectNameIsDefault();
    void setObjectNumber(int);
	//added this is for passing events from a GUI -ms
	virtual void processEvent(const std::string&){};

};

std::ostream &operator<<(std::ostream &out, Surface *s);

//typedef std::vector<Surface *> Surfaces;

class Surfaces : public std::vector<Surface *>
{
public:
	/** This pointer gives surfaces access to the particles
	 * to allow particles to represent e.g. RBF centers.
	 */
	ParticleSystems *psystems;

	void setParticleSystems(ParticleSystems *pss) { psystems = pss; }

	Surface *withName(std::string n) {
		iterator sit;
		for (sit = begin(); sit != end(); sit++)
			if ((*sit)->getObjectName() == n)
				return *sit;
		return NULL;
	}

	/** Refreshes all surface name references.
	 *  Call this attach after any (particles, surfaces) name change or deletion.
	 */
	void attach() {
//		if(psystems) {
//			std::string t = psystems->at(0)->name;
//		}
		for (size_t i = 0; i < size(); i++)
			at(i)->attach(this,psystems);
	}
};

std::ostream &operator<<(std::ostream &out, Surfaces *s);
std::ostream &operator<<(std::ostream &out, Surface *s);

std::istream &operator>>(std::istream &in, Surfaces *surfaces); /* defined in Implicit.cpp */

#endif // __SURFACE_H

