/**
Declaration of ParticleSystem
@file ParticleSystem.h
@author Wen Su
*/
#ifndef PARTICLESYSTEM_H
#define PARTICLESYSTEM_H

#include <vector>
#include <iostream>

#include "Particles/Particles.h"

class Surfaces;
class ParticleSystems;

/**
 * A ParticleSystem is a collection of Particles.
 */
class ParticleSystem
{
public:
	/// the parent particle systems
	ParticleSystems* particleSystems;
	
	/** The individual Particles objects.
	 * @todo We should eliminate this member and instead derive
	 *       ParticleSystem from std::vector<Particles *> instead,
	 *       as is done with ParticleSystems. -jch
	 */
	std::vector<Particles *> particles;

	/// The name of the ParticleSystem
	std::string name;

	/// Surfaces that can be used by the particle system
	Surfaces *surfaces;
	
	/// Create a new empty particle system with the
	/// name "A Particle System".
	ParticleSystem() : name("A Particle System") 
	{ 
		surfaces = 0; 
		particleSystems = NULL;
	}	
	/// Pass the pointer to the parent ParticleSystems.
	ParticleSystem(ParticleSystems* _ps) : name("A Particle System") 
	{ 
		surfaces = 0; 
		particleSystems = _ps;
	}

	///Virtual destructor to avoid warnings
	virtual ~ParticleSystem() {}
	
	/** Find a Particles object in the system by name.
	 * Returns 0 if not found.
	 */
	Particles *findParticles(std::string &name);

	int findanderase(Particles *p) {
		std::vector<Particles *>::iterator pi;
		int foundanddeleted = 0;
		for (pi = particles.begin(); pi != particles.end(); pi++) {
			if (*pi == p) {
				(*pi)->removeStuff();	// ought to put this in the destructor
				delete (*pi);
				particles.erase(pi);
				foundanddeleted = -1;
				break;
			}
		}
		return foundanddeleted;
	}

	virtual void fullUpdate() {
		for (size_t i = 0; i < particles.size(); i++)
			particles[i]->fullUpdate();
	}

	void attachAttributes() {
		std::vector<Particles *>::iterator pi;
		for (pi = particles.begin(); pi != particles.end(); pi++)
			(*pi)->attachAttributes();
	}

	void clearSelection()
	{ for (unsigned int i = 0; i < particles.size(); i++) particles[i]->selectedParticle = -1; }

	/** A temporary list of any errors that might occur.
	 * This will likely be cleared often.
	 */
	std::vector<std::string *> errors;

	void error(std::string e) { errors.push_back(new std::string(e)); }

	std::string profile();

	/** Returns a pointer to an attribute variable or a vector of perparticle variables.
	 *  \parameter ref String of the form "<particlesname>:<attributename>:<perparticleparametershortname>"
	 */
	void *getVariable(std::string ref);
	
	//set the xrot, yrot and zoom (GUI camera information)
	inline void setEulerAnglesAndZoom(const double x,const double y,const double z)
	{
		xrot = x; yrot = y; zoom = z;
	}
	//get the xrot, yrot and zoom (GUI camera information)
	inline void getEulerAnglesAndZoom(double* x,double* y,double* z)
	{
		*x = xrot; *y = yrot; *z = zoom;
	}

	inline void pause()
	{
		std::vector<Particles *>::iterator pi;
		for (pi = particles.begin(); pi != particles.end(); pi++)
			(*pi)->play = false;
	}

	inline void play()
	{
		std::vector<Particles *>::iterator pi;
		for (pi = particles.begin(); pi != particles.end(); pi++)
			(*pi)->play = true;
	}
	
private:
	//Euler angles and zoom (GUI camera information)
	double xrot, yrot, zoom;
};

std::ostream &operator<<(std::ostream &out, ParticleSystem *ps);
std::istream &operator>>(std::istream &in, ParticleSystem *ps);

class ParticleSystems : public std::vector<ParticleSystem *>
{
public:
	ParticleSystems(){animate = false;
	rotateX = true; rotateY = false; rotateZ = false; speed = 1; zoom = -5;}
	void attachAttributes();
	void setSurfaces(Surfaces *s);
	void surfacesAttach();
	void fullUpdate();

	void clearSelection()
	{ for (unsigned int j = 0; j < size(); j++) at(j)->clearSelection(); }

	/** A temporary list of any system-wide particle errors that might occur.
	 * This will likely be cleared often.
	 */
	std::vector<std::string *> errors;

	void error(std::string e) { errors.push_back(new std::string(e)); }
	bool iserror();
	std::string errorstring();

	/** Finds a collection of particles.
	 *  \param s A string of the form "<ParticleSystemName>:<ParticlesName>" or
	 *  "ParticlesName".
	 *
	 * If only the "<ParticlesName>" is provided, it searches all ParticleSystem objects for
	 * a Particles object of that name.
	 */
	Particles *findParticles(std::string s) {
		std::string psysname,pname;
		int colon = s.find(':');
		if (colon != -1) {
			psysname = s.substr(0,colon);
			pname = s.substr(colon+1);
		} else {
			pname = s;
		}

		for (size_t j = 0; j < size(); j++) {
			if (psysname.size() && at(j)->name != psysname) continue;
			for (size_t i = 0; i < at(j)->particles.size(); i++) {
				Particles *p = at(j)->particles[i];
				if (p->name == pname) return p;
			}
		}
		return NULL;
	}

	std::string profile();

	inline void setAnimation(bool _animate, bool _rotateX, bool _rotateY, bool _rotateZ, double _speed, double _zoom, double _xinit, double _yinit){
	animate = _animate;
	rotateX = _rotateX;
	rotateY = _rotateY;
	rotateZ = _rotateZ;
	speed = _speed;
	zoom = _zoom;
	xrot = _xinit;
	yrot = _yinit;
	
	}; //if animate == true, animate the GL canvas
	void animation();

	inline void pause()
	{
		std::vector<ParticleSystem *>::iterator pi;
		for (pi = begin(); pi != end(); pi++)
			(*pi)->pause();
	}

	inline void play()
	{
		std::vector<ParticleSystem *>::iterator pi;
		for (pi = begin(); pi != end(); pi++)
			(*pi)->play();
	}

private:

	bool animate;
	bool rotateX;
	bool rotateY;
	bool rotateZ;
	double speed;
	double xrot;
	double yrot;
	double zrot;
	double zoom;

};

std::ostream &operator<<(std::ostream &out, const ParticleSystems &psystems);
std::istream &operator>>(std::istream &in, ParticleSystems *psystems);

#endif
