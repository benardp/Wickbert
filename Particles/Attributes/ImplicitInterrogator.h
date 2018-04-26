/**
* Definition of ImplicitInterrogator.
* @file ImplicitInterrogator.h
* @author Ed Bachta some modifications by Elmar Eisemann
*/

#ifndef __IMPLICITINTERROGATOR_H__
#define __IMPLICITINTERROGATOR_H__

#include "Surface/Implicit/Implicit.h"
#include "ParticleAttribute.h"

#include <valarray>

class ParticlePosition;	// for caching implicit values

/**
* The ImplicitInterrogator attribute allows a particle system
* to interrogate an implicit function. It maintains
* a pointer to the Implicit object and an array
* of parameter change values.
*/
class ImplicitInterrogator : public ParticleAttribute
{
protected:
	/// An implicit function to interrogate.
	Implicit* _implicit;
	ParticlePosition *_position;

public:
	MAKE_PARTICLESTUFF_NAME();

	virtual double proc(unsigned int) const;
	virtual gmVector3 grad(unsigned int) const;
	virtual double Fx(unsigned int) const;
	virtual double Fy(unsigned int) const;
	virtual double Fz(unsigned int) const;
	virtual gmVector3 normal(unsigned int) const;

	//function exploits symmetry of the hessian
	virtual gmMatrix3 hess(unsigned int) const;
	virtual double Fxx(unsigned int) const;
	virtual double Fxy(unsigned int) const;
	virtual double Fxz(unsigned int) const;
	virtual double Fyy(unsigned int) const;
	virtual double Fyz(unsigned int) const;
	virtual double Fzz(unsigned int) const;


	//gaussian curvature
	virtual double kg(unsigned int i) const;
	//mean curvature
	virtual double km(unsigned int i) const;




	// temperary store the index points to the surface for file saving and loading
	int index;

	/** Need to change to find implicit by name */
	std::string impname;

	/// Change in implicit function parameters.
	std::valarray<double> dqdt;
	

	/** 
	* Assign the implicit interrogator attribute to a particle system.
	* @param ps 	  The owning particle system.
	* @param implicit The implicit function.
	* @param name	  The name of this interrogator.
	*/
	ImplicitInterrogator(Particles* ps=NULL, Implicit* implicit = NULL, const std::string &name = std::string("ImplicitInterrogator"));

	virtual void attachAttributes();
	virtual void setParticleSystem(Particles *new_ps);

	/** Set the implicit to use.
	*/
	void setImplicit(Implicit* implicit);
	void setImplicit(std::string s);
	void empty();
	
	/** Get the current implicit.
	*/
	Implicit* getImplicit()
	{
		return _implicit;
	}
	
};

class PSParamImplicitName : public PSParamString
{
public:
	virtual void set(std::string s) {
		*pv = s;
		ImplicitInterrogator *ii = dynamic_cast<ImplicitInterrogator *>(stuff);
		if (ii) ii->setImplicit(s);
	}

	PSParamImplicitName(ParticleStuff *parent, std::string *v, std::string init,const std::string& sn,const std::string& n,const std::string& d = "")
	: PSParamString(parent,v,init,sn,n,d) {}
};

#endif

