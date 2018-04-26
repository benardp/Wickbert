/**
* Definition of ImplicitInterrogatorCached.
* @file ImplicitInterrogatorCached.h
* @author Elmar Eisemann
*/

#ifndef __IMPLICITINTERROGATORCACHED_H__
#define __IMPLICITINTERROGATORCACHED_H__

#include "ImplicitInterrogator.h"

/**
* The ImplicitInterrogatorCached attribute allows a particle system
* to interrogate an implicit function. It maintains
* a pointer to the Implicit object nevertheless this pointer should NOT be used for the 
* particles having this attribute.
* Values are cached and thus can be reused efficiently.
*/
class ImplicitInterrogatorCached : public ImplicitInterrogator
{
protected:
	/// An implicit function to interrogate.
	bool _calculateProc;
	std::vector<double> _proc;	///< cache of per-particle implicit value
	bool _calculateGrad;
	std::vector<gmVector3> _grad;	///< cache of per-particle gradient
	bool _calculateHessian;
	std::vector<double> _fxx;		///< cache of per-particle second partial deriv
	std::vector<double> _fxy;		///< cache of per-particle second partial deriv
	std::vector<double> _fxz;		///< cache of per-particle second partial deriv
	std::vector<double> _fyy;		///< cache of per-particle second partial deriv
	std::vector<double> _fyz;		///< cache of per-particle second partial deriv
	std::vector<double> _fzz;		///< cache of per-particle second partial deriv

	bool _calculateKg;
	std::vector<double> _kg;		///< cache of per-particle gaussian curvature
	bool _calculateKm;
	std::vector<double> _km;		///< cache of per-particle mean curvature


public:
	MAKE_PARTICLESTUFF_NAME();


	//values are populated by the Implicit Interrogator itself.
	//the idea is that this class should cache the values and then all other 
	//classes just recover them. In fact the implict in the interrogator should NOT be accessible!
	//This was a design flaw in Wickbert... mostly this is fixed now.- ELMAR

	virtual double proc(unsigned int) const;

	virtual gmVector3 grad(unsigned int) const;
	virtual double Fx(unsigned int) const;
	virtual double Fy(unsigned int) const;
	virtual double Fz(unsigned int) const;

	virtual double kg(unsigned int) const;
	virtual double km(unsigned int) const;

	//function exploits symmetry of the hessian
	virtual gmMatrix3 hess(unsigned int) const;

	virtual double Fxx(unsigned int) const;
	virtual double Fxy(unsigned int) const;
	virtual double Fxz(unsigned int) const;
	virtual double Fyy(unsigned int) const;
	virtual double Fyz(unsigned int) const;
	virtual double Fzz(unsigned int) const;




	/** 
	* Assign the implicit interrogator attribute to a particle system.
	* @param ps 	  The owning particle system.
	* @param implicit The implicit function.
	* @param name	  The name of this interrogator.
	*/
	ImplicitInterrogatorCached(Particles* ps=NULL, Implicit* implicit = NULL, const std::string &name = std::string("ImplicitInterrogatorCached"));

	virtual void setParticleSystem(Particles *new_ps);

	virtual void prepare();

	virtual void particleRemoved(unsigned int i);
	virtual void particleAdded();
};

#endif

