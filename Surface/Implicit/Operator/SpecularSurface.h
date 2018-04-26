/**
* Declaration of a specular implicit 
 * @file SpecularSurface.h
 * @date 16 March. 2006
 * @author Matei N. Stroila
 * @remarks
 */

#ifndef SPECULARSURFACE_H
#define SPECULARSURFACE_H

#include "Surface/Implicit/Implicit.h"

class SpecularSurface : public Implicit
{
public:
	SpecularSurface();	
	SpecularSurface(double _shine, Particles* _p, Implicit* _myF = NULL) : shine(_shine), p(_p), myF(_myF){};	
	double proc(const gmVector3 & x);
//	gmVector3 grad(const gmVector3 & x);
	
	void setImplicit(Implicit* _myF){myF = _myF;};
	void setParticles(Particles* _p){p = _p;};
	~SpecularSurface(){}; 
	
	MAKE_NAME();

private:
	
	double shine;
	Particles* p;
	Implicit* myF; 
	
};

#endif

