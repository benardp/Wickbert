/**
* Declaration of an implicit 
 * @file Icosahedral.h
 * @date 16 March. 2006
 * @author Matei N. Stroila
 * @remarks
 */

#ifndef Icosahedral_h
#define Icosahedral_h

#include "Surface/Implicit/Implicit.h"

class Icosahedral : public Implicit
{
public:
	Icosahedral();	
	double proc(const gmVector3 & x);
	gmVector3 grad(const gmVector3 & x);
	gmMatrix3 hess(const gmVector3 & x); 
	
	~Icosahedral(){}; 
	
	MAKE_NAME();
		
};

#endif

