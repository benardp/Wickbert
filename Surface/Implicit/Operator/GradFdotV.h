/**
* Declaration of a new implicit defined as the zero set of \f$ \nabla(F) \dot V \f$, 
 * where \f$ F \f$ is an implicit and  \f$ V \f$ is vector.
 * @file GradFdotV.h
 *@date 19 May. 2005
 * @author Matei N. Stroila
 * @remarks
 */

#ifndef GRADFDOTV_H
#define GRADFDOTV_H

#include "Surface/Implicit/Implicit.h"

class GradFdotV : public Implicit
{
public:
	GradFdotV();
	~GradFdotV(){}; 
	
	double proc(const gmVector3 & x);
	gmVector3 grad(const gmVector3 & x);
		
	MAKE_NAME();
		
private:
		  
	Particles* p;
	Implicit *myF; ///< The input implicit.	
	
};

#endif

