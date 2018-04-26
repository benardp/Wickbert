/**
* Declaration of a new implicit defined as the zero set of \f$ \nabla(F) \dot V \f$, 
 * where \f$ F \f$ is an implicit and  \f$ V \f$ is vector.
 * @file DirectionalDerivative.h
 *@date 19 May. 2005
 * @author Matei N. Stroila
 * @remarks Updated to DirectionalDerivative 2 Dec 2005 -jch
 */

#ifndef DIRECTIONALDERIVATIVE_H
#define DIRECTIONALDERIVATIVE_H

#include "UnaryOp.h"
#include "libgm/gm.h"

class DirectionalDerivative : public UnaryOp
{
public:
	gmVector3 m_v;	///< direction in which to differentiate

	DirectionalDerivative(Implicit *f = NULL, gmVector3 v = gmVector3(0,0,1)); ///< Explicit constructor.
	double proc(const gmVector3 & x);
	gmVector3 grad(const gmVector3 & x);
	virtual ~DirectionalDerivative(){}; 
	
	MAKE_NAME();
};

#endif

