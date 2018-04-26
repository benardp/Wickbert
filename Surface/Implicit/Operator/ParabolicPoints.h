/**
* Declaration of a new implicit defined as the zero set of \f$ \nabla(F) \dot V \f$, 
 * where \f$ F \f$ is an implicit and  \f$ V \f$ is vector.
 * @file ParabolicPoints.h
 *@date 19 May. 2005
 * @author Matei N. Stroila
 * @remarks Updated to DirectionalDerivative 2 Dec 2005 -jch
 */

#ifndef PARABOLICPOINTS_H
#define PARABOLICPOINTS_H

#include "UnaryOp.h"

class ParabolicPoints : public UnaryOp
{
public:
	ParabolicPoints(Implicit *f = NULL); ///< Explicit constructor.
	double proc(const gmVector3 & x);
	gmVector3 grad(const gmVector3 & x) { return Implicit::grad(x); }
	gmMatrix3 hess(const gmVector3 & x) { return Implicit::hess(x); }
	virtual ~ParabolicPoints(){}; 
	
	MAKE_NAME();
};

#endif

