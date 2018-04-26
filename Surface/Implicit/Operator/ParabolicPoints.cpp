/**
 * Implementation of a new implicit defined as the zero set of \f$ \nabla(F) \dot V \f$, 
 * where \f$ F \f$ is an implicit and  \f$ V \f$ is vector.
 * @file DirectionalDerivative.cpp
 * @date 19 May. 2005
 * @author Matei N. Stroila
 * @remarks
 */

#include "ParabolicPoints.h"

REGISTER_IMPLICIT(ParabolicPoints,"UnaryOp:ParabolicPoints");

ParabolicPoints::ParabolicPoints(Implicit *f)
: UnaryOp(f)
{
}

double ParabolicPoints::proc(const gmVector3 & x)
{
	if (m_f) 
		return m_f->gaussianCurvature(x);
	else
		return 0;
}
