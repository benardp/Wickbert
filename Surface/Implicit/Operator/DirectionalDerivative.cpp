/**
 * Implementation of a new implicit defined as the zero set of \f$ \nabla(F) \dot V \f$, 
 * where \f$ F \f$ is an implicit and  \f$ V \f$ is vector.
 * @file DirectionalDerivative.cpp
 * @date 19 May. 2005
 * @author Matei N. Stroila
 * @remarks
 */

#include "DirectionalDerivative.h"
#include "SurfParam.h"

REGISTER_IMPLICIT(DirectionalDerivative,"UnaryOp:DirectionalDerivative");

DirectionalDerivative::DirectionalDerivative(Implicit *f,  gmVector3 v)
{
	m_f = f;

	new SurfParamgmVector3(this, &m_v, v, "v", "direction", "Derivative direction.");
}

double DirectionalDerivative::proc(const gmVector3 & x)
{
	if (m_f) 
		return dot(m_f->grad(x), m_v);
	else
		return 0;
}

/** Computed via product rule: grad(grad f . v) = Hf v + 0 (since grad v = 0)
 */
gmVector3 DirectionalDerivative::grad(const gmVector3 & x)
{
	if (m_f)
		return m_f->hess(x) * m_v;
	else
		return gmVector3();
}