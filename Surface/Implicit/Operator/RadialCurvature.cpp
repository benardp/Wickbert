/**
 * Implementation of a new implicit defined as the zero set of w^T H w/(w^T w)
 * where w = v - v.n n is the tangent component of the suppied view vector.
 * @file RadialCurvature.h
 * @date 22 Dec. 2005.
 * @author John C. Hart based on the thesis of Matei Stroila.
 */

#include "RadialCurvature.h"

REGISTER_IMPLICIT(RadialCurvature,"UnaryOp:RadialCurvature");

RadialCurvature::RadialCurvature(Implicit *f,  gmVector3 v)
{
	m_f = f;

	new SurfParamgmVector3(this, &m_v, v, "v", "direction", "View direction.");
	new SurfParamBool(this, &m_den, false, "den", "denominator",
		"Scale result to report actual curvature.");
	new SurfParamBool(this, &m_proj, true, "proj", "project",
		"Project view vector into surface tangent plane.");
}

double RadialCurvature::proc(const gmVector3 & x)
{
	if (!m_f) return 0;

	gmVector3 gf = m_f->grad(x);
	double gfmag2 = gf.lengthSquared();
	gmVector3 w = (m_proj && gfmag2 == 0.0) ? m_v : m_v - dot(m_v,gf)*gf/gfmag2;
	double wmag2 = m_den ? w.lengthSquared() : 1.0;
	double kr = (wmag2 == 0.0) ? 0.0 : dot(m_f->hess(x)*w,w)/wmag2;

	return kr;
}

