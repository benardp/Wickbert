/**
 * Declaration of a new implicit defined as the zero set of w^T H w/(w^T w)
 * where w = v - v.n n is the tangent component of the suppied view vector.
 * @file RadialCurvature.h
 * @date 22 Dec. 2005.
 * @author John C. Hart based on the thesis of Matei Stroila.
 */

#ifndef RADIALCURVATURE_H
#define RADIALCURVATURE_H

#include "UnaryOp.h"

/* An implicit operator that defines as the zero set of w^T H w/(w^T w)
 * where w = v - v.n n is the tangent component of the suppied view vector.
 */
class RadialCurvature: public UnaryOp
{
public:
	gmVector3 m_v;	///< direction in which to differentiate
	bool m_den;		///< divide by ||w||^2
	bool m_proj;	///< project view vector onto surface tangent plane

	RadialCurvature(Implicit *f = NULL, gmVector3 v = gmVector3(0,0,1)); ///< Explicit constructor.
	double proc(const gmVector3 & x);
	gmVector3 grad(const gmVector3 & x) { return Implicit::grad(x); }
	gmMatrix3 hess(const gmVector3 & x) { return Implicit::hess(x); }
	
	MAKE_NAME();
};

#endif
