//*****************************************************************************
//	CSGIntersection.h: Header file for the CSGIntersection class.
//	Author: Xinlai Ni
//	1/13/2006
//*****************************************************************************

#ifndef	XNTOOLS_IMPLICIT_CSGINTERSECTION_H
#define	XNTOOLS_IMPLICIT_CSGINTERSECTION_H

#include "Surface/Implicit/Implicit.h"
#include "Surface/Implicit/Operator/BinaryOp.h"

class CSGIntersection : public BinaryOp
{
public:
	double h(double f, double g) {
		return f > g ? f : g;
	}
	double hf(double f, double g) {
		return f > g ? 1.0 : 0.0;
	}
	double hg(double f, double g) {
		return g > f ? 1.0 : 0.0;
	}
	double hff(double f, double g) {
		return 0.0;
	}
	double hfg(double f, double g) {
		return 0.0;
	}
	double hgg(double f, double g) {
		return 0.0;
	}

#ifndef INTERVAL_EVAL_ONLY
    virtual double proc(const gmVector3 & x);
    virtual gmVector3 grad(const gmVector3 & x);
	virtual gmMatrix3 hess(const gmVector3 & x);
#endif

    void procq(const gmVector3 & x, double* q);
	MAKE_NAME();
};

#endif	//	XNTOOLS_IMPLICIT_CSGINTERSECTION_H

