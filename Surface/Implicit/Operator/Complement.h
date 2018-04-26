/**
 * @file Complement.h
 * Represents the intersection of two implicit surfaces.
 * @author Steve Zelinka, Ed Bachta
 */

#ifndef _COMPLEMENT_H
#define _COMPLEMENT_H

#include "UnaryOp.h"
#include "Surface/Implicit/Algebraic/Algebraic.h"

class Complement : public UnaryOp
{
  public:
    /// Constructors
    Complement() { m_f = NULL; }
    Complement(Implicit *f) { m_f = f; }

#ifndef INTERVAL_EVAL_ONLY
    virtual double proc(const gmVector3 & x);
    virtual gmVector3 grad(const gmVector3 & x);
    virtual gmMatrix3 hess(const gmVector3 & x);
#endif

    virtual Intervald proc(const Box<double>& x);
    virtual Box3d grad(const Box<double>& x);
    virtual IMatrix3d hess(const Box<double>& x);

    virtual void procq(const gmVector3& p, double * q);

    MAKE_NAME();
};

#endif

