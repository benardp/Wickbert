/**
 * A translation adapter.
 * @file Mover.h
 * @date 7 Nov. 2001
 * @author John C. Hart
 */

#ifndef MOVER_H
#define MOVER_H

#include "UnaryOp.h"

class Mover : public UnaryOp
{
  private:
    gmVector3 m_o;   ///< Translation
    Box3d m_olrp(Intervald);

  public:
    Mover();

    virtual gmVector3 dinv(const gmVector3 & x);
    virtual Box3d dinv(Box3d x);
    virtual Box4d dinv(Box4d x);

    virtual gmMatrix3 dinvjac(const gmVector3 & x);
    virtual IMatrix3d dinvjac(const Box<double>& x);

#ifndef INTERVAL_EVAL_ONLY
    virtual double proc(const gmVector3 & x);
    virtual gmVector3 grad(const gmVector3 & x);
    virtual gmMatrix3 hess(const gmVector3 & x);
#endif

    virtual Intervald proc(const Box<double>& );
    virtual Box3d grad(const Box<double>& );
    virtual IMatrix3d hess(const Box<double>& );

    virtual Intervald proct(const Box<double>& );
    virtual Box3d gradt(const Box<double>& );

    virtual void procq(const gmVector3 & x, double *q);
    virtual void getq(double *q);
    virtual void _setq(double *q);
    virtual unsigned int qlen();
    virtual unsigned int plen();

    virtual void getqname(char** qn);

    MAKE_NAME();

  private:
    bool solve();
    bool impdiff();
};

#endif

