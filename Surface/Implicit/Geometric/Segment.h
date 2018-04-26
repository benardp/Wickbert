/*
 * @file Segment.h
 *
 * This class defines a line segment. Functions are provided for evaluating
 * the implicit surface function at a point in space as well as the gradient
 * and Hessian.
 *
 * @author Ed Bachta
 */

#ifndef __SEGMENT_H__EHB_
#define __SEGMENT_H__EHB_

#include "Geometric.h"

class Segment : public Geometric
{
  private:
    gmVector3 m_a; // endpoint a
    gmVector3 m_b; // endpoint b

  private:
    void init(const gmVector3&,gmVector3);
    Box3d m_alrp(Intervald);
    Box3d m_blrp(Intervald);

  public:
    // Constructors
    Segment();
    Segment(const gmVector3& a, const gmVector3& b);

#ifndef FORCE_NO_INTERVAL_EVAL
    virtual double proc(const gmVector3 & x);    // evaluation of function
    virtual gmVector3 grad(const gmVector3 & x); // evaluation of gradient
    virtual gmMatrix3 hess(const gmVector3 & x); // evaluation of Hessian
#endif

    virtual Intervald proc(const Box<double>& );
    virtual Box3d grad(const Box<double>& );
    virtual IMatrix3d hess(const Box<double>& ); 
    
    virtual Intervald proct(const Box<double>& );
    virtual Box3d gradt(const Box<double>& );

    virtual void getq(double* q);
    virtual void _setq(double *q);
    virtual unsigned int qlen() { return 6; }

    virtual void getqname(char** qn);
    
    MAKE_NAME();
};

#endif

