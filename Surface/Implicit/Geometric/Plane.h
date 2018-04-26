/**
 * @file Plane.h
 *
 * This class defines a plane. Functions are provided for evaluating the
 * implicit surface function at a point in space as well as the gradient
 * and Hessian.
 *
 * @author: Ed Bachta
 */

#ifndef __PLANE_H__EHB_
#define __PLANE_H__EHB_

#include "Geometric.h"

class Plane : public Geometric
{
  private:
    gmVector3 m_n; // plane's unit normal
    double m_d;    // shortest distance from origin

  private:
    void init(const gmVector3&,double);
    Box3d m_nlrp(Intervald);
    Intervald m_dlrp(Intervald);

  public:
    // Constructors
    Plane();
    Plane(double a, double b, double c, double d);
    Plane(const gmVector3& n, double d);
    Plane(const gmVector3& n, const gmVector3& p);

#ifndef INTERVAL_EVAL_ONLY
    virtual double proc(const gmVector3 & x);  
    virtual gmVector3 grad(const gmVector3 & x);
    virtual gmMatrix3 hess(const gmVector3 & x); 
#endif
    
    virtual Intervald proc(const Box<double>& x);
    virtual Box3d grad(const Box<double>& x);
    virtual IMatrix3d hess(const Box<double>& x);

    virtual Intervald proct(const Box<double>& x);
    virtual Box3d gradt(const Box<double>& x);

    virtual void procq(const gmVector3&, double*);
    virtual void getq(double*);
    virtual void _setq(double*);
    virtual unsigned int qlen() { return 4; }

    virtual void getqname(char** qn);

    MAKE_NAME();

    virtual const char ** getPixmapXPM(const int& size) const;
};

#endif

