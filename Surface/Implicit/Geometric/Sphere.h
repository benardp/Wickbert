/**
 * Declaration of the geometric sphere.
 * @file Geometric/Sphere.h
 * @date Fall 2000
 * @author Ed Bachta
 */

#ifndef __ASL_SPHERE_H__
#define __ASL_SPHERE_H__

#include "Geometric.h"

/**
 * This class defines a sphere. Functions are provided for evaluating the
 * implicit surface function at a point in space as well as the gradient 
 * and Hessian.
 */
class Sphere : public Geometric
{
  protected:
    gmVector3 m_x;  ///< Center.
    double m_r;     ///< Radius.
    void init(const gmVector3&,double);
    Box3d m_xlrp(Intervald);
    Intervald m_rlrp(Intervald);

  public:
    // Constructors
    Sphere();
    Sphere(double r);
    Sphere(const gmVector3 & x, double r);

#ifndef INTERVAL_EVAL_ONLY
    virtual double proc(const gmVector3 & x);     // evaluation of function
    virtual gmVector3 grad(const gmVector3 & x);  // evaluation of gradient
    virtual gmMatrix3 hess(const gmVector3 & x);  // evaluation of Hessian
#endif

    virtual Intervald proc(const Box<double>& );
    virtual Box3d grad(const Box<double>& );
    virtual IMatrix3d hess(const Box<double>& );

    virtual Intervald proct(const Box<double>& );
    virtual Box3d gradt(const Box<double>& );

    virtual void procq(const gmVector3&, double*);
    virtual void getq(double*);
    virtual void _setq(double*);
    virtual unsigned int qlen() { return 4; }

    virtual void getqname(char** qn);

	virtual double area() { return 4.0*gmPI*m_r*m_r; }

    virtual const char ** getPixmapXPM(const int& size) const;

    MAKE_NAME();
};

#endif

