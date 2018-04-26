/**
 * @file Subtraction.h
 * Represents the subtraction of two implicit surfaces.
 *
 * @author Steve Zelinka
 */

#ifndef _SUBTRACTION_H
#define _SUBTRACTION_H

#include "Intersection.h"
#include "Complement.h"

class Subtraction : public Blend
{
  private:
    /**
     * The intersection which represents the subtraction.  Any subtraction
     * can simply be represented as the intersection with the complement of
     * the function being subtracted.
     */
    Intersection m_f;

    /**
     * Storage for the complement of g.
     */
    Complement m_comp_g;
    
    /// Put all constructor initialization in one place.
    void init(Implicit f, Implicit g, int cont)
    {
      m_comp_g = Complement(g);
      m_f = Intersection(f,m_comp_g,cont);
    }

  public:
    /// Default constructor
    Subtraction() { init(NULL,NULL,2); }

    /**
    * Constructor for a Subtraction.
    *
    * @param f    The function from which we subtract.
    * @param g    The function which is subtracted.
    * @param cont  The continuity desired.
    */
    Subtraction(Implicit f, Implicit g, int cont = 2) { init(f,g,cont); }

#ifndef INTERVAL_EVAL_ONLY
    virtual double proc(const gmVector3 & x)    { return m_f.proc(x); }
    virtual gmVector3 grad(const gmVector3 & x) { return m_f.grad(x); }
    virtual gmMatrix3 hess(const gmVector3 & x) { return m_f.hess(x); }
#endif

    virtual Intervald proc(const Box<double>& x) { return m_f.proc(x); }
    virtual Box3d grad(const Box<double>& x)     { return m_f.grad(x); }
    virtual IMatrix3d hess(const Box<double>& x) { return m_f.hess(x); }

    /**
    * Get the length of the q vector.
    * @return    The length of the q vector.
    */
    virtual int qlen() { return m_f.qlen(); }

    /**
    * Get the q vector.
    * @param q  On return, the current q vector.  Must have at
    *           least this->qlen() elements.
    */
    virtual void getq(double * q) { m_f.getq(q); }

    /**
    * Set the current q vector.
    * @param q The new q vector.  Must have this->qlen() elements.
    */
    virtual void _setq(double * q) { m_f._setq(q); }

    /**
    * Evaluate the derivative of proc with respect to q at a given point.
    * @param p  The point at which to evaluate.
    * @param q  The derivatives with respect to q vector elements,
    *           on return.  Must have this->qlen() elements.
    */
    virtual void procq(const gmVector3& p, double * q) { m_f.procq(p, q); }

    /**
    * Get the parameter names.
    * @param qn  An array of pointers to hold the resulting parameter
    *            names.  Must have this->qlen() elements.
    */
    virtual void getqname(char **qn) { m_f.getqname(qn); }
};

#endif

