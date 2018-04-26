/**
 * Declaration of the unary operation.
 * @file Offset.h
 * @date Fall 2000
 * @author Ed Bachta
 */

#ifndef __ASL_OFFSET_H__
#define __ASL_OFFSET_H__

#include "UnaryOp.h"

/**
 * This class defines an offset, which acts as a node in our scene graph. 
 * Any implicit surface can be connected below an Offset node. Functions 
 * are provided for evaluating the implicit surface function at a point in 
 * space as well as the gradient and Hessian.
 *
 * The implicit surface function for an offset can be described as:
 *
 * o(x,y,z) = f(x,y,z) - R
 *
 * where R is the offset value. Because the R term does not depend on 
 * (x,y,z), it will drop out of the gradient, leaving the original surface 
 * gradient. The Hessian will also be the original surface Hessian.
 */
class Offset : public UnaryOp
{
  private:
    double m_r;                          ///< The numerical offset.
    void init(Implicit*,double);
    Intervald m_rlrp(Intervald);
    double m_rdiff(void);

  public:
    Offset();                              ///< Default constructor.
    Offset(Implicit *f, double r);          ///< Explicit constructor.

#ifndef INTERVAL_EVAL_ONLY
    virtual double proc(const gmVector3 & x);       ///< Evaluation of function.
    virtual gmVector3 grad(const gmVector3 & x);   ///< Evaluation of gradient.
    virtual gmMatrix3 hess(const gmVector3 & x);   ///< Evaluation of hessian.
#endif

    virtual Intervald proc(const Box<double>& );
    virtual Box3d grad(const Box<double>& );
    virtual IMatrix3d hess(const Box<double>& );

    virtual Intervald proct(const Box<double>& );
    virtual Box3d gradt(const Box<double>& );

    virtual void procq(const gmVector3&, double*); ///< Evaluation of dFdq.
    virtual void getq(double*);             ///< Retreives parameters.
    virtual void _setq(double*);            ///< Assigns parameters.
    virtual unsigned int qlen();                    ///< Returns the # of parameters.

    virtual void getqname(char** qn);      ///< Retreives parameter names.

    MAKE_NAME();

    virtual const char ** getPixmapXPM(const int& size) const;
};

#endif

