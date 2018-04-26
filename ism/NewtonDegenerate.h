/**
 * This file contains the class declaration for a NewtonDegenerate object,
 * which redefines the NewtonInput method to use the gradient and interval
 * Hessian.
 *
 * @file NewtonDegenerate.h
 * @author Terry Fleury
 * @date October 18, 2002
 */

#ifndef NEWTONDEGENERATE_H
#define NEWTONDEGENERATE_H

#include "Newton.h"

/** 
 * NewtonDegenerate is a class for finding critical points of an Implicit.
 * The class basically redefines the NewtonInput to set the system of
 * equations such that grad(x) == 0 and hess(x)->determinate() == 0.
 */
class NewtonDegenerate : public virtual Newton
{
  protected:
    Implicit* surface;  ///< Find critical points for this Implicit surface.

  public:
    /// Default constructor.
    NewtonDegenerate();
    /// Constructs an object which finds critical points for an Implicit.
    NewtonDegenerate(Implicit*);

    /// Changes the surface to search.
    virtual void setSurface(Implicit* surf);
    /// Returns the current surface being searched.
    virtual Implicit* getSurface();

  protected:
    /// Overridden NewtonEquation to solve A(Y-Xc) = b for critical points.
    virtual bool NewtonEquation(Box<double> &X, IMatrix &A, Box<double> &b);
    /// Overridden Newton solver method to stop subdividing.
    virtual bool NewtonSubdivisionBreak(Box<double> & X);
};

#endif

