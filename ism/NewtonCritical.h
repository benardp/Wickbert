/**
 * This file contains the class declaration for a NewtonCritical object,
 * which redefines the NewtonInput method to use the gradient and interval
 * Hessian.
 *
 * @file NewtonCritical.h
 * @author John C. Hart
 * @author Terry Fleury
 * @date 11 October 2001
 */

#ifndef NEWTONCRITICAL_H
#define NEWTONCRITICAL_H

#include "Newton.h"

/** 
 * NewtonCritical is a class for finding critical points of an Implicit.
 * The class basically redefines the NewtonEquation() method to set up the
 * system of equations for A*(Y-Xc) = b.  See Newton.h for more information
 * on the basics of how it goes about first subdiving space and then calling
 * Newton's Method to find the roots.
 *
 * The equation in NewtonEquation() comes from the SIGGRAPH 97 paper
 * "Guaranteeing the Polygonization of an Implicit Surface Polygonization
 * for Interactive Modeling".  In that paper, we try to find the critical
 * points of an implicit surface.  A critical point is where the gradient of
 * the implicit equation goes to zero.  Skipping the details, we can do this
 * by first using a simple subdivision technique to eliminate most of the
 * empty space (since there really aren't that many critical points for a
 * given implicit).  Once we have a small interval to search, we can switch
 * to Newton's method to find the critical points a little quicker.
 *
 * In the context of the paper, we have the following.  V(X) is the interval
 * Hessian matrix of the implicit surface. Vc = midpts(V(X)) is a scalar
 * matrix consisting of the midpoints of the intervals of the Hessian.  Xc =
 * midpts(X) is a scalar vector conisting of the midpoints of the interval
 * to be searched.  So our equation we are solving is:
 *
 *    Vc^-1 * V(X) * (Y-Xc) + Vc^-1 * grad(Xc) = 0
 *
 * From this we can see that A = Vc^-1 * V(X) and b= Vc^-1 * grad(Xc). 
 * The reason for including the Vc^-1 term is to try to make the system
 * converge faster.
 *
 * For the subdivision part of the searching process, we call
 * NewtonSubdivisionBreak() on a Box to see if we can conclude that there is
 * no root in that Box.  We do this by looking at the gradient of the Box.
 * If this interval vector does not contain 0, then there's no roots to be
 * found.
 */
class NewtonCritical : public virtual Newton
{
  protected:
    Implicit* surface;  ///< Find critical points for this Implicit surface.

  public:
    /// Default constructor.
    NewtonCritical();
    /// Constructs an object which finds critical points for an Implicit.
    NewtonCritical(Implicit*);

    /// Changes the surface to search.
    virtual void setSurface(Implicit*);
    /// Returns the current surface being searched.
    virtual Implicit* getSurface();

  protected:
    /// Overridden NewtonEquation to solve A(Y-Xc) = b for critical points.
    virtual bool NewtonEquation(Box<double> &X, IMatrix &A, Box<double> &b);
    /// Overridden Newton solver method to stop subdividing.
    virtual bool NewtonSubdivisionBreak(Box<double> & X);

  private:
    /// Using the LU decomposition solve for the inverse of a matrix.  
    bool Invert(TNT::Matrix<double> aorig, TNT::Matrix<double> &ainv);

};

#endif

