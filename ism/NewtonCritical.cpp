/** 
 * This file contains the method definitions for a NewtonCritical object,
 * whose main purpose is to redefine NewtonInput to use the gradient and
 * interval Hessian.
 *
 * @file NewtonCritical.cpp
 * @author John C. Hart
 * @author Terry Fleury (tfleury@uiuc.edu)
 * @date 11 October 2001
 */

#include "NewtonCritical.h"
#include "Surface/gmTNTconvert.h"
#include "GaussJordan.h"

/** 
 * Default constructor.  This method isn't very useful since it sets the
 * surface to be NULL which means that there is no surface for which to find
 * critical points. 
 */
NewtonCritical::NewtonCritical()
{
  setSurface(NULL);
}

/**
 * Constructs an object which finds critical points for an Implicit.  You
 * MUST pass in an Implicit surface or Newton's method won't do anything.
 * @param surf The implicit surface to be used when finding critical points.
 */
NewtonCritical::NewtonCritical(Implicit* surf)
{
  setSurface(surf);
}

/**
 * Changes the surface to search.  This method is called by the constructors
 * to set up the private surface variable.  It is also virtual and should be
 * overridden (along with getSurface()) to keep the surface variable
 * consistent.
 * @param surf The new Implicit surface to search for critical points.
 * @see getSurface()
 */
void NewtonCritical::setSurface(Implicit* surf)
{
  surface = surf;
}

/**
 * Returns the current surface being searched.  This method is virtual and
 * should be overridden (along with setSurface(Implicit*)) to keep the
 * surface variable consistent.
 * @return The Implicit surface being searched for critical points.
 * @see setSurface(Implicit*)
 */
Implicit* NewtonCritical::getSurface()
{
  return surface;
}

/**
 * Overridden NewtonEquation to solve A(Y-Xc) = b for critical points.
 * In the basic finding of critical points, we look for points in the Box
 * where the gradient is zero.  To do this, we need to use the Hessian of
 * the function so Newton's Method knows which way to go at each
 * iteration.
 *
 * As stated in the SIGGRAPH 97 paper "Guaranteeing the Polygonization of an
 * Implicit Surface Polygonization for Interactive Modeling", this comes
 * down to solving the following equation for Y:
 *
 *    Vc^-1 * V(X) * (Y-Xc) + Vc^-1 * grad(Xc) = 0
 *
 * where V(X) is the interval Hessian matrix of the implicit surface, Vc =
 * midpts(V(X)) is a scalar matrix consisting of the midpoints of the
 * intervals of the Hessian, and Xc = midpts(X) is a scalar vector conisting
 * of the midpoints of the interval to be searched.  
 *
 * So in this case, we have A = Vc^-1 * V(X) and b = -Vc^-1 * grad(Xc). 
 * The reason for including the Vc^-1 term is to try to make the system
 * converge faster.
 */
bool NewtonCritical::NewtonEquation(Box<double> &X, IMatrix &A, Box<double> &b)
{
  if (getSurface() == NULL)
    return false;

  TNTPoint gradXc = convert(getSurface()->grad(convert(X.center())));
  IMatrix3d VofX = getSurface()->hess(X);
  TNT::Matrix<double> Vc = VofX.center();
  TNT::Matrix<double> VcInv(Vc.num_rows(),Vc.num_cols());
  bool okay = GaussJordan::Invert(Vc,VcInv);

  if (okay)
    {
      A = IMatrix(VcInv) * VofX;
      b = (-VcInv) * gradXc;
    }
  else // Vc was singular - couldn't invert it - don't use it
    {
      A = VofX;
      b = -gradXc;
    }
  return true;
}

/**
 * Overridden Newton solver method to stop checking the current box for
 * roots during the subdivision phase.  If the gradient of the Box being 
 * searched does not contain 0, then there aren't any roots to be found 
 * in it and we can quit looking in this Box.  
 * @param X The Box region to check for a 0 gradient.
 * @return False if the gradient of the region contains a 0 value,
 *         indicating that we should check the region for a critical point.
 *         True if we CAN quit early since the region doesn't contain a 0
 *         gradient point.
 */
bool NewtonCritical::NewtonSubdivisionBreak(Box<double> & X)
{
  bool retval = false;
  Implicit* surf = getSurface();
  int i;

  if (surf != NULL)
    {
      Box3d X1 = Box3d(X);
      Box3d temp = surf->grad(X1);
      for (i = 0 ; i < X.size(); i++)
        if (temp[i].isNegative() || temp[i].isPositive())
          retval =  true;   // No root in the region - we can quit early
    }

  return retval;
}


