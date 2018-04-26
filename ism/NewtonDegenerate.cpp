/** 
 * This file contains the method definitions for a NewtonDegenerate object,
 * whose main purpose is to redefine NewtonInput to use the gradient and
 * interval Hessian.
 *
 * @file NewtonDegenerate.cpp
 * @author Terry Fleury (tfleury@uiuc.edu)
 * @date October 18, 2002
 */

#include "NewtonDegenerate.h"
#include "Surface/gmTNTconvert.h"
#include "GaussJordan.h"

/** 
 * Default constructor.  This method isn't very useful since it sets the
 * surface to be NULL which means that there is no surface for which to find
 * critical points. 
 */
NewtonDegenerate::NewtonDegenerate()
{
  setSurface(NULL);
}

/**
 * Constructs an object which finds critical points for an Implicit.  You
 * MUST pass in an Implicit surface or Newton's method won't do anything.
 * @param surf The implicit surface to be used when finding critical points.
 */
NewtonDegenerate::NewtonDegenerate(Implicit* surf)
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
void NewtonDegenerate::setSurface(Implicit* surf)
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
Implicit* NewtonDegenerate::getSurface()
{
  return surface;
}

/**
 * Overridden NewtonEquation to solve A(Y-Xc) = b for degenerate critical
 * points. In the basic finding of degenerate critical points, we look for
 * points in the Box where the function value is zero AND the gradient is
 * zero over a 4D space (x,y,z,t).  Since we are now dealing with a "time"
 * element of the Box vector, we need to slightly alter our concept of A and
 * b in the equation Ax + y = b.  
 *
 * As stated in the 1978 article "Interval Forms of Newtons (sic) Method" 
 * by Eldon Hansen in "Computing" Volume 20, pp 153-163, we are now actually
 * trying to solve two equations with two unknowns.  We want to find
 * proc(x,t)=0 and grad(X,t)=0.  So we have the following:
 *
 * Our "function" is now a combination of two functions.  When searching for
 * critial points, our "function" was grad(), ie. we were looking for X
 * where grad(X) = 0.  Now our function is a "vector" of two functions:
 *
 * @verbatim
 *    Let p = proc(X,t) 
 *
 *                                        / dp \
 *                                        | -- |
 *                                        | dx |
 *                                        |    |
 *                      | grad(X,t) |     | dp |
 *    New "function" =  | proc(X,t) |  =  | -- |
 *                                        | dy |
 *                                        |    |
 *                                        | dp |
 *                                        | -- |
 *                                        | dy |
 *                                        |    |
 *                                        \ p  /
 *
 *
 *                                                  /  dp    dp    dp    dp  \
 *                                                  | ----  ----  ----  ---- |
 *                                                  | dx^2  dxdy  dxdz  dxdt |
 *                                                  |                        |
 *                     / dgrad(X,t)  dgrad(X,t) \   |  dp    dp    dp    dp  |
 *                     | ----------  ---------- |   | ----  ----  ----  ---- |
 *                     |     dX          dt     |   | dydx  dy^2  dydz  dydt |
 *    New "Jacobian" = |                        | = |                        |
 *                     | dproc(X,t)  dproc(X,t) |   |  dp    dp    dp    dp  |
 *                     | ----------  ---------- |   | ----  ----  ----  ---- |
 *                     \     dX          dt     /   | dzdx  dzdy  dz^2  dzdt |
 *                                                  |                        |
 *                                                  |  dp    dp    dp    dp  |
 *                                                  | ----  ----  ----  ---- |
 *                                                  \  dx    dy    dz    dt  /
 * @endverbatim
 *
 * The upper 3x3 matrix in the 4x4 Jacobian above is simply the Hessian, and
 * the first three elements of the bottom row is simply the grad().  The
 * last column is the tricky part.   Let's start with dp/dt.  
 *
 * When calculating proc(X,t), we use a LRP to linearly interpolate between
 * a set of "old" Qs and "current" Qs.  Eg. Using a simple function like
 * proc=ax, we would have Intervald A = t*(a_current - a_old) + a_old, and
 * thus proc(X,t) = AX = [t*(a_current - a_old) + a_old] * X.  Taking the
 * partial derivative of this wrt t is very simple since everything but the
 * (a_current - a_old) term is a constant.  So we get dproc(X,t)/dt =
 * (a_current - a_old) * X.
 *
 * A similar calculation can be made for dgrad(X,t)/dt.  
 *
 * Practically speaking, we can use procq and gradq (since these two
 * functions don't really care if q is a function of t or not) along with
 * dqdt to calculate dproc/dt and dgrad/dt, using the chain rule. 
 */
bool NewtonDegenerate::NewtonEquation(Box<double> &X,IMatrix &A,Box<double> &b)
{
  Implicit* surf = getSurface();
  int i,j;

  if (getSurface() == NULL)
    return false;

  // First, get the various parts of the "function" vector and construct it.
  TNTPoint gradXc = convert(surf->grad(convert(X.center())));
  TNT::Vector<double> thefunc(4);
  for (i = 0; i < 3; i++)
    thefunc[i] = gradXc[i];
  thefunc[3] = surf->proc(X.center());

  // Then, get the parts of the "Jacobian" matrix and construct it.
  IMatrix3d hessX = surf->hess(X);
  Box3d gradX = surf->grad(X);
  Box3d gradtX = surf->gradt(X);
  IMatrix VofX(4,4);
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      VofX[i][j] = hessX[i][j];
  for (i = 0; i < 3; i++)
    VofX[i][3] = gradtX[i];
  for (j = 0; j < 3; j++)
    VofX[3][j] = gradX[i];
  VofX[3][3] = surf->proct(X);

  TNT::Matrix<double> Vc = VofX.center();
  TNT::Matrix<double> VcInv(Vc.num_rows(),Vc.num_cols());
  bool okay = GaussJordan::Invert(Vc,VcInv);

  if (okay)
    {
      A = IMatrix(VcInv) * VofX;
      b = (-VcInv) * thefunc;
    }
  else // Vc was singular - couldn't invert it - don't use it
    {
      A = VofX;
      b = -thefunc;
    }
  return true;
}

/**
 * Overridden Newton solver method to stop checking the current box for
 * degenerate roots during the subdivision phase.  If the proc and grad of
 * the Box being searched does not contain 0, then there aren't any
 * degenerate roots to be found in it and we can quit looking in this Box.  
 * @param X The Box region to check for a 0 gradient.
 * @return False if the gradient of the region contains a 0 value,
 *         indicating that we should check the region for a critical point.
 *         True if we CAN quit early since the region doesn't contain a 0
 *         gradient point.
 */
bool NewtonDegenerate::NewtonSubdivisionBreak(Box<double> & X)
{
  bool retval = false;
  Implicit* surf = getSurface();

  if (surf != NULL)
    {
      // Since proc is easier to calculate than grad, we first calculate
      // proc and see if the interval contains zero.  If not, then there's
      // no degenerate root in the region and we can stop.  Otherwise, we
      // need to check the grad for the same condition.
      Intervald tempproc = surf->proc((Box4d)X);
      if (tempproc.isNegative() || tempproc.isPositive())
        retval = true;
      else
        { // There's zero proc in there.  Check for a zero grad.
          Box<double> tempgrad = surf->grad((Box4d)X);
          for (int i = 0 ; i < X.size(); i++)
            if (tempgrad[i].isNegative() || tempgrad[i].isPositive())
              retval =  true;   // No root in the region - we can quit early
        }
    }

  return retval;
}

