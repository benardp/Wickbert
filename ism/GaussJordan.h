/** This file contains a class to solve the equation Ax = b where A is an
 * NxN-dimensional matrix, and x and b are n-dimensional vectors.  Note that
 * the methods in this class are all static so that you can simply call them
 * as if they were global functions.
 *
 * @file GaussJordan.h
 * @author Barton Stander
 * @author Terry Fleury
 * @date March 31, 2002
 */

#ifndef GAUSSJORDAN_H
#define GAUSSJORDAN_H

#include "tnt/cmat.h"
#include "tnt/vec.h"

class GaussJordan
{
public:
  /// Solves the equation Ax=b using Gauss-Jordan method.
  static bool Solve(TNT::Matrix<double> &A, TNT::Vector<double> &b);
  
  /// Inverts the matrix A which is returned in AInv.
  static bool Invert(TNT::Matrix<double> &A, TNT::Matrix<double> &AInv);
};

#endif

