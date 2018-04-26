/**
 * Singular Value Decomposition
 * @file svd.h
 * @date July 6, 2001
 * @author Ed Bachta
 * @remarks Based on code from Numerical Recipies
 */

#ifndef __ASL_SVD_H__
#define __ASL_SVD_H__

#include <math.h>
#include "tnt/tnt.h"
#include "tnt/cmat.h"
#include "tnt/vec.h"

// Macros
static double maxarg1,maxarg2; 
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2)) 

static int iminarg1,iminarg2; 
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double sqrarg; 
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

#define CONVERGE_LIMIT 30

/**
 * The SVD class performs singular value decomposition on TNT matricies.
 * It can also solve systems of linear equations using the SVD method.
 * In doing so it eliminates small singular values to improve matrix
 * conditioning.
 * @todo Perhaps the SVD class could keep a local copy of the decomposition 
 *       for the last matrix decomposed and repeatedly use it in a solver
 *       function which only passes the vector to be solved for.
 */
class SVD {

 private:
 
  int limit; ///< Maximum number of iterations for convergance.

  /// Computes sqrt(a^2 + b^2 ).
  double pythag(double a, double b);

  /// Fixes ill-conditioning by removing small singular values.
  void condition(TNT::Vector<double> &w);

 public:

  /**
   * Constructor.
   * @param l Iteration limit for convergance.
   */
  SVD(int l=CONVERGE_LIMIT) { limit = l; }

  /// Computes SVD using double arrays. (unsupported)
  //bool compute(double **a, int m, int n, double w[], double **v);

  /// Computes SVD using TNT Matrices.
  bool compute(TNT::Matrix<double> &A,
               TNT::Vector<double> &W,
               TNT::Matrix<double> &V);

  /// Solves a system of equations using SVD components.
  void solve(TNT::Matrix<double> &U,
             TNT::Vector<double> &W,
             TNT::Matrix<double> &V,
             TNT::Vector<double> b,
             TNT::Vector<double> &x);

  /// Solves a system of equations.
  bool solve(TNT::Matrix<double> A,
             TNT::Vector<double> &x,
             TNT::Vector<double> b);

};

#endif
