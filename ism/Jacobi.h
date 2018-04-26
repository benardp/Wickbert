/** 
 * This code was originally from Numerical Recipies, but Wayne Cochran fixed
 * it up a bit, and Terry Fleury then made it more object oriented to contain
 * the algorithm in the Jacobi class.
 * IMPORTANT: the original jacobi code returned the eigen vectors in COLUMN
 * MAJOR order, but extra code was added to swap the components.  So now the
 * eigen vectors are returned in the more standard ROW MAJOR order.  
 * @file Jacobi.h
 * @author Bart Stander 
 * @author Wayne Cochran
 * @author Terry Fleury
 * @date (doxygenated 11 June 2001, -jch)
 * @date (OOPed 03 July 2001 -TGF)
 */

#ifndef JACOBI_H
#define JACOBI_H

#include "libgm/gm.h"
#include "tnt/cmat.h"
#include "tnt/vec.h"

/**
 * The Jacobi class contains a method to solve for the eigen vectors and
 * eigen values of a 3x3 matrix.  To use this class, you call one of the
 * constructors using either a gmMatrix3, double[3][3], or
 * TNT::Matrix<double> as the input matrix to solve.  Next, call the solve()
 * method to actually solve for eigen vectors/values.  The resulting eigen
 * vectors are stored in the gmMatrix3 member 'vectors', and the eigen
 * values are stored in the gmVector3 member 'values'.  The original matrix
 * is preserved in the gmMatrix3 member 'orig' just in case you need that
 * information after solving. 
 *
 * Note that the resulting eigen vectors are stored in row-major order for
 * easy access.  (This was not the case with the original jacobi code which
 * required that you swap elements along the diagonal before use).  Also,
 * the resulting eigen vectors/values are NOT sorted by default.  If you
 * want your eigen vectors/values sorted in ascending order, you must call 
 * sort() after solve().
 */
class Jacobi
{
  private:
    /// Called by constructors to initialize the internal data fields.
    void init(gmMatrix3);
    /// Changes the matrix from column-major order to row-major order.
    void colToRow();     

  public:
    static const int MAXSWEEP; ///< Max # of times to sweep through the matrix
    // Note: the above const is defined at the top of Jacobi.cpp
    gmMatrix3 orig;     ///< The original input matrix to solve.
    gmMatrix3 vectors;  ///< The three eigen vectors solved.
    gmVector3 values;   ///< The three eigen values solved.
    
    /// Constructor which takes in a gmMatrix3 as the matrix to be solved.
    Jacobi(gmMatrix3);
    /// Constructor which takes in a double[3][3] as the matrix to be solved.
    Jacobi(double[3][3]);
    /// Constructor with takes in a TNT::Matrix as the matrix to be solved.
    Jacobi(TNT::Matrix<double>);

    /// Uses the input matrix to solve for eigen vectors and values.
    int solve();
    /// Sorts the eigen vectors based on eigen values, in ascending order.
    void sort();
};

#endif /* JACOBI_H */

