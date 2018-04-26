/** 
 * This file contains the method definitions for a Jacobi object.
 * Basically, this object takes in a 3x3 matrix (by calling the 'solve'
 * method), and finds the eigen values and eigen vectors of that matrix.
 * The eigen values/vectors are stored in public members.  This code was
 * originally taken from Numerical Recipies, but Wayne Cochran fixed it up a
 * bit, and Terry Fleury made it an object.  As a result, the eigen vectors
 * are stored in row-major order.
 * @file Jacobi.cpp
 * @author Bart Stander
 * @author Wayne Cochran
 * @author Terry Fleury
 * @date (doxygenated 11 June 2001, -jch)
 * @date (OOPed 03 July 2001 -tgf)
 */

#include <algorithm>
#include <math.h>
#include "Jacobi.h"

// Definition of Jacobi's class constants
const int Jacobi::MAXSWEEP = 50;

/**
 * Called by constructors to initialize the internal data fields.
 * @param mat The matrix to run the Jacobi method on to solve for eigen
 *            vectors and eigen values.
 */
void Jacobi::init(gmMatrix3 mat)
{
  orig = mat;
}

/**
 * Constructor which takes in a gmMatrix3 as the matrix to be solved.
 * @param mat The matrix to run the Jacobi method on to solve for eigen
 *            vectors and eigen values.
 */
Jacobi::Jacobi(gmMatrix3 mat)
{
  init(mat);
}

/**
 * Constructor which takes in a double[3][3] as the matrix to be solved.
 * @param mat The matrix to run the Jacobi method on to solve for eigen
 *            vectors and eigen values.
 */
Jacobi::Jacobi(double mat[3][3])
{
  init(gmMatrix3(mat[0][0],mat[0][1],mat[0][2],
                 mat[1][0],mat[1][1],mat[1][2],
                 mat[2][0],mat[2][1],mat[2][2]));
}

/**
 * Constructor which takes in a TNT::Matrix as the matrix to be solved.
 * @param mat The matrix to run the Jacobi method on to solve for eigen
 *            vectors and eigen values.
 */
Jacobi::Jacobi(TNT::Matrix<double> mat)
{
  init(gmMatrix3(mat[0][0],mat[0][1],mat[0][2],
                 mat[1][0],mat[1][1],mat[1][2],
                 mat[2][0],mat[2][1],mat[2][2]));
}

/**
 * Uses the input matrix to solve for eigen vectors and values.  Given a
 * symetric 3x3 matrix (passed in during object construction), determine its
 * eigenvalues (saved in values) and the eigenvectors (saved in vectors).
 * The number of Jacobi rotations used is returned.  See section 11.1 of
 * Numerical Recipes in Pascal for more information.
 * @return The number of Jacobi rotations used in solving.
 */
int Jacobi::solve()
{
  gmMatrix3 inp;
  int rot = 0;
  register int i;
  gmVector3 B;
  gmVector3 Z;

  // Save the input matrix in orig, use new matrix inp
  inp = orig;
  // Set vectors to the identity matrix
  vectors = gmMatrix3::identity();
  // Set B and values to the diagonal of the input matrix 
  for (i = 0; i < 3; i++)
    B[i] = values[i] = inp[i][i];

  // Rotate until off diagonal elements of input matrix are zero
  for (int sweep = 0; sweep++ < MAXSWEEP;)
    {
      double sum = fabs(inp[0][1]) + fabs(inp[0][2]) + fabs(inp[1][2]);
      double thresh;
      register int p, q;

      if (fabs(sum) < gmEPSILON)  // Normal return relies on quadratic
        {                         // convergence and machine underflow.
          colToRow();  // Put vectors in row-major order
          return rot;  // Return number of rotations
        }

      thresh = (sweep < 4) ? sum * 0.2 / 9.0 : 0.0;  // First three sweeps?

      for (p =  0; p < 2; p++)
        for (q = p+1; q < 3; q++)
          {
            register double g = 100.0 * fabs(inp[p][q]);

            // After 4 sweeps, skip the rotation if the 
            // off-diagonal element is small.
            if ((sweep > 4) && (fabs(g) < gmEPSILON))
              inp[p][q] = 0.0;
            else if (fabs(inp[p][q]) > thresh)
              {
                register double h = values[q] - values[p];
                double c, s, t;  // cosine, sine, tangent of rotation angle
                double tau;
                register int j;

                if (fabs(g) < gmEPSILON)
                  t = inp[p][q] / h;
                else 
                  {
                    double theta = 0.5 * h / inp[p][q];
                    t = 1.0 / (fabs(theta) + sqrt(1.0 + theta*theta));
                    if (theta < 0.0) 
                      t = -t;
                  }

                c = 1.0 / sqrt(1.0 + t*t);   // cosine of rotation angle
                s = t*c;                     // sine of rotation angle
                tau = s / (1.0 + c);

                h = t * inp[p][q];
                Z[p] -= h;
                Z[q] += h;
                values[p] -= h;
                values[q] += h;
                inp[p][q] = 0.0;

                // case of rotations 0 <= j < p-1 
                for (j = 0; j <= p-1; j++) 
                  {
                    g = inp[j][p];
                    h = inp[j][q];
                    inp[j][p] = g - s*(h + g*tau);
                    inp[j][q] = h + s*(g - h*tau);
                  }

                // case of rotations p < j < q 
                for (j = p + 1;  j < q; j++) 
                  {
                    g = inp[p][j];
                    h = inp[j][q];
                    inp[p][j] = g - s*(h - g*tau);
                    inp[j][q] = h + s*(g - h*tau);
                  }

                // case of rotations q < j < 3 
                for (j = q + 1; j < 3; j++) 
                  {
                    g = inp[p][j];
                    h = inp[q][j];
                    inp[p][j] = g - s*(h + g*tau);
                    inp[q][j] = h + s*(g - h*tau);
                  }

                // Set the eigen vectors
                for (j = 0; j < 3; j++) 
                  {
                    g = vectors[j][p];
                    h = vectors[j][q];
                    vectors[j][p] = g - s*(h + g*tau);
                    vectors[j][q] = h + s*(g - h*tau);
                  }
                rot++;
              } // end if (fabs(inp[p][q]) > thresh)
          } // end 'for p' and 'for q' loops

      // Set the eigen values
      B += Z;
      values = B;
      Z.assign(0.0,0.0,0.0);
    } // end 'for sweep' loop

  colToRow(); // Put vectors in row-major order
  return -1;  // Non-normal return - too many rotations
}

/**
 * Sorts the eigen vectors based on eigen values, in ascending order.  This
 * is basically an expanded bubblesort for 3 values.  The resulting eigen
 * values/vectors will be sorted such that values[0]/vectors[0] has the
 * smallest eigen value, and values[2]/vectors[2] has the largest eigen
 * value.
 */
void Jacobi::sort()
{
  int i;

  if (values[0] > values[1])
    {
      std::swap(values[0],values[1]);
      for (i = 0; i < 3; i++)
        std::swap(vectors[0][i],vectors[1][i]);
    }
  if (values[1] > values[2])
    {
      std::swap(values[1],values[2]);
      for (i = 0; i < 3; i++)
        std::swap(vectors[1][i],vectors[2][i]);
    }
  if (values[0] > values[1])
    {
      std::swap(values[0],values[1]);
      for (i = 0; i < 3; i++)
        std::swap(vectors[0][i],vectors[1][i]);
    }
}

/**
 * Changes the matrix from column-major order to row-major order.  This
 * method is called from within solve just before the method returns to put
 * the eigen vectors in row-major order.
 */
void Jacobi::colToRow()
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      std::swap(vectors[i][j],vectors[j][i]);
}

