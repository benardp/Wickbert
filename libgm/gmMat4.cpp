/**
 * This file contains method definitions for the 4x4 matrix class.
 * @file  gmMat4.cpp
 * @date  15 June 1994
 * @author Ferdi Scheepers
 * @author Stephen F May
 */

#include "gmMat4.h"

#include "gmVec3.h"
#include "gmVec4.h"

// private function: RCD
// - dot product of row i of matrix A and row j of matrix B

/**
 * This (private) utility function returns the dot product of 'row i of
 * matrix A' and 'column j of matrix B'.  It is used when multiplying two
 * matrices together.
 * @param A The first matrix.
 * @param B The second matrix.
 * @param i The row of matrix A to use.
 * @param j The column of matrix B to use.
 * @return The dot product of 'matrix A, row i' and 'matrix B, column j'.
 */
inline double RCD(const gmMatrix4& A, const gmMatrix4& B, int i, int j)
{
  return A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j] +
         A[i][3] * B[3][j];
}

// private function: MINOR
// - MINOR of M[r][c]; r in {0,1,2,3}-{r0,r1,r2}; c in {0,1,2,3}-{c0,c1,c2}

/**
 * This (private) utility function calculate the minor product difference of
 * a matrix.  It is used when calculating the adjoint and determinant of a
 * matrix.   
 * @param M The matrix to use.
 * @param r0 The first row.
 * @param r1 The second row.
 * @param r2 The third row.
 * @param c0 The first column.
 * @param c1 The second column.
 * @param c2 The third column.
 * @return The MINOR of the matrix
 * @remarks Each of the values for r0, r1, r2, c0, c1, and c2 must lie
 * within [0,3] since no error checking is performed.
 */
inline double MINOR(const gmMatrix4& M,
                    int r0, int r1, int r2, int c0, int c1, int c2)
{
  return M[r0][c0] * (M[r1][c1] * M[r2][c2] - M[r2][c1] * M[r1][c2]) -
	 M[r0][c1] * (M[r1][c0] * M[r2][c2] - M[r2][c0] * M[r1][c2]) +
	 M[r0][c2] * (M[r1][c0] * M[r2][c1] - M[r2][c0] * M[r1][c1]);
}

// CONSTRUCTORS

/** 
 * Default constructor. Sets all elements of the matrix to be 0.0.
 */
gmMatrix4::gmMatrix4()
{
  assign(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0);
}

/** 
 * Copy constructor.  Creates a new matrix by using the values of another
 * matrix.
 * @param M The matrix to copy.
 */
gmMatrix4::gmMatrix4(const gmMatrix4& M)
{
  assign(M[0][0], M[0][1], M[0][2], M[0][3],
	 M[1][0], M[1][1], M[1][2], M[1][3],
	 M[2][0], M[2][1], M[2][2], M[2][3],
	 M[3][0], M[3][1], M[3][2], M[3][3]);
}

/**
 * Constructs a new matrix object using sixteen doubles.  
 * @param a00 The value for element in first row, first column.
 * @param a01 The value for element in first row, second column.
 * @param a02 The value for element in first row, third column.
 * @param a03 The value for element in first row, fourth column.
 * @param a10 The value for element in second row, first column.
 * @param a11 The value for element in second row, second column.
 * @param a12 The value for element in second row, third column.
 * @param a13 The value for element in second row, fourth column.
 * @param a20 The value for element in third row, first column.
 * @param a21 The value for element in third row, second column.
 * @param a22 The value for element in third row, third column.
 * @param a23 The value for element in third row, fourth column.
 * @param a30 The value for element in fourth row, first column.
 * @param a31 The value for element in fourth row, second column.
 * @param a32 The value for element in fourth row, third column.
 * @param a33 The value for element in fourth row, fourth column.
 */
gmMatrix4::gmMatrix4(double a00, double a01, double a02, double a03,
		     double a10, double a11, double a12, double a13,
		     double a20, double a21, double a22, double a23,
		     double a30, double a31, double a32, double a33)
{
  assign(a00, a01, a02, a03,
         a10, a11, a12, a13,
	 a20, a21, a22, a23,
	 a30, a31, a32, a33);
}

// ASSIGNMENT

/**
 * Sets the matrix's elements to the given values.
 * @param a00 The value for element in first row, first column.
 * @param a01 The value for element in first row, second column.
 * @param a02 The value for element in first row, third column.
 * @param a03 The value for element in first row, fourth column.
 * @param a10 The value for element in second row, first column.
 * @param a11 The value for element in second row, second column.
 * @param a12 The value for element in second row, third column.
 * @param a13 The value for element in second row, fourth column.
 * @param a20 The value for element in third row, first column.
 * @param a21 The value for element in third row, second column.
 * @param a22 The value for element in third row, third column.
 * @param a23 The value for element in third row, fourth column.
 * @param a30 The value for element in fourth row, first column.
 * @param a31 The value for element in fourth row, second column.
 * @param a32 The value for element in fourth row, third column.
 * @param a33 The value for element in fourth row, fourth column.
 */
gmMatrix4& gmMatrix4::assign(double a00, double a01, double a02, double a03,
		             double a10, double a11, double a12, double a13,
		             double a20, double a21, double a22, double a23,
		             double a30, double a31, double a32, double a33)
{
  m_[0][0] = a00; m_[0][1] = a01; m_[0][2] = a02; m_[0][3] = a03;
  m_[1][0] = a10; m_[1][1] = a11; m_[1][2] = a12; m_[1][3] = a13;
  m_[2][0] = a20; m_[2][1] = a21; m_[2][2] = a22; m_[2][3] = a23;
  m_[3][0] = a30; m_[3][1] = a31; m_[3][2] = a32; m_[3][3] = a33;
  return *this;
}

/**
 * Sets the value of the matrix to another matrix.
 * @param M The other matrix to use for setting values.
 * @return A matrix with new values.
 */
gmMatrix4& gmMatrix4::operator =(const gmMatrix4& M)
{
  assign(M[0][0], M[0][1], M[0][2], M[0][3],
	 M[1][0], M[1][1], M[1][2], M[1][3],
	 M[2][0], M[2][1], M[2][2], M[2][3],
	 M[3][0], M[3][1], M[3][2], M[3][3]);
  return *this;
}

// MATH OPERATORS

/**
 * Addition/assignment operator.
 * @param M The matrix to add to the current matrix.
 * @return The matrix += M.
 */
gmMatrix4& gmMatrix4::operator +=(const gmMatrix4& M)
{
  m_[0][0] += M[0][0]; m_[0][1] += M[0][1];
  m_[0][2] += M[0][2]; m_[0][3] += M[0][3];
  
  m_[1][0] += M[1][0]; m_[1][1] += M[1][1];
  m_[1][2] += M[1][2]; m_[1][3] += M[1][3];
  
  m_[2][0] += M[2][0]; m_[2][1] += M[2][1];
  m_[2][2] += M[2][2]; m_[2][3] += M[2][3];
  
  m_[3][0] += M[3][0]; m_[3][1] += M[3][1];
  m_[3][2] += M[3][2]; m_[3][3] += M[3][3];
  return *this;
}

/**
 * Subtraction/assignment operator.
 * @param M The matrix to subtract from the current matrix.
 * @return The matrix -= M.
 */
gmMatrix4& gmMatrix4::operator -=(const gmMatrix4& M)
{
  m_[0][0] -= M[0][0]; m_[0][1] -= M[0][1];
  m_[0][2] -= M[0][2]; m_[0][3] -= M[0][3];
  
  m_[1][0] -= M[1][0]; m_[1][1] -= M[1][1];
  m_[1][2] -= M[1][2]; m_[1][3] -= M[1][3];
  
  m_[2][0] -= M[2][0]; m_[2][1] -= M[2][1];
  m_[2][2] -= M[2][2]; m_[2][3] -= M[2][3];
  
  m_[3][0] -= M[3][0]; m_[3][1] -= M[3][1];
  m_[3][2] -= M[3][2]; m_[3][3] -= M[3][3];
  return *this;
}

/**
 * Multiplication/assignment operator.
 * @param M The matrix to multiply with the current matrix.
 * @return The matrix *= M.
 */
gmMatrix4& gmMatrix4::operator *=(const gmMatrix4& M)
{
  assign(RCD(*this, M, 0, 0), RCD(*this, M, 0, 1), 
	 RCD(*this, M, 0, 2), RCD(*this, M, 0, 3),
	 RCD(*this, M, 1, 0), RCD(*this, M, 1, 1),
	 RCD(*this, M, 1, 2), RCD(*this, M, 1, 3),
	 RCD(*this, M, 2, 0), RCD(*this, M, 2, 1),
	 RCD(*this, M, 2, 2), RCD(*this, M, 2, 3),
	 RCD(*this, M, 3, 0), RCD(*this, M, 3, 1),
	 RCD(*this, M, 3, 2), RCD(*this, M, 3, 3));
  return *this;
}

/**
 * Multiplication/assignment operator with a scalar.
 * @param d The scalar value multiplied with each element of the matrix.
 * @return d*M.
 */
gmMatrix4& gmMatrix4::operator *=(double d)
{
  m_[0][0] *= d; m_[0][1] *= d; m_[0][2] *= d; m_[0][3] *= d;
  m_[1][0] *= d; m_[1][1] *= d; m_[1][2] *= d; m_[1][3] *= d;
  m_[2][0] *= d; m_[2][1] *= d; m_[2][2] *= d; m_[2][3] *= d;
  m_[3][0] *= d; m_[3][1] *= d; m_[3][2] *= d; m_[3][3] *= d;
  return *this;
}

/**
 * Division/assignment operator with a scalar.
 * @param d The scalar value divided into each element of the matrix.
 * @return 1/d * M.
 */
gmMatrix4& gmMatrix4::operator /=(double d)
{
  double di = 1 / d;
  m_[0][0] *= di; m_[0][1] *= di; m_[0][2] *= di; m_[0][3] *= di;
  m_[1][0] *= di; m_[1][1] *= di; m_[1][2] *= di; m_[1][3] *= di;
  m_[2][0] *= di; m_[2][1] *= di; m_[2][2] *= di; m_[2][3] *= di;
  m_[3][0] *= di; m_[3][1] *= di; m_[3][2] *= di; m_[3][3] *= di;
  return *this;
}

/**
 * Addition operator.  Adds two matrices together.
 * @param M The matrix to add to the current matrix.
 * @return A new matrix resulting from the addition of two matrices.
 */
gmMatrix4 gmMatrix4::operator +(const gmMatrix4& M) const
{
  return gmMatrix4(m_[0][0] + M[0][0], m_[0][1] + M[0][1],
		   m_[0][2] + M[0][2], m_[0][3] + M[0][3],
		   m_[1][0] + M[1][0], m_[1][1] + M[1][1],
		   m_[1][2] + M[1][2], m_[1][3] + M[1][3],
		   m_[2][0] + M[2][0], m_[2][1] + M[2][1],
		   m_[2][2] + M[2][2], m_[2][3] + M[2][3],
		   m_[3][0] + M[3][0], m_[3][1] + M[3][1],
		   m_[3][2] + M[3][2], m_[3][3] + M[3][3]);
}

/**
 * Subtraction operator.  Subtracts one matrix from another.
 * @param M The matrix to subtract from the current matrix.
 * @return A new matrix resulting from the difference of two matrices.
 */
gmMatrix4 gmMatrix4::operator -(const gmMatrix4& M) const
{
  return gmMatrix4(m_[0][0] - M[0][0], m_[0][1] - M[0][1],
		   m_[0][2] - M[0][2], m_[0][3] - M[0][3],
		   m_[1][0] - M[1][0], m_[1][1] - M[1][1],
		   m_[1][2] - M[1][2], m_[1][3] - M[1][3],
		   m_[2][0] - M[2][0], m_[2][1] - M[2][1],
		   m_[2][2] - M[2][2], m_[2][3] - M[2][3],
		   m_[3][0] - M[3][0], m_[3][1] - M[3][1],
		   m_[3][2] - M[3][2], m_[3][3] - M[3][3]);
}

/**
 * Unary negation operator. 
 * @return A new matrix with all elements negated.
 */
gmMatrix4 gmMatrix4::operator -() const
{
  return gmMatrix4(-m_[0][0], -m_[0][1], -m_[0][2], -m_[0][3],
		   -m_[1][0], -m_[1][1], -m_[1][2], -m_[1][3],
		   -m_[2][0], -m_[2][1], -m_[2][2], -m_[2][3],
		   -m_[3][0], -m_[3][1], -m_[3][2], -m_[3][3]);
}

/**
 * Multiplication operator.  Multiplies two matrices together.
 * @param M The matrix to multiply with the current matrix.
 * @return A new matrix resulting from the product of two matrices.
 */
gmMatrix4 gmMatrix4::operator *(const gmMatrix4& M) const
{
  return gmMatrix4(RCD(*this, M, 0, 0), RCD(*this, M, 0, 1), 
		   RCD(*this, M, 0, 2), RCD(*this, M, 0, 3),
		   RCD(*this, M, 1, 0), RCD(*this, M, 1, 1),
		   RCD(*this, M, 1, 2), RCD(*this, M, 1, 3),
		   RCD(*this, M, 2, 0), RCD(*this, M, 2, 1),
		   RCD(*this, M, 2, 2), RCD(*this, M, 2, 3),
		   RCD(*this, M, 3, 0), RCD(*this, M, 3, 1),
		   RCD(*this, M, 3, 2), RCD(*this, M, 3, 3));
}

/**
 * Scalar multiplication operator.  Multiplies a scalar with the matrix.
 * @param d The scalar value multiplied with each element of the matrix.
 * @return A new matrix with d*M.
 */
gmMatrix4 gmMatrix4::operator *(double d) const
{
  return gmMatrix4(m_[0][0] * d, m_[0][1] * d, m_[0][2] * d, m_[0][3] * d,
		   m_[1][0] * d, m_[1][1] * d, m_[1][2] * d, m_[1][3] * d,
		   m_[2][0] * d, m_[2][1] * d, m_[2][2] * d, m_[2][3] * d,
		   m_[3][0] * d, m_[3][1] * d, m_[3][2] * d, m_[3][3] * d);
}

/**
 * Scalar division operator.  Divides a scalar into the matrix.
 * @param d The scalar value divided into each element of the matrix.
 * @return A new matrix with 1/d * M.
 */
gmMatrix4 gmMatrix4::operator /(double d) const
{
  assert(!gmIsZero(d));
  double di = 1 / d;
  return gmMatrix4(m_[0][0] * di, m_[0][1] * di, m_[0][2] * di, m_[0][3] * di,
		   m_[1][0] * di, m_[1][1] * di, m_[1][2] * di, m_[1][3] * di,
		   m_[2][0] * di, m_[2][1] * di, m_[2][2] * di, m_[2][3] * di,
		   m_[3][0] * di, m_[3][1] * di, m_[3][2] * di, m_[3][3] * di);
}

/**
 * Scalar multiplication operator.  Multiplies a scalar with the matrix.
 * @param d The scalar value multiplied with each element of the matrix.
 * @param M The matrix to multiply.
 * @return A new matrix with d*M.
 */
gmMatrix4 operator *(double d, const gmMatrix4& M)
{
  return gmMatrix4(M[0][0] * d, M[0][1] * d, M[0][2] * d, M[0][3] * d,
		   M[1][0] * d, M[1][1] * d, M[1][2] * d, M[1][3] * d,
		   M[2][0] * d, M[2][1] * d, M[2][2] * d, M[2][3] * d,
		   M[3][0] * d, M[3][1] * d, M[3][2] * d, M[3][3] * d);
}

/**
 * Vector multiplication operator.  Multiplies a vector with the matrix.
 * @param v The vector to multiply.
 * @return A new vector resulting from M*v.
 */
gmVector4 gmMatrix4::operator *(const gmVector4& v) const
{
  return gmVector4(
    m_[0][0] * v[0] + m_[0][1] * v[1] + m_[0][2] * v[2] + m_[0][3] * v[3],
    m_[1][0] * v[0] + m_[1][1] * v[1] + m_[1][2] * v[2] + m_[1][3] * v[3],
    m_[2][0] * v[0] + m_[2][1] * v[1] + m_[2][2] * v[2] + m_[2][3] * v[3],
    m_[3][0] * v[0] + m_[3][1] * v[1] + m_[3][2] * v[2] + m_[3][3] * v[3]);
}

/**
 * Vector multiplication operator.  Multiplies a vector with the matrix.
 * @param v The vector to multiply.
 * @param M The matrix to multiply.
 * @return A new vector resulting from M*v.
 */
gmVector4 operator *(const gmVector4& v, const gmMatrix4& M)
{
  return gmVector4(
    v[0] * M[0][0] + v[1] * M[1][0] + v[2] * M[2][0] + v[3] * M[3][0],
    v[0] * M[0][1] + v[1] * M[1][1] + v[2] * M[2][1] + v[3] * M[3][1],
    v[0] * M[0][2] + v[1] * M[1][2] + v[2] * M[2][2] + v[3] * M[3][2],
    v[0] * M[0][3] + v[1] * M[1][3] + v[2] * M[2][3] + v[3] * M[3][3]);
}

/**
 * Equivalence operator.  Checks to see if two matrices are equal.
 * @param v The matrix to check for equivalence.
 * @return True if all elements are fuzzy equal, false otherwise.
 */
bool gmMatrix4::operator ==(const gmMatrix4& M) const
{
  return (gmFuzEQ(m_[0][0], M[0][0]) && gmFuzEQ(m_[0][1], M[0][1]) &&
	  gmFuzEQ(m_[0][2], M[0][2]) && gmFuzEQ(m_[0][3], M[0][3]) &&
	 
	  gmFuzEQ(m_[1][0], M[1][0]) && gmFuzEQ(m_[1][1], M[1][1]) &&
	  gmFuzEQ(m_[1][2], M[1][2]) && gmFuzEQ(m_[1][3], M[1][3]) &&

 	  gmFuzEQ(m_[2][0], M[2][0]) && gmFuzEQ(m_[2][1], M[2][1]) &&
	  gmFuzEQ(m_[2][2], M[2][2]) && gmFuzEQ(m_[2][3], M[2][3]) &&

	  gmFuzEQ(m_[3][0], M[3][0]) && gmFuzEQ(m_[3][1], M[3][1]) &&
	  gmFuzEQ(m_[3][2], M[3][2]) && gmFuzEQ(m_[3][3], M[3][3]));
}

/**
 * Non-equivalence operator.  Checks to see if two matrices are not equal.
 * @param v The matrix to check for non-equivalence.
 * @return False if all elements are fuzzy equal, true otherwise.
 */
bool gmMatrix4::operator !=(const gmMatrix4& M) const
{
  return (!(*this == M));
}

// OPERATIONS

/**
 * Returns the inverse of the matrix.  This method makes sure that the
 * matrix can be inverted by verifying that the matrix is not singular.
 * @return A new matrix which is the inverse of the original.
 */
gmMatrix4 gmMatrix4::inverse() const
{
  assert(!isSingular());
  return adjoint() * gmInv(determinant());
}

/**
 * Returns the transpose of the matrix.  This simply swaps elements across
 * the diagonal.
 * @return A new matrix with elements transposed.
 */
gmMatrix4 gmMatrix4::transpose() const
{
  return gmMatrix4(m_[0][0], m_[1][0], m_[2][0], m_[3][0],
		   m_[0][1], m_[1][1], m_[2][1], m_[3][1],
		   m_[0][2], m_[1][2], m_[2][2], m_[3][2],
		   m_[0][3], m_[1][3], m_[2][3], m_[3][3]);
}

/**
 * Returns the adjoint of the matrix.
 * @return The adjoint of the matrix.
 */
gmMatrix4 gmMatrix4::adjoint() const
{
  gmMatrix4 A;
  
  A[0][0] =  MINOR(*this, 1, 2, 3, 1, 2, 3);
  A[0][1] = -MINOR(*this, 0, 2, 3, 1, 2, 3);
  A[0][2] =  MINOR(*this, 0, 1, 3, 1, 2, 3);
  A[0][3] = -MINOR(*this, 0, 1, 2, 1, 2, 3);
  A[1][0] = -MINOR(*this, 1, 2, 3, 0, 2, 3);
  A[1][1] =  MINOR(*this, 0, 2, 3, 0, 2, 3);
  A[1][2] = -MINOR(*this, 0, 1, 3, 0, 2, 3);
  A[1][3] =  MINOR(*this, 0, 1, 2, 0, 2, 3);
  A[2][0] =  MINOR(*this, 1, 2, 3, 0, 1, 3);
  A[2][1] = -MINOR(*this, 0, 2, 3, 0, 1, 3);
  A[2][2] =  MINOR(*this, 0, 1, 3, 0, 1, 3);
  A[2][3] = -MINOR(*this, 0, 1, 2, 0, 1, 3);
  A[3][0] = -MINOR(*this, 1, 2, 3, 0, 1, 2);
  A[3][1] =  MINOR(*this, 0, 2, 3, 0, 1, 2);
  A[3][2] = -MINOR(*this, 0, 1, 3, 0, 1, 2);
  A[3][3] =  MINOR(*this, 0, 1, 2, 0, 1, 2);

  return A;
}

/**
 * Returns the determinant of the matrix.
 * @return The determinant of the matrix.
 */
double gmMatrix4::determinant() const
{
  return m_[0][0] * MINOR(*this, 1, 2, 3, 1, 2, 3) -
	 m_[0][1] * MINOR(*this, 1, 2, 3, 0, 2, 3) +
	 m_[0][2] * MINOR(*this, 1, 2, 3, 0, 1, 3) -
	 m_[0][3] * MINOR(*this, 1, 2, 3, 0, 1, 2);
}

/**
 * Returns if the matrix is singular or not.  This method checks to see if
 * the determinant of the matrix equals 0.0.
 * @return True if the matrix is singular (determinant == 0.0), false
 * otherwise.
 */
bool gmMatrix4::isSingular() const
{
  return gmIsZero(determinant());
}

/**
 * Transforms (multiplies) a gmVector3 by the current matrix.
 * @param v The gmVector3 to multiply.
 * @return A new vector reflecting v * the first three columns of M.
 */
gmVector3 gmMatrix4::transform(const gmVector3& v) const
{
  return gmVector3(
    v[0] * m_[0][0] + v[1] * m_[1][0] + v[2] * m_[2][0] + m_[3][0],
    v[0] * m_[0][1] + v[1] * m_[1][1] + v[2] * m_[2][1] + m_[3][1],
    v[0] * m_[0][2] + v[1] * m_[1][2] + v[2] * m_[2][2] + m_[3][2]);
}

// TRANSFORMATION MATRICES

/**
 * Returns the identity matrix.
 * @return A new matrix with ones along the diagonal and zeros elsewhere.
 */
gmMatrix4 gmMatrix4::identity()
{
  return gmMatrix4(1, 0, 0, 0,
		   0, 1, 0, 0,
		   0, 0, 1, 0,
		   0, 0, 0, 1);
}

/**
 * Returns a matrix corresponding to degrees rotation about an axis.
 * @param angle The angle in degrees to rotate in the xy-plane.
 * @param axis A vector representing the axis of rotation.
 * @return A new rotation matrix.
 */
gmMatrix4 gmMatrix4::rotate(double angle, const gmVector3& axis)
{
  gmMatrix4 m_;
  double length = axis.length();
  double a = axis[0] / length;
  double b = axis[1] / length;
  double c = axis[2] / length;
  double aa = a * a;
  double bb = b * b;
  double cc = c * c;
  double sine = sin(gmRadians(angle));
  double cosine = cos(gmRadians(angle));
  double omcos = 1 - cosine;

  m_[0][0] = aa + (1 - aa) * cosine;
  m_[1][1] = bb + (1 - bb) * cosine;
  m_[2][2] = cc + (1 - cc) * cosine;
  m_[0][1] = a * b * omcos + c * sine;
  m_[0][2] = a * c * omcos - b * sine;
  m_[1][0] = a * b * omcos - c * sine;
  m_[1][2] = b * c * omcos + a * sine;
  m_[2][0] = a * c * omcos + b * sine;
  m_[2][1] = b * c * omcos - a * sine;
  m_[0][3] = m_[1][3] = m_[2][3] = m_[3][0] = m_[3][1] = m_[3][2] = 0;
  m_[3][3] = 1;
   
  return m_;
}

/**
 * Returns a matrix corresponding to scaling by (x,y,z) factors.
 * @param x Scaling factor for x direction.
 * @param y Scaling factor for y direction.
 * @param z Scaling factor for y direction.
 * @return A new scaling matrix.
 */
gmMatrix4 gmMatrix4::scale(double x, double y, double z)
{
  return gmMatrix4(x, 0, 0, 0,
		   0, y, 0, 0,
		   0, 0, z, 0,
		   0, 0, 0, 1);
}

/**
 * Returns a matrix corresponding to translation by (x,y).
 * @param x Amount to translate in x direction.
 * @param y Amount to translate in y direction.
 * @param z Amount to translate in z direction.
 * @return A new translation matrix.
 * @remark Note that this assumes you are working in 3-D since this returns
 * a matrix using homogeneous coordinate space.
 */
gmMatrix4 gmMatrix4::translate(double x, double y, double z)
{
  return gmMatrix4(1, 0, 0, 0,
		   0, 1, 0, 0,
		   0, 0, 1, 0,
		   x, y, z, 1);
}

// CUBIC BASIS MATRICES

/**
 * Returns a Bezier spline basis matrix.
 * @return A Bezier spline basis matrix.
 */
gmMatrix4 gmMatrix4::bezierBasis()
{
  return gmMatrix4(-1,  3, -3,  1,
	 	    3, -6,  3,  0,
		   -3,  3,  0,  0,
		    1,  0,  0,  0);
}

/**
 * Returns a b-spline basis matrix.
 * @return A b-spline basis matrix.
 */
gmMatrix4 gmMatrix4::bsplineBasis()
{  
  return gmMatrix4(-1,  3, -3,  1,
		    3, -6,  3,  0,
		   -3,  0,  3,  0,
		    1,  4,  1,  0) / 6;
}

/**
 * Returns a Catmull-Rom spline basis matrix.
 * @return A Catmull-Rom spline basis matrix.
 */
gmMatrix4 gmMatrix4::catmullromBasis()
{
  return gmMatrix4(-1,  3, -3,  1,
		    2, -5,  4, -1,
		   -1,  0,  1,  0,
		    0,  2,  0,  0) / 2;
}

/**
 * Returns a Hermite spline basis matrix.
 * @return A Hermite spline basis matrix.
 */
gmMatrix4 gmMatrix4::hermiteBasis()
{
  return gmMatrix4( 2,  1, -2,  1,
		   -3, -2,  3, -1,
		    0,  1,  0,  0,
		    1,  0,  0,  0);
}

/**
 * Returns a Tensed b-spline basis matrix.
 * @param tension The tension parameter for the b-spline.
 * @return A Tensed b-spline basis matrix.
 */
gmMatrix4 gmMatrix4::tensedBSplineBasis(double tension)
{
  gmMatrix4 m;
  double sixth = 1.0 / 6.0;

  m[0][0] = sixth * (-tension);
  m[0][1] = sixth * (12.0 - 9.0 * tension);
  m[0][2] = sixth * (9.0 * tension - 12.0);
  m[0][3] = sixth * tension;

  m[1][0] = sixth * 3.0 * tension;
  m[1][1] = sixth * (12.0 * tension - 18.0);
  m[1][2] = sixth * (18.0 - 15.0 * tension);
  m[1][3] = 0.0;

  m[2][0] = sixth * -3.0 * tension;
  m[2][1] = 0.0;
  m[2][2] = sixth * 3.0 * tension;
  m[2][3] = 0.0;

  m[3][0] = sixth * tension;
  m[3][1] = sixth * (6.0 - 2.0 * tension);
  m[3][2] = sixth * tension;
  m[3][3] = 0.0;

  return m;
}

/**
 * Returns a Cardinal spline basis matrix.
 * @param tension The tension parameter for the spline.
 * @return A Cardinal spline basis matrix.
 */
gmMatrix4 gmMatrix4::cardinalBasis(double tension)
{
  gmMatrix4 m;
  
  m[0][0] = -tension;
  m[0][1] = 2.0 - tension;
  m[0][2] = tension - 2.0;
  m[0][3] = tension;

  m[1][0] = 2.0 * tension;
  m[1][1] = tension - 3.0;
  m[1][2] = 3 - 2.0 * tension;
  m[1][3] = -tension;

  m[2][0] = -tension;
  m[2][1] = 0.0;
  m[2][2] = tension;
  m[2][3] = 0.0;

  m[3][0] = 0.0;
  m[3][1] = 1.0;
  m[3][2] = 0.0;
  m[3][3] = 0.0;
  
  return m;
}

/**
 * Returns a Tau spline basis matrix.
 * @param bias The bias parameter for the Tau spline.
 * @param tension The tension parameter for the Tau spline.
 * @return A Tau spline basis matrix.
 */
gmMatrix4 gmMatrix4::tauBasis(double bias, double tension)
{
  gmMatrix4 m;
  
  m[0][0] = tension * (bias - 1.0);
  m[0][1] = 2.0 - tension * bias;
  m[0][2] = tension * (1.0 - bias) - 2.0;
  m[0][3] = tension * bias;

  m[1][0] = tension * 2.0 * (1.0 - bias);
  m[1][1] = tension * (3.0 * bias - 1.0) - 3.0;
  m[1][2] = 3.0 - tension;
  m[1][3] = -tension * bias;

  m[2][0] = tension * (bias - 1.0);
  m[2][1] = tension * (1.0 - 2.0 * bias);
  m[2][2] = tension * bias;
  m[2][3] = 0.0;

  m[3][0] = 0.0;
  m[3][1] = 1.0;
  m[3][2] = 0.0;
  m[3][3] = 0.0;
  
  return m;
}

/**
 * Returns a beta-spline basis matrix.
 * @param bias The bias parameter for the beta-spline.
 * @param tension The tension parameter for the beta-spline.
 * @return A beta-spline basis matrix.
 */
gmMatrix4 gmMatrix4::betaSplineBasis(double bias, double tension)
{
  gmMatrix4 m;
  double bias2, bias3, d;

  bias2 = bias * bias;
  bias3 = bias * bias2;
  d = 1.0 / (tension + 2.0 * bias3 + 4.0 * bias2 + 4.0 * bias + 2.0);

  m[0][0] = -d * 2.0 * bias3;
  m[0][1] = d * 2.0 * (tension + bias3 + bias2 + bias);
  m[0][2] = -d * 2.0 * (tension + bias2 + bias + 1.0);
  m[0][3] = d * 2.0;

  m[1][0] = d * 6.0 * bias3;
  m[1][1] = -d * 3.0 * (tension + 2.0 * bias3 + 2.0 * bias2);
  m[1][2] = d * 3.0 * (tension + 2 * bias2);
  m[1][3] = 0.0;
  
  m[2][0] = -d * 6.0 * bias3;
  m[2][1] = d * 6.0 * (bias3 - bias);
  m[2][2] = d * 6.0 * bias;
  m[2][3] = 0.0;

  m[3][0] = d * 2.0 * bias3;
  m[3][1] = d * (tension + 4 * (bias2 + bias));
  m[3][2] = d * 2.0;
  m[3][3] = 0.0;

  return m;
}

