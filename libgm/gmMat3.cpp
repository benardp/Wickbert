/**
* This file contains method definitions for the 3x3 matrix class.
* @file  gmMat3.cpp
* @date  15 June 1994
* @author Ferdi Scheepers
* @author Stephen F May
*/

#include "gmMat3.h"

#include "gmVec2.h"
#include "gmVec3.h"

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
inline double RCD(const gmMatrix3& A, const gmMatrix3& B, int i, int j)
{
	return A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j];
}

/**
* This (private) utility function calculate the minor product difference of
* a matrix.  Given two rows and two columns, we form an "X" across these
* values, multiply the "ends" of the "X", and take their difference.  It is
* used when calculating the adjoint and determinant of a matrix.   
* @param M The matrix to use.
* @param r0 The first row.
* @param r1 The second row.
* @param c0 The first column.
* @param c1 The second column.
* @return The MINOR of the matrix
* @remarks Each of the values for r0, r1, c0, and c1 must lie within [0,2]
* since no error checking is performed.
*/
inline double MINOR(const gmMatrix3& M, int r0, int r1, int c0, int c1)
{
	return M[r0][c0] * M[r1][c1] - M[r1][c0] * M[r0][c1];
}

// CONSTRUCTORS

/** 
* Default constructor. Sets all elements of the matrix to be 0.0.
*/
gmMatrix3::gmMatrix3()
{
	assign(0,0,0, 0,0,0, 0,0,0);
}

/** 
* Copy constructor.  Creates a new matrix by using the values of another
* matrix.
* @param M The matrix to copy.
*/
gmMatrix3::gmMatrix3(const gmMatrix3& M)
{
	assign(M[0][0], M[0][1], M[0][2],
		M[1][0], M[1][1], M[1][2],
		M[2][0], M[2][1], M[2][2]);
}

/**
* Constructs a new matrix object using nine doubles.  
* @param a00 The value for element in first row, first column.
* @param a01 The value for element in first row, second column.
* @param a02 The value for element in first row, third column.
* @param a10 The value for element in second row, first column.
* @param a11 The value for element in second row, second column.
* @param a12 The value for element in second row, third column.
* @param a20 The value for element in third row, first column.
* @param a21 The value for element in third row, second column.
* @param a22 The value for element in third row, third column.
*/
gmMatrix3::gmMatrix3(double a00, double a01, double a02,
					 double a10, double a11, double a12,
					 double a20, double a21, double a22)
{
	assign(a00, a01, a02,
		a10, a11, a12,
		a20, a21, a22);
}

// ASSIGNMENT

/** 
* Sets the matrix's elements to the given values.
* @param a00 The value for element in first row, first column.
* @param a01 The value for element in first row, second column.
* @param a02 The value for element in first row, third column.
* @param a10 The value for element in second row, first column.
* @param a11 The value for element in second row, second column.
* @param a12 The value for element in second row, third column.
* @param a20 The value for element in third row, first column.
* @param a21 The value for element in third row, second column.
* @param a22 The value for element in third row, third column.
* @return A matrix with new values.
*/
gmMatrix3& gmMatrix3::assign(double a00, double a01, double a02,
							 double a10, double a11, double a12,
							 double a20, double a21, double a22)
{
	m_[0][0] = a00; m_[0][1] = a01; m_[0][2] = a02;
	m_[1][0] = a10; m_[1][1] = a11; m_[1][2] = a12;
	m_[2][0] = a20; m_[2][1] = a21; m_[2][2] = a22;
	return *this;
}

/**
* Sets the value of the matrix to another matrix.
* @param M The other matrix to use for setting values.
* @return A matrix with new values.
*/
gmMatrix3& gmMatrix3::operator =(const gmMatrix3& M)
{
	assign(M[0][0], M[0][1], M[0][2],
		M[1][0], M[1][1], M[1][2],
		M[2][0], M[2][1], M[2][2]);
	return *this;
}

// MATH OPERATORS

/**
* Addition/assignment operator.
* @param M The matrix to add to the current matrix.
* @return The matrix += M.
*/
gmMatrix3& gmMatrix3::operator +=(const gmMatrix3& M)
{
	m_[0][0] += M[0][0]; m_[0][1] += M[0][1]; m_[0][2] += M[0][2];
	m_[1][0] += M[1][0]; m_[1][1] += M[1][1]; m_[1][2] += M[1][2];
	m_[2][0] += M[2][0]; m_[2][1] += M[2][1]; m_[2][2] += M[2][2];
	return *this;
}

/**
* Subtraction/assignment operator.
* @param M The matrix to subtract from the current matrix.
* @return The matrix -= M.
*/
gmMatrix3& gmMatrix3::operator -=(const gmMatrix3& M)
{
	m_[0][0] -= M[0][0]; m_[0][1] -= M[0][1]; m_[0][2] -= M[0][2];
	m_[1][0] -= M[1][0]; m_[1][1] -= M[1][1]; m_[1][2] -= M[1][2];
	m_[2][0] -= M[2][0]; m_[2][1] -= M[2][1]; m_[2][2] -= M[2][2];
	return *this;
}

/**
* Multiplication/assignment operator.
* @param M The matrix to multiply with the current matrix.
* @return The matrix *= M.
*/
gmMatrix3& gmMatrix3::operator *=(const gmMatrix3& M)
{
	assign(RCD(*this, M, 0, 0), RCD(*this, M, 0, 1), RCD(*this, M, 0, 2),
		RCD(*this, M, 1, 0), RCD(*this, M, 1, 1), RCD(*this, M, 1, 2),
		RCD(*this, M, 2, 0), RCD(*this, M, 2, 1), RCD(*this, M, 2, 2));
	return *this;
}

/**
* Multiplication/assignment operator with a scalar.
* @param d The scalar value multiplied with each element of the matrix.
* @return d*M.
*/
gmMatrix3& gmMatrix3::operator *=(double d)
{
	m_[0][0] *= d; m_[0][1] *= d; m_[0][2] *= d;
	m_[1][0] *= d; m_[1][1] *= d; m_[1][2] *= d;
	m_[2][0] *= d; m_[2][1] *= d; m_[2][2] *= d;
	return *this;
}

/**
* Division/assignment operator with a scalar.
* @param d The scalar value divided into each element of the matrix.
* @return 1/d * M.
*/
gmMatrix3& gmMatrix3::operator /=(double d)
{
	double di = 1 / d;
	m_[0][0] *= di; m_[0][1] *= di; m_[0][2] *= di;
	m_[1][0] *= di; m_[1][1] *= di; m_[1][2] *= di;
	m_[2][0] *= di; m_[2][1] *= di; m_[2][2] *= di;
	return *this;
}

/**
* Addition operator.  Adds two matrices together.
* @param M The matrix to add to the current matrix.
* @return A new matrix resulting from the addition of two matrices.
*/
gmMatrix3 gmMatrix3::operator +(const gmMatrix3& M) const
{
	return gmMatrix3(m_[0][0] + M[0][0], m_[0][1] + M[0][1], m_[0][2] + M[0][2],
		m_[1][0] + M[1][0], m_[1][1] + M[1][1], m_[1][2] + M[1][2],
		m_[2][0] + M[2][0], m_[2][1] + M[2][1], m_[2][2] + M[2][2]);
}

/**
* Subtraction operator.  Subtracts one matrix from another.
* @param M The matrix to subtract from the current matrix.
* @return A new matrix resulting from the difference of two matrices.
*/
gmMatrix3 gmMatrix3::operator -(const gmMatrix3& M) const
{
	return gmMatrix3(m_[0][0] - M[0][0], m_[0][1] - M[0][1], m_[0][2] - M[0][2],
		m_[1][0] - M[1][0], m_[1][1] - M[1][1], m_[1][2] - M[1][2],
		m_[2][0] - M[2][0], m_[2][1] - M[2][1], m_[2][2] - M[2][2]);
}

/**
* Unary negation operator. 
* @return A new matrix with all elements negated.
*/
gmMatrix3 gmMatrix3::operator -() const
{
	return gmMatrix3(-m_[0][0], -m_[0][1], -m_[0][2],
		-m_[1][0], -m_[1][1], -m_[1][2],
		-m_[2][0], -m_[2][1], -m_[2][2]);
}

/**
* Multiplication operator.	Multiplies two matrices together.
* @param M The matrix to multiply with the current matrix.
* @return A new matrix resulting from the product of two matrices.
*/
gmMatrix3 gmMatrix3::operator *(const gmMatrix3& M) const
{
	return
		gmMatrix3(RCD(*this, M, 0, 0), RCD(*this, M, 0, 1), RCD(*this, M, 0, 2),
		RCD(*this, M, 1, 0), RCD(*this, M, 1, 1), RCD(*this, M, 1, 2),
		RCD(*this, M, 2, 0), RCD(*this, M, 2, 1), RCD(*this, M, 2, 2));
}

/**
* Scalar multiplication operator.  Multiplies a scalar with the matrix.
* @param d The scalar value multiplied with each element of the matrix.
* @return A new matrix with d*M.
*/
gmMatrix3 gmMatrix3::operator *(double d) const
{
	return gmMatrix3(m_[0][0] * d, m_[0][1] * d, m_[0][2] * d,
		m_[1][0] * d, m_[1][1] * d, m_[1][2] * d,
		m_[2][0] * d, m_[2][1] * d, m_[2][2] * d);
}

/**
* Scalar division operator.  Divides a scalar into the matrix.
* @param d The scalar value divided into each element of the matrix.
* @return A new matrix with 1/d * M.
*/
gmMatrix3 gmMatrix3::operator /(double d) const
{
	assert(!gmIsZero(d));
	double di = 1 / d;
	return gmMatrix3(m_[0][0] * di, m_[0][1] * di, m_[0][2] * di,
		m_[1][0] * di, m_[1][1] * di, m_[1][2] * di,
		m_[2][0] * di, m_[2][1] * di, m_[2][2] * di);
}

/**
* Scalar multiplication operator.  Multiplies a scalar with the matrix.
* @param d The scalar value multiplied with each element of the matrix.
* @param M The matrix to multiply.
* @return A new matrix with d*M.
*/
gmMatrix3 operator *(double d, const gmMatrix3& M)
{
	return gmMatrix3(M[0][0] * d, M[0][1] * d, M[0][2] * d,
		M[1][0] * d, M[1][1] * d, M[1][2] * d,
		M[2][0] * d, M[2][1] * d, M[2][2] * d);
}

/**
* Vector multiplication operator.  Multiplies a vector with the matrix.
* @param v The vector to multiply.
* @return A new vector resulting from M*v.
*/
gmVector3 gmMatrix3::operator *(const gmVector3& v) const
{
	return gmVector3(m_[0][0] * v[0] + m_[0][1] * v[1] + m_[0][2] * v[2],
		m_[1][0] * v[0] + m_[1][1] * v[1] + m_[1][2] * v[2],
		m_[2][0] * v[0] + m_[2][1] * v[1] + m_[2][2] * v[2]);
}

/**
* Vector multiplication operator.  Multiplies a vector with the matrix.
* @param v The vector to multiply.
* @param M The matrix to multiply.
* @return A new vector resulting from M*v.
*/
gmVector3 operator *(const gmVector3& v, const gmMatrix3& M)
{
	return gmVector3(v[0] * M[0][0] + v[1] * M[1][0] + v[2] * M[2][0],
		v[0] * M[0][1] + v[1] * M[1][1] + v[2] * M[2][1],
		v[0] * M[0][2] * v[1] * M[1][2] + v[2] * M[2][2]);
}

/**
* Equivalence operator.  Checks to see if two matrices are equal.
* @param v The matrix to check for equivalence.
* @return True if all elements are fuzzy equal, false otherwise.
*/
bool gmMatrix3::operator ==(const gmMatrix3& M) const
{
	return(gmFuzEQ(m_[0][0], M[0][0]) &&
		gmFuzEQ(m_[0][1], M[0][1]) &&
		gmFuzEQ(m_[0][2], M[0][2]) &&
		
		gmFuzEQ(m_[1][0], M[1][0]) &&
		gmFuzEQ(m_[1][1], M[1][1]) &&
		gmFuzEQ(m_[1][2], M[1][2]) &&
		
		gmFuzEQ(m_[2][0], M[2][0]) &&
		gmFuzEQ(m_[2][1], M[2][1]) &&
		gmFuzEQ(m_[2][2], M[2][2]));
}

/**
* Non-equivalence operator.  Checks to see if two matrices are not equal.
* @param v The matrix to check for non-equivalence.
* @return False if all elements are fuzzy equal, true otherwise.
*/
bool gmMatrix3::operator !=(const gmMatrix3& M) const
{
	return (!(*this == M));
}

// OPERATIONS

/**
* Returns the inverse of the matrix.  This method makes sure that the
* matrix can be inverted by verifying that the matrix is not singular.
* @return A new matrix which is the inverse of the original.
*/
gmMatrix3 gmMatrix3::inverse() const
{
	assert(!isSingular());
	return adjoint() * gmInv(determinant());
}

/**
* Returns the transpose of the matrix.	This simply swaps elements across
* the diagonal.
* @return A new matrix with elements transposed.
*/
gmMatrix3 gmMatrix3::transpose() const
{
	return gmMatrix3(m_[0][0], m_[1][0], m_[2][0],
		m_[0][1], m_[1][1], m_[2][1],
		m_[0][2], m_[1][2], m_[2][2]);
}

/**
* Returns the adjoint of the matrix.
* @return The adjoint of the matrix.
*/
gmMatrix3 gmMatrix3::adjoint() const
{
	gmMatrix3 A;
	
	A[0][0] =  MINOR(*this, 1, 2, 1, 2);
	A[0][1] = -MINOR(*this, 0, 2, 1, 2);
	A[0][2] =  MINOR(*this, 0, 1, 1, 2);
	A[1][0] = -MINOR(*this, 1, 2, 0, 2);
	A[1][1] =  MINOR(*this, 0, 2, 0, 2);
	A[1][2] = -MINOR(*this, 0, 1, 0, 2);
	A[2][0] =  MINOR(*this, 1, 2, 0, 1);
	A[2][1] = -MINOR(*this, 0, 2, 0, 1);
	A[2][2] =  MINOR(*this, 0, 1, 0, 1);
	
	return A;
}

/**
* Returns the determinant of the matrix.
* @return The determinant of the matrix.
*/
double gmMatrix3::determinant() const
{
	return m_[0][0] * MINOR(*this, 1, 2, 1, 2) -
		m_[0][1] * MINOR(*this, 1, 2, 0, 2) +
		m_[0][2] * MINOR(*this, 1, 2, 0, 1);
}

/**
* Returns if the matrix is singular or not.  This method checks to see if
* the determinant of the matrix equals 0.0.
* @return True if the matrix is singular (determinant == 0.0), false
* otherwise.
*/
bool gmMatrix3::isSingular() const
{
	return gmIsZero(determinant());
}

/**
* Transforms (multiplies) a gmVector2 by the current matrix.
* @param v The gmVector2 to multiply.
* @return A new vector reflecting v * the first two columns of M.
*/
gmVector2 gmMatrix3::transform(const gmVector2& v) const
{
	return gmVector2(v[0] * m_[0][0] + v[1] * m_[1][0] + m_[2][0],
		v[0] * m_[0][1] + v[1] * m_[1][1] + m_[2][1]);
}

// TRANSFORMATION MATRICES

/**
* Returns the identity matrix.
* @return A new matrix with ones along the diagonal and zeros elsewhere.
*/
gmMatrix3 gmMatrix3::identity()
{
	return gmMatrix3(1, 0, 0,
		0, 1, 0,
		0, 0, 1);
}

/**
* Returns a matrix corresponding to degrees rotation in the xy-plane.
* @param angle The angle in degrees to rotate in the xy-plane.
* @return A new rotation matrix.
*/
gmMatrix3 gmMatrix3::rotate(double angle)
{
	double sine = sin(gmRadians(angle));
	double cosine = cos(gmRadians(angle));
	
	return gmMatrix3( cosine, sine,   0,
		-sine,	 cosine, 0,
		0,		0,		1);
}

/**
* Returns a matrix corresponding to scaling by (x,y) factors.
* @param x Scaling factor for x direction.
* @param y Scaling factor for y direction.
* @return A new scaling matrix.
*/
gmMatrix3 gmMatrix3::scale(double x, double y)
{
	return gmMatrix3(x, 0, 0,
		0, y, 0,
		0, 0, 1);
}

/**
* Returns a matrix corresponding to translation by (x,y).
* @param x Amount to translate in x direction.
* @param y Amount to translate in y direction.
* @return A new translation matrix.
* @remark Note that this assumes you are working in 2-D since this returns
* a matrix using homogeneous coordinate space.
*/
gmMatrix3 gmMatrix3::translate(double x, double y)
{
	return gmMatrix3(1, 0, 0,
		0, 1, 0,
		x, y, 1);
}

