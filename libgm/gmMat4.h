/**
 * This file contains the class declaration of the 4x4 matrix class.
 * @file  gmMat4.h
 * @date  15 June 1994
 * @author Ferdi Scheepers
 * @author Stephen F May
 */

#ifndef GMMATRIX4_H
#define GMMATRIX4_H

#include "gmUtils.h"

class gmVector3;
class gmVector4;

/**
 * The gmMatrix4 class provides methods and operators for creating and
 * working with 4x4 matrices.  The matrices are defined in row-major order,
 * which makes it easy to construct a new matrix when passing in values for
 * the elements.
 * @nosubgrouping
 */
class gmMatrix4 {

protected:
  double m_[4][4];   ///< Private storage for the elements.

public:
  /** @name Constructors
   * @{
   */
  /// Default constructor.
  gmMatrix4();
  /// Copy constructor.
  gmMatrix4(const gmMatrix4&);
  /// Constructs a new matrix object using nine doubles.
  gmMatrix4(double, double, double, double,
  	    double, double, double, double,
 	    double, double, double, double,
	    double, double, double, double);
  /** @} end of constructors */

  /** @name Array Access
   * @{
   */
  /// Returns a row of the matrix.
  double* operator [](int);
  /// Returns a const row of the matrix.
  const double* operator [](int) const;
  /** @} end of array access */

  /** @name Assignment
   * @{
   */
  /// Sets the matrix's elements to the given values.
  gmMatrix4& assign(double, double, double, double,
		    double, double, double, double,
		    double, double, double, double,
		    double, double, double, double);
  /// Sets the value of the matrix to another matrix.
  gmMatrix4& operator =(const gmMatrix4&);
  /** @} end of assignment */

  /** @name Math Operators
   * @{
   */
  /// Addition/assignment operator.
  gmMatrix4& operator +=(const gmMatrix4&);
  /// Subtraction/assignment operator.
  gmMatrix4& operator -=(const gmMatrix4&);
  /// Multiplication/assignment operator.
  gmMatrix4& operator *=(const gmMatrix4&);
  /// Multiplication/assignment with a scalar.
  gmMatrix4& operator *=(double);
  /// Division/assignment with a scalar.
  gmMatrix4& operator /=(double);

  /// Addition operator.
  gmMatrix4 operator +(const gmMatrix4&) const;
  /// Subtraction operator.
  gmMatrix4 operator -(const gmMatrix4&) const;
  /// Unary negation operator.
  gmMatrix4 operator -() const;
  /// Multiplication operator.
  gmMatrix4 operator *(const gmMatrix4&) const;
  /// Scalar multiplication operator.
  gmMatrix4 operator *(double) const;
  /// Scalar multiplication operator.
  friend gmMatrix4 operator *(double, const gmMatrix4&);
  /// Scalar division operator.
  gmMatrix4 operator /(double) const;
  /// Vector multiplication operator.
  gmVector4 operator *(const gmVector4&) const;
  /// Vector multiplication operator.
  friend gmVector4 operator *(const gmVector4&, const gmMatrix4&);

  /// Equivalence operator.
  bool operator ==(const gmMatrix4&) const;
  /// Non-equivalence operator.
  bool operator !=(const gmMatrix4&) const;
  /** @} end of math operations */

  /** @name Operations
   * @{
   */
  /// Returns the inverse of the matrix.
  gmMatrix4 inverse() const;
  /// Returns the transpose of the matrix.
  gmMatrix4 transpose() const;
  /// Returns the adjoint of the matrix.
  gmMatrix4 adjoint() const;
  
  /// Returns the determinant of the matrix.
  double determinant() const;
  /// Returns if the matrix is singular or not.
  bool isSingular() const;

  /// Transforms (multiplies) a gmVector3 by the current matrix.
  gmVector3 transform(const gmVector3&) const;
  
  /// Puts the elements of the matrix in a 2-D array of floats.
  void copyTo(float [4][4]) const;
  /// Puts the elements of the matrix in a 2-D array of doubles.
  void copyTo(double [4][4]) const;
  /// Puts the elements of the matrix in a 1-D array of floats
  void copyTo(float[16]) const;
  /// Puts the elements of the matrix in a 1-D array of doubles
  void copyTo(double[16]) const;
  /** @} end of operations */

  /** @name Transformation Matrices
   * Note that transformation matrices are formed such that any point should
   * be represented as a row, and thus premultiplied with the matrix 
   * (ie. v' = v * M).
   * @{
   */
  /// Returns the identity matrix.
  static gmMatrix4 identity();
  /// Returns a matrix corresponding to degrees rotation about an axis.
  static gmMatrix4 rotate(double, const gmVector3& axis);
  /// Returns a matrix corresponding to scaling by (x,y,z) factors.
  static gmMatrix4 scale(double, double, double);
  /// Returns a matrix corresponding to translation by (x,y,z).
  static gmMatrix4 translate(double, double, double);
  /** @} end of transformation matrices */

  /** @name Cubic Spline Basis Matrices
   * @{
   */
  /// Returns a Bezier spline basis matrix.
  static gmMatrix4 bezierBasis();
  /// Returns a B-spline basis matrix.
  static gmMatrix4 bsplineBasis();
  /// Returns a Catmull-Rom spline basis matrix.
  static gmMatrix4 catmullromBasis();
  /// Returns a Hermite spline basis matrix.
  static gmMatrix4 hermiteBasis();
  
  /// Returns a Tensed b-spline basis matrix.
  static gmMatrix4 tensedBSplineBasis(double);
  /// Returns a Cardinal spline basis matrix.
  static gmMatrix4 cardinalBasis(double);
  /// Returns a Tau spline basis matrix.
  static gmMatrix4 tauBasis(double, double);
  /// Returns a beta-spline basis matrix.
  static gmMatrix4 betaSplineBasis(double, double);
  /** @} end of cubic basis matrices */

  /** @name I/O
  * @{
   */
  /// Outputs the matrix to a stream.
  friend std::ostream & operator << ( std::ostream &, const gmMatrix4 & );
  /** @} end of I/O operations */
};

// ARRAY ACCESS
/**
 * Returns a row of the matrix.  The method uses an 'assert' call to make
 * sure that the value of the array index is valid.
 * @param i The row number of the matrix to access (zero indexed).
 * @return A pointer to an array containing the row.
 */
inline double* gmMatrix4::operator [](int i)
{
  assert(i == 0 || i == 1 || i == 2 || i == 3);
  return &m_[i][0];
}

/**
 * Returns a const row of the matrix.  The method uses an 'assert' call to
 * make sure that the value of the array index is valid.
 * @param i The row number of the matrix to access (zero indexed).
 * @return A const pointer to an array containing the row.
 */
inline const double* gmMatrix4::operator [](int i) const
{
  assert(i == 0 || i == 1 || i == 2 || i == 3);
  return &m_[i][0];
}

/**
 * Puts the elements of the matrix in a 2-D array of floats.  This is useful
 * if you need to pass the matrix as a standard 2-D C array.
 * @param f A 2-D array of floats to copy the matrix to.
 */
inline void gmMatrix4::copyTo(float f[4][4]) const
{
  f[0][0] = (float)m_[0][0]; f[0][1] = (float)m_[0][1]; 
  f[0][2] = (float)m_[0][2]; f[0][3] = (float)m_[0][3];
  f[1][0] = (float)m_[1][0]; f[1][1] = (float)m_[1][1]; 
  f[1][2] = (float)m_[1][2]; f[1][3] = (float)m_[1][3];
  f[2][0] = (float)m_[2][0]; f[2][1] = (float)m_[2][1]; 
  f[2][2] = (float)m_[2][2]; f[2][3] = (float)m_[2][3];
  f[3][0] = (float)m_[3][0]; f[3][1] = (float)m_[3][1]; 
  f[3][2] = (float)m_[3][2]; f[3][3] = (float)m_[3][3];
}

/**
 * Puts the elements of the matrix in a 2-D array of doubles.  This is useful
 * if you need to pass the matrix as a standard 2-D C array.
 * @param f A 2-D array of doubles to copy the matrix to.
 */
inline void gmMatrix4::copyTo(double f[4][4]) const
{
  f[0][0] = m_[0][0]; f[0][1] = m_[0][1]; 
  f[0][2] = m_[0][2]; f[0][3] = m_[0][3];
  f[1][0] = m_[1][0]; f[1][1] = m_[1][1]; 
  f[1][2] = m_[1][2]; f[1][3] = m_[1][3];
  f[2][0] = m_[2][0]; f[2][1] = m_[2][1]; 
  f[2][2] = m_[2][2]; f[2][3] = m_[2][3];
  f[3][0] = m_[3][0]; f[3][1] = m_[3][1]; 
  f[3][2] = m_[3][2]; f[3][3] = m_[3][3];
}

/**
 * Puts the elements of the matrix in a 1-D array of floats.  This is useful
 * if you need to pass the matrix as a standard GL array.
 * @param f A 1-D array of floats to copy the matrix to.
 */
inline void gmMatrix4::copyTo(float f[16]) const
{
  f[0] = (float)m_[0][0]; f[8] = (float)m_[0][2];
  f[1] = (float)m_[1][0]; f[9] = (float)m_[1][2];
  f[2] = (float)m_[2][0]; f[10] = (float)m_[2][2];
  f[3] = (float)m_[3][0]; f[11] = (float)m_[3][2];
  f[4] = (float)m_[0][1]; f[12] = (float)m_[0][3];
  f[5] = (float)m_[1][1]; f[13] = (float)m_[1][3];
  f[6] = (float)m_[2][1]; f[14] = (float)m_[2][3];
  f[7] = (float)m_[3][1]; f[15] = (float)m_[3][3];
}

/**
 * Puts the elements of the matrix in a 1-D array of doubles.  This is useful
 * if you need to pass the matrix as a standard GL array.
 * @param f A 1-D array of doubles to copy the matrix to.
 */
inline void gmMatrix4::copyTo(double f[16]) const
{
  f[0] = m_[0][0]; f[8] = m_[0][2];
  f[1] = m_[1][0]; f[9] = m_[1][2];
  f[2] = m_[2][0]; f[10] = m_[2][2];
  f[3] = m_[3][0]; f[11] = m_[3][2];
  f[4] = m_[0][1]; f[12] = m_[0][3];
  f[5] = m_[1][1]; f[13] = m_[1][3];
  f[6] = m_[2][1]; f[14] = m_[2][3];
  f[7] = m_[3][1]; f[15] = m_[3][3];
}

inline std::ostream & operator << ( std::ostream& os, const gmMatrix4& m )
{
  os << "{ < " << m[0][0] << " " << m[0][1] << " "
               << m[0][2] << " " << m[0][3] << " >" << std::endl;
  os << "  < " << m[1][0] << " " << m[1][1] << " "
               << m[1][2] << " " << m[1][3] << " >" << std::endl;
  os << "  < " << m[2][0] << " " << m[2][1] << " "
               << m[2][2] << " " << m[2][3] << " >" << std::endl;
  os << "  < " << m[3][0] << " " << m[3][1] << " "
               << m[3][2] << " " << m[3][3] << " > }";
  return os;
}

#endif // GMMATRIX4_H

