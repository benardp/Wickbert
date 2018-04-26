/**
* This file contains the class declaration of the 3x3 matrix class.
* @file  gmMat3.h
* @date  15 June 1994
* @author Ferdi Scheepers
* @author Stephen F May
*/

#ifndef GMMATRIX3_H
#define GMMATRIX3_H

#include "gmUtils.h"

class gmVector2;
class gmVector3;

/**
* The gmMatrix3 class provides methods and operators for creating and
* working with 3x3 matrices.  The matrices are defined in row-major order,
* which makes it easy to construct a new matrix when passing in values for
* the elements.
* @nosubgrouping
*/
class gmMatrix3 {
  
protected:
  double m_[3][3];    ///< Private storage for the elements.
  
public:
/** @name Constructors
* @{
  */
  /// Default constructor.
  gmMatrix3();
  /// Copy constructor.
  gmMatrix3(const gmMatrix3&);
  /// Constructs a new matrix object using nine doubles.
  gmMatrix3(double, double, double,
    double, double, double,
    double, double, double);
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
  gmMatrix3& assign(double, double, double,
    double, double, double,
    double, double, double);
  /// Sets the value of the matrix to another matrix.
  gmMatrix3& operator =(const gmMatrix3&);
  /** @} end of assignment */
  
  /** @name Math Operators
  * @{
  */
  /// Addition/assignment operator.
  gmMatrix3& operator +=(const gmMatrix3&);
  /// Subtraction/assignment operator.
  gmMatrix3& operator -=(const gmMatrix3&);
  /// Multiplication/assignment operator.
  gmMatrix3& operator *=(const gmMatrix3&);
  /// Multiplication/assignment with a scalar.
  gmMatrix3& operator *=(double);
  /// Division/assignment with a scalar.
  gmMatrix3& operator /=(double);
  
  /// Addition operator.
  gmMatrix3 operator +(const gmMatrix3&) const;
  /// Subtraction operator.
  gmMatrix3 operator -(const gmMatrix3&) const;
  /// Unary negation operator.
  gmMatrix3 operator -() const;
  /// Multiplication operator.
  gmMatrix3 operator *(const gmMatrix3&) const;
  /// Scalar multiplication operator.
  gmMatrix3 operator *(double) const;
  /// Scalar multiplication operator.
  friend gmMatrix3 operator *(double, const gmMatrix3&);
  /// Scalar division operator.
  gmMatrix3 operator /(double) const;
  /// Vector multiplication operator.
  gmVector3 operator *(const gmVector3&) const;
  /// Vector multiplication operator.
  friend gmVector3 operator *(const gmVector3&, const gmMatrix3&);
  
  /// Equivalence operator.
  bool operator ==(const gmMatrix3&) const;
  /// Non-equivalence operator.
  bool operator !=(const gmMatrix3&) const;
  /** @} end of math operations */
  
  /** @name Operations
  * @{
  */
  /// Returns the inverse of the matrix.
  gmMatrix3 inverse() const;
  /// Returns the transpose of the matrix.
  gmMatrix3 transpose() const;
  /// Returns the adjoint of the matrix.
  gmMatrix3 adjoint() const;
  
  /// Returns the determinant of the matrix.
  double determinant() const;
  /// Returns if the matrix is singular or not.
  bool isSingular() const;
  
  /// Transforms (multiplies) a gmVector2 by the current matrix.
  gmVector2 transform(const gmVector2&) const;
  
  /// Puts the elements of the matrix in a 2-D array of floats.
  void copyTo(float [3][3]) const;
  /// Puts the elements of the matrix in a 2-D array of doubles.
  void copyTo(double [3][3]) const;
  /** @} end of operations */
  
  /** @name Transformation Matrices
        * Note that transformation matrices are formed such that any point
        * should be represented as a row, and thus premultiplied with the
        * matrix * (ie. v' = v * M).
  * @{
  */
  /// Returns the identity matrix.
  static gmMatrix3 identity();
  /// Returns a matrix corresponding to degrees rotation in xy-plane.
  static gmMatrix3 rotate(double);
  /// Returns a matrix corresponding to scaling by (x,y) factors.
  static gmMatrix3 scale(double, double);
  /// Returns a matrix corresponding to translation by (x,y).
  static gmMatrix3 translate(double, double);
  /** @} end of transformation matrices */

        /** @name I/O
  * @{
         */
        /// Outputs the matrix to a stream.
        friend std::ostream & operator << ( std::ostream &, const gmMatrix3 & );
        /** @} end of I/O operations */
};

// ARRAY ACCESS

/**
* Returns a row of the matrix.  The method uses an 'assert' call to make
* sure that the value of the array index is valid.
* @param i The row number of the matrix to access (zero indexed).
* @return A pointer to an array containing the row.
*/
inline double* gmMatrix3::operator [](int i)
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
inline const double* gmMatrix3::operator [](int i) const
{
  assert(i == 0 || i == 1 || i == 2 || i == 3);
  return &m_[i][0];
}

/**
* Puts the elements of the matrix in a 2-D array of floats.  This is useful
* if you need to pass the matrix as a standard 2-D C array.
* @param f A 2-D array of floats to copy the matrix to.
*/  
inline void gmMatrix3::copyTo(float f[3][3]) const
{
  f[0][0] = (float)m_[0][0]; f[0][1] = (float)m_[0][1]; f[0][2] = (float)m_[0][2];
  f[1][0] = (float)m_[1][0]; f[1][1] = (float)m_[1][1]; f[1][2] = (float)m_[1][2];
  f[2][0] = (float)m_[2][0]; f[2][1] = (float)m_[2][1]; f[2][2] = (float)m_[2][2];
}

/**
* Puts the elements of the matrix in a 2-D array of doubles.  This is useful
* if you need to pass the matrix as a standard 2-D C array.
* @param f A 2-D array of floats to copy the matrix to.
*/  
inline void gmMatrix3::copyTo(double f[3][3]) const
{
  f[0][0] = m_[0][0]; f[0][1] = m_[0][1]; f[0][2] = m_[0][2];
  f[1][0] = m_[1][0]; f[1][1] = m_[1][1]; f[1][2] = m_[1][2];
  f[2][0] = m_[2][0]; f[2][1] = m_[2][1]; f[2][2] = m_[2][2];
}

inline std::ostream & operator << ( std::ostream& os, const gmMatrix3& m )
{
  os << "{ < " << m[0][0] <<" "<< m[0][1] <<" "<< m[0][2] << " >" << std::endl;
  os << "  < " << m[1][0] <<" "<< m[1][1] <<" "<< m[1][2] << " >" << std::endl;
  os << "  < " << m[2][0] <<" "<< m[2][1] <<" "<< m[2][2] << " > }";
  return os;
}

#endif

