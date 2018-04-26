/**
 * This file contains functions for 3-element vectors.
 * @file  gmVec3.h
 * @date  15 June 1994
 * @author Ferdi Scheepers
 * @author Stephen F May
 */

#ifndef GMVECTOR3_H
#define GMVECTOR3_H

#include "gmUtils.h"

/**
 * The gmVector3 class provides methods and operators for creating and
 * working with vectors containing three elements.  You can access the
 * elements of the vector as if the vector was an array.
 * @nosubgrouping
 */
class gmVector3 {

protected:
  double v_[3];  ///< Private storage for the elements

public:
  /** @name Constructors
   * @{
   */
  /// Default constructor.
  gmVector3();
  /// Copy constructor.
  gmVector3(const gmVector3&);
  /// Constructs a new vector object using three doubles.
  gmVector3(double, double, double);
  /** @} end of constructors */

  /** @name Array Access
   * @{
   */
  /// Returns an element of the vector.
  double& operator [](int);
  /// Returns a const element of the vector.
  const double& operator [](int) const;
  /** @} end of array access */

  /** @name Assignment
   * @{
   */
  /// Sets the vector's elements to the given values.
  gmVector3& assign(double, double, double);
  /// Sets the value of the vector to another vector.
  gmVector3& operator =(const gmVector3&);
  /** @} end of assignment */

  /** @name Math Operators
   * @{
   */
  /// Addition/assignment operator.
  gmVector3& operator +=(const gmVector3&);
  /// Subtraction/assignment operator.
  gmVector3& operator -=(const gmVector3&);
  /// Multiplication/assignment operator.
  gmVector3& operator *=(double);
  /// Division/assignment operator.
  gmVector3& operator /=(double);

  /// Addition operator.
  gmVector3 operator +(const gmVector3&) const;
  /// Subtraction operator.
  gmVector3 operator -(const gmVector3&) const;
  /// Unary negation operator.
  gmVector3 operator -() const;
  /// Scalar multiplication operator.
  gmVector3 operator *(double) const;
  /// Scalar division operator.
  gmVector3 operator /(double) const;

  /// Scalar multiplication operator.
friend gmVector3 operator *(double, const gmVector3&);

  /// Equivalence operator.
  bool operator ==(const gmVector3&) const;
  /// Non-equivalence operator.
  bool operator !=(const gmVector3&) const;
  /** @} end of math operations */

  /** @name Operations
   * @{
   */
  /// Clamps each element of the vector to a range.
  gmVector3& clamp(double, double);
  /// Returns the length of the vector.
  double length() const;
  /// Returns the length of the vector squared.
  double lengthSquared() const;
  /// Returns a normalized copy of the vector.
  gmVector3& normalize();
  /// Returns a normalized copy of the vector.
  gmVector3& normalizeNonZero();

  /// Puts the elements of the vector in an array of floats.
  void copyTo(float [3]) const;
  /// Puts the elements of the vector in an array of doubles.
  void copyTo(double [3]) const;

  /// Returns the cross product of two vectors.
friend gmVector3 cross(const gmVector3&, const gmVector3&);
  /// Returns the distance between two vectors.
friend double distance(const gmVector3&, const gmVector3&);
  /// Returns the distance between two vectors squared.
friend double distanceSquared(const gmVector3&, const gmVector3&);
  /// Returns the dot product of two vectors.
friend double dot(const gmVector3&, const gmVector3&);
  /// Returns a linear interpolation between two vectors.
friend gmVector3 lerp(double, const gmVector3&, const gmVector3&);
  /** @} end operations */

  /** @name Output
   * @{
   */
  /// Output stream concatenation operator.
friend std::ostream & operator << ( std::ostream &, const gmVector3 & );
  /** @} end output */
};

// CONSTRUCTORS

/**
 * Default constructor. Sets all elements of the vector to be 0.0.
 */
inline gmVector3::gmVector3()
{
  v_[0] = v_[1] = v_[2] = 0;
}

/**
 * Copy constructor.  Creates a new vector by using the values of another
 * vector.
 * @param v The vector to copy.
 */
inline gmVector3::gmVector3(const gmVector3& v)
{

	v_[0] = v.v_[0]; v_[1] = v.v_[1]; v_[2] = v.v_[2];
}

/**
 * Constructs a new vector object using three doubles.
 * @param x The value to set for vector[0].
 * @param y The value to set for vector[1].
 * @param z The value to set for vector[2].
 */
inline gmVector3::gmVector3(double x, double y, double z)
{
  v_[0] = x; v_[1] = y; v_[2] = z;
}

// ARRAY ACCESS

/**
 * Returns an element of the vector.  The method uses an 'assert' call to
 * make sure that the value of the array index is valid.
 * @param i The vector element to access (zero indexed).
 * @return The value of the vector element.
 */
inline double& gmVector3::operator [](int i) 
{
  assert(i == 0 || i == 1 || i == 2);
  return v_[i];
}

/**
 * Returns a const element of the vector.  The method uses an 'assert' call
 * to make sure that the value of the array index is valid.  
 * @param i The vector element to access (zero indexed).
 * @return The const value of the vector element.
 */
inline const double& gmVector3::operator [](int i) const
{
  assert(i == 0 || i == 1 || i == 2);
  return v_[i];
}

// ASSIGNMENT

/**
 * Sets the vector's elements to the given values.
 * @param x The value to set for vector[0].
 * @param y The value to set for vector[1].
 * @param z The value to set for vector[2].
 * @return The vector with new values.
 */
inline gmVector3& gmVector3::assign(double x, double y, double z)
{
  v_[0] = x; v_[1] = y; v_[2] = z;
  return *this;
}

/**
 * Sets the value of the vector to another vector.
 * @param v The other vector to use for setting values.
 * @return The vector with new values.
 */
inline gmVector3& gmVector3::operator =(const gmVector3& v)
{
  v_[0] = v[0]; v_[1] = v[1]; v_[2] = v[2];
  return *this;
}

// MATH OPERATORS

/**
 * Addition/assignment operator.
 * @param v The vector to add to the current vector.
 * @return The vector += v.
 */
inline gmVector3& gmVector3::operator +=(const gmVector3& v)
{
  v_[0] += v[0]; v_[1] += v[1]; v_[2] += v[2];
  return *this;
}

/**
 * Subtraction/assignment operator.
 * @param v The vector to subtraction from the current vector.
 * @return The vector -= v.
 */
inline gmVector3& gmVector3::operator -=(const gmVector3& v)
{
  v_[0] -= v[0]; v_[1] -= v[1]; v_[2] -= v[2];
  return *this;
}

/**
 * Multiplication/assignment operator with a scalar.
 * @param c The scalar to multiply with the current vector.
 * @return The vector *= c.
 */
inline gmVector3& gmVector3::operator *=(double c)
{
  v_[0] *= c; v_[1] *= c; v_[2] *= c;
  return *this;
}

/**
 * Division/assignment operator with a scalar.
 * @param c The scalar to divide by the current vector.
 * @return The vector /= c.
 * @note Modified code to be more efficient by inverting c once and
 * multiplying each term.
 */
inline gmVector3& gmVector3::operator /=(double c)
{
  assert(!gmIsZero(c));
  double ci = 1.0/c;
  v_[0] *= ci; v_[1] *= ci; v_[2] *= ci;
  return *this;
}

/**
 * Addition operator.  Adds two vectors together.
 * @param v The vector to add to the current vector.
 * @return A new vector resulting from the addition of two vectors.
 */
inline gmVector3 gmVector3::operator +(const gmVector3& v) const
{
  return gmVector3(v_[0] + v[0], v_[1] + v[1], v_[2] + v[2]);
}

/**
 * Subtraction operator.  Subtracts one vector from another.
 * @param v The vector to subtract from the current vector.
 * @return A new vector resulting from the difference of two vectors.
 */
inline gmVector3 gmVector3::operator -(const gmVector3& v) const
{
  return gmVector3(v_[0] - v[0], v_[1] - v[1], v_[2] - v[2]);
}

/**
 * Unary negation operator. 
 * @return A new vector with all elements negated.
 */
inline gmVector3 gmVector3::operator -() const
{
  return gmVector3(-v_[0], -v_[1], -v_[2]);
}

/**
 * Scalar multiplication operator.  Multiplies a vector with a double.
 * @param c The double to multiply with the vector.
 * @return A new vector resulting from the multiplication of c with each
 * element of the vector.
 */
inline gmVector3 gmVector3::operator *(double c) const
{
  return gmVector3(v_[0] * c, v_[1] * c, v_[2] * c);
}

/**
 * Scalar division operator.  Divides a vector by a double.
 * @param c The double to divide through the vector.
 * @return A new vector resulting from the division of c with each
 * element of the vector.
 * @note Modified code to be more efficient by inverting c once and
 * multiplying each term.
 */
inline gmVector3 gmVector3::operator /(double c) const
{
  assert(!gmIsZero(c));
  double ci = 1.0/c;
  return gmVector3(v_[0] * ci, v_[1] * ci, v_[2] * ci);
}

/**
 * Scalar multiplication operator.  Multiplies a vector with a double.
 * @param c The double to multiply with the vector.
 * @param v The vector to multiply with the double.
 * @return A new vector resulting from the multiplication of c with each
 * element of the vector v.
 */
inline gmVector3 operator *(double c, const gmVector3& v)
{
  return gmVector3(c * v[0], c * v[1], c * v[2]);
}

/**
 * Equivalence operator.  Checks to see if two vectors are equal.
 * @param v The vector to check for equivalence.
 * @return True if all elements are fuzzy equal, false otherwise.
 */
inline bool gmVector3::operator ==(const gmVector3& v) const
{
  return
    (gmFuzEQ(v_[0], v[0]) && gmFuzEQ(v_[1], v[1]) && gmFuzEQ(v_[2], v[2]));
}

/**
 * Non-equivalence operator.  Checks to see if two vectors are not equal.
 * @param v The vector to check for non-equivalence.
 * @return False if all elements are fuzzy equal, true otherwise.
 */
inline bool gmVector3::operator !=(const gmVector3& v) const
{
  return (!(*this == v));
}

// OPERATIONS

/** 
 * Clamps each element of the vector to a range.
 * @param lo The low value of the range.
 * @param hi The high value of the range.
 * @return A new vector with its elements clamped to the range [lo,hi].
 * @remarks This method does not check to verify lo < hi.
 */
inline gmVector3& gmVector3::clamp(double lo, double hi)
{
  gmClamp(v_[0], lo, hi); gmClamp(v_[1], lo, hi); gmClamp(v_[2], lo, hi);
  return *this;
}

/**
 * Returns the length of the vector.
 * @return The length of the vector.
 */
inline double gmVector3::length() const
{
  return sqrt(gmSqr(v_[0]) + gmSqr(v_[1]) + gmSqr(v_[2]));
}

/**
 * Returns the length of the vector squared.
 * @return The length of the vector squared.
 */
inline double gmVector3::lengthSquared() const
{
  return gmSqr(v_[0]) + gmSqr(v_[1]) + gmSqr(v_[2]);
}

/**
 * Returns a normalized copy of the vector.  This method tests to make sure
 * the length of the vector is not zero.  If the length of the vector is
 * zero, it simply returns the vector without doing the normalization.
 * @return A new vector with length = 1.0 unless the length of the vector
 * was zero.
 */
inline gmVector3& gmVector3::normalize()
{
  double lensqr = lengthSquared();
  if (!gmIsZero(lensqr))
    *this /= sqrt(lensqr);
  return *this;
}

/**
 * Returns a normalized copy of the vector.  This method uses an 'assert' to
 * make sure that the length of the vector is not zero (instead of always
 * doing an 'if' test).
 * @return A new vector with length = 1.0.
 */
inline gmVector3& gmVector3::normalizeNonZero()
{
  double lensqr = lengthSquared();
  assert(!gmIsZero(lensqr));
  *this /= sqrt(lensqr);
  return *this;
}

/**
 * Puts the elements of the vector in an array of floats.  This is useful
 * if you need to pass the vector as a standard C array.
 * @param f An array of floats to copy the vector to.
 */
inline void gmVector3::copyTo(float f[3]) const
{
  f[0] = (float)v_[0]; f[1] = (float)v_[1]; f[2] = (float)v_[2];
}

/**
 * Puts the elements of the vector in an array of doubles.  This is useful
 * if you need to pass the vector as a standard C array.
 * @param f An array of doubles to copy the vector to.
 */
inline void gmVector3::copyTo(double f[3]) const
{
  f[0] = v_[0]; f[1] = v_[1]; f[2] = v_[2];
}

/**
 * Returns the cross product of two vectors.
 * @param v1 The first vector (point).
 * @param v2 The second vector (point).
 * @param return A new vector with v1 X v2.
 */
inline gmVector3 cross(const gmVector3& v1, const gmVector3& v2)
{
  return gmVector3(v1[1] * v2[2] - v1[2] * v2[1],
                   v1[2] * v2[0] - v1[0] * v2[2],
                   v1[0] * v2[1] - v1[1] * v2[0]);
}

/**
 * Returns the distance between two vectors.  
 * @param v1 The first vector (point).
 * @param v2 The second vector (point).
 * @return The distance between the two vectors (points).
 */
inline double distance(const gmVector3& v1, const gmVector3& v2)
{
  return
    sqrt(gmSqr(v1[0] - v2[0]) + gmSqr(v1[1] - v2[1]) + gmSqr(v1[2] - v2[2]));
}

/**
 * Returns the distance between two vectors squared.  
 * @param v1 The first vector (point).
 * @param v2 The second vector (point).
 * @return The distance squared between the two vectors (points).
 */
inline double distanceSquared(const gmVector3& v1, const gmVector3& v2)
{
  return gmSqr(v1[0] - v2[0]) + gmSqr(v1[1] - v2[1]) + gmSqr(v1[2] - v2[2]);
}

/** 
 * Returns the dot product of two vectors.
 * @param v1 The first vector (point).
 * @param v2 The second vector (point).
 * @return v1 dot v2.
 */
inline double dot(const gmVector3& v1, const gmVector3& v2)
{
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

/** 
 * Returns a linear interpolation between two vectors. 
 * @param f The number at which to evaluate the linear interpolation.
 * @param v1 The vector at which f = 0.
 * @param v2 The vector at which f = 1.
 * @return The vector representing the linear interpolation between v1 and
 * v2 evaluated at f.
 */
inline gmVector3 lerp(double f, const gmVector3& v1, const gmVector3& v2)
{
  return v1 + ((v2 - v1) * f);
}

// OUTPUT

/**
 * Output stream concatenation operator.  This allows the text output of a
 * vector object using the << operator.
 * @param os The output stream to print to.
 * @param v The vector to output.
 * @return The resulting output stream.
 */
inline std::ostream & operator << ( std::ostream& os, const gmVector3& v)
{
  os << "< " << v[0] << " " << v[1] << " " << v[2] << " >";
  return os;
}

#endif 

