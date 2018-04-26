/**
 * This file contains functions for 2-element vectors.
 * @file  gmVec2.h
 * @date  15 June 1994
 * @author Ferdi Scheepers
 * @author Stephen F May
 */

#ifndef GMVECTOR2_H
#define GMVECTOR2_H

#include "gmUtils.h"

/**
 * The gmVector2 class provides methods and operators for creating and
 * working with vectors containing two elements.  You can access the
 * elements of the vector as if the vector was an array.
 * @nosubgrouping
 */
class gmVector2 {

protected:
  double v_[2];  ///< Private storage for the elements

public:
  /** @name Constructors
   * @{
   */
  /// Default constructor.
  gmVector2();
  /// Copy constructor.
  gmVector2(const gmVector2&);
  /// Constructs a new vector object using two doubles.
  gmVector2(double, double);
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
  gmVector2& assign(double, double);
  /// Sets the value of the vector to another vector.
  gmVector2& operator =(const gmVector2&);
  /** @} end of assignment */

  /** @name Math Operators
   * @{
   */
  /// Addition/assignment operator.
  gmVector2& operator +=(const gmVector2&);
  /// Subtraction/assignment operator.
  gmVector2& operator -=(const gmVector2&);
  /// Multiplication/assignment operator.
  gmVector2& operator *=(double);
  /// Division/assignment operator.
  gmVector2& operator /=(double);

  /// Addition operator.
  gmVector2 operator +(const gmVector2&) const;
  /// Subtraction operator.
  gmVector2 operator -(const gmVector2&) const;
  /// Unary negation operator.
  gmVector2 operator -() const;
  /// Scalar multiplication operator.
  gmVector2 operator *(double) const;
  /// Scalar division operator.
  gmVector2 operator /(double) const;

  /// Scalar multiplication operator.
friend gmVector2 operator *(double, const gmVector2&);

  /// Equivalence operator.
  bool operator ==(const gmVector2&) const;
  /// Non-equivalence operator.
  bool operator !=(const gmVector2&) const;
  /** @} end of math operations */

  /** @name Operations
   * @{
   */
  /// Clamps each element of the vector to a range.
  gmVector2& clamp(double, double);
  /// Returns the length of the vector.
  double length() const;
  /// Returns the length of the vector squared.
  double lengthSquared() const;
  /// Returns a normalized copy of the vector.
  gmVector2& normalize();
  /// Returns a normalized copy of the vector.
  gmVector2& normalizeNonZero();

  /// Puts the elements of the vector in an array of floats.
  void copyTo(float [2]) const;
  /// Puts the elements of the vector in an array of doubles.
  void copyTo(double [2]) const;

  /// Returns the distance between two vectors.
friend double distance(const gmVector2&, const gmVector2&);
  /// Returns the distance between two vectors squared.
friend double distanceSquared(const gmVector2&, const gmVector2&);
  /// Returns the dot product of two vectors.
friend double dot(const gmVector2&, const gmVector2&);
  /// Returns a linear interpolation between two vectors.
friend gmVector2 lerp(double, const gmVector2&, const gmVector2&);
  /** @} end operations */

  /** @name Output
   * @{
   */
  /// Output stream concatenation operator.
friend std::ostream & operator << ( std::ostream &, const gmVector2 & );
  /** @} end output */
};

// CONSTRUCTORS

/**
 * Default constructor. Sets all elements of the vector to be 0.0.
 */
inline gmVector2::gmVector2()
{
  v_[0] = v_[1] = 0;
}

/**
 * Copy constructor.  Creates a new vector by using the values of another
 * vector.
 * @param v The vector to copy.
 */
inline gmVector2::gmVector2(const gmVector2& v)
{
  v_[0] = v.v_[0]; v_[1] = v.v_[1];
}

/**
 * Constructs a new vector object using two doubles.
 * @param x The value to set for vector[0].
 * @param y The value to set for vector[1].
 */
inline gmVector2::gmVector2(double x, double y)
{
  v_[0] = x; v_[1] = y;
}

// ARRAY ACCESS

/**
 * Returns an element of the vector.  The method uses an 'assert' call to
 * make sure that the value of the array index is valid.
 * @param i The vector element to access (zero indexed).
 * @return The value of the vector element.
 */
inline double& gmVector2::operator [](int i) 
{
  assert(i == 0 || i == 1);
  return v_[i];
}

/**
 * Returns a const element of the vector.  The method uses an 'assert' call
 * to make sure that the value of the array index is valid.  
 * @param i The vector element to access (zero indexed).
 * @return The const value of the vector element.
 */
inline const double& gmVector2::operator [](int i) const
{
  assert(i == 0 || i == 1);
  return v_[i];
}

// ASSIGNMENT

/**
 * Sets the vector's elements to the given values.
 * @param x The value to set for vector[0].
 * @param y The value to set for vector[1].
 * @return The vector with new values.
 */
inline gmVector2& gmVector2::assign(double x, double y)
{
  v_[0] = x; v_[1] = y;
  return *this;
}

/**
 * Sets the value of the vector to another vector.
 * @param v The other vector to use for setting values.
 * @return The vector with new values.
 */
inline gmVector2& gmVector2::operator =(const gmVector2& v)
{
  v_[0] = v[0]; v_[1] = v[1];
  return *this;
}

// MATH OPERATORS

/**
 * Addition/assignment operator.
 * @param v The vector to add to the current vector.
 * @return The vector += v.
 */
inline gmVector2& gmVector2::operator +=(const gmVector2& v)
{
  v_[0] += v[0]; v_[1] += v[1];
  return *this;
}

/**
 * Subtraction/assignment operator.
 * @param v The vector to subtraction from the current vector.
 * @return The vector -= v.
 */
inline gmVector2& gmVector2::operator -=(const gmVector2& v)
{
  v_[0] -= v[0]; v_[1] -= v[1];
  return *this;
}

/**
 * Multiplication/assignment operator.
 * @param v The vector to multiply with the current vector.
 * @return The vector *= v.
 */
inline gmVector2& gmVector2::operator *=(double c)
{
  v_[0] *= c; v_[1] *= c;
  return *this;
}

/**
 * Division/assignment operator.
 * @param v The vector to divide by the current vector.
 * @return The vector /= v.
 * @note Modified code to be more efficient by inverting c once and
 * multiplying each term.
 */
inline gmVector2& gmVector2::operator /=(double c)
{
  assert(!gmIsZero(c));
  double ci = 1.0/c;
  v_[0] *= ci; v_[1] *= ci;
  return *this;
}

/**
 * Addition operator.  Adds two vectors together.
 * @param v The vector to add to the current vector.
 * @return A new vector resulting from the addition of two vectors.
 */
inline gmVector2 gmVector2::operator +(const gmVector2& v) const
{
  return gmVector2(v_[0] + v[0], v_[1] + v[1]);
}

/**
 * Subtraction operator.  Subtracts one vector from another.
 * @param v The vector to subtract from the current vector.
 * @return A new vector resulting from the difference of two vectors.
 */
inline gmVector2 gmVector2::operator -(const gmVector2& v) const
{
  return gmVector2(v_[0] - v[0], v_[1] - v[1]);
}

/**
 * Unary negation operator. 
 * @return A new vector with all elements negated.
 */
inline gmVector2 gmVector2::operator -() const
{
  return gmVector2(-v_[0], -v_[1]);
}

/**
 * Scalar multiplication operator.  Multiplies a vector with a double.
 * @param c The double to multiply with the vector.
 * @return A new vector resulting from the multiplication of c with each
 * element of the vector.
 */
inline gmVector2 gmVector2::operator *(double c) const
{
  return gmVector2(v_[0] * c, v_[1] * c);
}

/**
 * Scalar division operator.  Divides a vector by a double.
 * @param c The double to divide through the vector.
 * @return A new vector resulting from the division of c with each
 * element of the vector.
 * @note Modified code to be more efficient by inverting c once and
 * multiplying each term.
 */
inline gmVector2 gmVector2::operator /(double c) const
{
  assert(!gmIsZero(c));
  double ci = 1.0/c;
  return gmVector2(v_[0] * ci, v_[1] * ci);
}

/**
 * Scalar multiplication operator.  Multiplies a vector with a double.
 * @param c The double to multiply with the vector.
 * @param v The vector to multiply with the double.
 * @return A new vector resulting from the multiplication of c with each
 * element of the vector v.
 */
inline gmVector2 operator *(double c, const gmVector2& v)
{
  return gmVector2(c * v[0], c * v[1]);
}

/**
 * Equivalence operator.  Checks to see if two vectors are equal.
 * @param v The vector to check for equivalence.
 * @return True if all elements are fuzzy equal, false otherwise.
 */
inline bool gmVector2::operator ==(const gmVector2& v) const
{
  return (gmFuzEQ(v_[0], v[0]) && gmFuzEQ(v_[1], v[1]));
}

/**
 * Non-equivalence operator.  Checks to see if two vectors are not equal.
 * @param v The vector to check for non-equivalence.
 * @return False if all elements are fuzzy equal, true otherwise.
 */
inline bool gmVector2::operator !=(const gmVector2& v) const
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
inline gmVector2& gmVector2::clamp(double lo, double hi)
{
  gmClamp(v_[0], lo, hi); gmClamp(v_[1], lo, hi);
  return *this;
}

/**
 * Returns the length of the vector.
 * @return The length of the vector.
 */
inline double gmVector2::length() const
{
  return sqrt(gmSqr(v_[0]) + gmSqr(v_[1]));
}

/**
 * Returns the length of the vector squared.
 * @return The length of the vector squared.
 */
inline double gmVector2::lengthSquared() const
{
  return gmSqr(v_[0]) + gmSqr(v_[1]);
}

/**
 * Returns a normalized copy of the vector.  This method tests to make sure
 * the length of the vector is not zero.  If the length of the vector is
 * zero, it simply returns the vector without doing the normalization.
 * @return A new vector with length = 1.0 unless the length of the vector
 * was zero.
 */
inline gmVector2& gmVector2::normalize()
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
inline gmVector2& gmVector2::normalizeNonZero()
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
inline void gmVector2::copyTo(float f[2]) const
{
  f[0] = (float)v_[0]; f[1] = (float)v_[1];
}

/**
 * Puts the elements of the vector in an array of doubles.  This is useful
 * if you need to pass the vector as a standard C array.
 * @param f An array of doubles to copy the vector to.
 */
inline void gmVector2::copyTo(double f[2]) const
{
  f[0] = v_[0]; f[1] = v_[1];
}

/**
 * Returns the distance between two vectors.  
 * @param v1 The first vector (point).
 * @param v2 The second vector (point).
 * @return The distance between the two vectors (points).
 */
inline double distance(const gmVector2& v1, const gmVector2& v2)
{
  return sqrt(gmSqr(v1[0] - v2[0]) + gmSqr(v1[1] - v2[1]));
}

/**
 * Returns the distance between two vectors squared.  
 * @param v1 The first vector (point).
 * @param v2 The second vector (point).
 * @return The distance squared between the two vectors (points).
 */
inline double distanceSquared(const gmVector2& v1, const gmVector2& v2)
{
  return gmSqr(v1[0] - v2[0]) + gmSqr(v1[1] - v2[1]);
}

/** 
 * Returns the dot product of two vectors.
 * @param v1 The first vector (point).
 * @param v2 The second vector (point).
 * @return v1 dot v2.
 */
inline double dot(const gmVector2& v1, const gmVector2& v2)
{
  return v1[0] * v2[0] + v1[1] * v2[1];
}

/** 
 * Returns a linear interpolation between two vectors. 
 * @param f The number at which to evaluate the linear interpolation.
 * @param v1 The vector at which f = 0.
 * @param v2 The vector at which f = 1.
 * @return The vector representing the linear interpolation between v1 and
 * v2 evaluated at f.
 */
inline gmVector2 lerp(double f, const gmVector2& v1, const gmVector2& v2)
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
inline std::ostream & operator << ( std::ostream& os, const gmVector2& v)
{
  os << "< " << v[0] << " " << v[1] << " >";
  return os;
}

#endif

