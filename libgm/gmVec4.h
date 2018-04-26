/**
 * This file contains functions for 4-element vectors.
 * @file  gmVec4.h
 * @date  15 June 1994
 * @author Ferdi Scheepers
 * @author Stephen F May
 */

#ifndef GMVECTOR4_H
#define GMVECTOR4_H

#include "gmUtils.h"

/**
 * The gmVector4 class provides methods and operators for creating and
 * working with vectors containing four elements.  You can access the
 * elements of the vector as if the vector was an array.
 * @nosubgrouping
 */
class gmVector4 {

protected:
  double v_[4];   ///< Private storage for the elements

public:
  /** @name Constructors
   * @{
   */
  /// Default constructor.
  gmVector4();
  /// Copy constructor.
  gmVector4(const gmVector4&);
  /// Constructs a new vector object using four doubles.
  gmVector4(double, double, double, double);
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
  gmVector4& assign(double, double, double, double);
  /// Sets the value of the vector to another vector.
  gmVector4& operator =(const gmVector4&);
  /** @} end of assignment */

  /** @name Math Operators
   * @{
   */
  /// Addition/assignment operator.
  gmVector4& operator +=(const gmVector4&);
  /// Subtraction/assignment operator.
  gmVector4& operator -=(const gmVector4&);
  /// Multiplication/assignment operator.
  gmVector4& operator *=(double);
  /// Division/assignment operator.
  gmVector4& operator /=(double);

  /// Addition operator.
  gmVector4 operator +(const gmVector4&) const;
  /// Subtraction operator.
  gmVector4 operator -(const gmVector4&) const;
  /// Unary negation operator.
  gmVector4 operator -() const;
  /// Scalar multiplication operator.
  gmVector4 operator *(double) const;
  /// Scalar division operator.
  gmVector4 operator /(double) const;

  /// Scalar multiplication operator.
friend gmVector4 operator *(double, const gmVector4&);

  /// Equivalence operator.
  bool operator ==(const gmVector4&) const;
  /// Non-equivalence operator.
  bool operator !=(const gmVector4&) const;
  /** @} end of math operations */

  /** @name Operations
   * @{
   */
  /// Clamps each element of the vector to a range.
  gmVector4& clamp(double, double);
  /// Returns the length of the vector.
  double length() const;
  /// Returns the length of the vector squared.
  double lengthSquared() const;
  /// Returns a normalized copy of the vector.
  gmVector4& normalize();
  /// Returns a normalized copy of the vector.
  gmVector4& normalizeNonZero();

  /// Puts the elements of the vector in an array of floats.
  void copyTo(float [4]) const;
  /// Puts the elements of the vector in an array of doubles.
  void copyTo(double [4]) const;

  /// Returns the distance between two vectors.
friend double distance(const gmVector4&, const gmVector4&);
  /// Returns the distance between two vectors squared.
friend double distanceSquared(const gmVector4&, const gmVector4&);
  /// Returns the dot product of two vectors.
friend double dot(const gmVector4&, const gmVector4&);
  /// Returns a linear interpolation between two vectors.
friend gmVector4 lerp(double, const gmVector4&, const gmVector4&);
  /** @} end operations */

  /** @name Output
   * @{
   */
  /// Output stream concatenation operator.
friend std::ostream & operator << ( std::ostream &, const gmVector4 & );
  /** @} end output */
};

// CONSTRUCTORS

/**
 * Default constructor. Sets all elements of the vector to be 0.0.
 */
inline gmVector4::gmVector4()
{
  v_[0] = v_[1] = v_[2] = v_[3] = 0;
}

/**
 * Copy constructor.  Creates a new vector by using the values of another
 * vector.
 * @param v The vector to copy.
 */
inline gmVector4::gmVector4(const gmVector4& v)
{
  v_[0] = v.v_[0]; v_[1] = v.v_[1]; v_[2] = v.v_[2]; v_[3] = v.v_[3];
}

/**
 * Constructs a new vector object using four doubles.
 * @param x The value to set for vector[0].
 * @param y The value to set for vector[1].
 * @param z The value to set for vector[2].
 * @param a The value to set for vector[3].
 */
inline gmVector4::gmVector4(double x, double y, double z, double a)
{
  v_[0] = x; v_[1] = y; v_[2] = z; v_[3] = a;
}

// ARRAY ACCESS

/**
 * Returns an element of the vector.  The method uses an 'assert' call to
 * make sure that the value of the array index is valid.
 * @param i The vector element to access (zero indexed).
 * @return The value of the vector element.
 */
inline double& gmVector4::operator [](int i) 
{
  assert(i == 0 || i == 1 || i == 2 || i == 3);
  return v_[i];
}

/**
 * Returns a const element of the vector.  The method uses an 'assert' call
 * to make sure that the value of the array index is valid.  
 * @param i The vector element to access (zero indexed).
 * @return The const value of the vector element.
 */
inline const double& gmVector4::operator [](int i) const
{
  assert(i == 0 || i == 1 || i == 2 || i == 3);
  return v_[i];
}

// ASSIGNMENT

/**
 * Sets the vector's elements to the given values.
 * @param x The value to set for vector[0].
 * @param y The value to set for vector[1].
 * @param z The value to set for vector[2].
 * @param a The value to set for vector[3].
 * @return The vector with new values.
 */
inline gmVector4& gmVector4::assign(double x, double y, double z, double a)
{
  v_[0] = x; v_[1] = y; v_[2] = z; v_[3] = a;
  return *this;
}

/**
 * Sets the value of the vector to another vector.
 * @param v The other vector to use for setting values.
 * @return The vector with new values.
 */
inline gmVector4& gmVector4::operator =(const gmVector4& v)
{
  v_[0] = v[0]; v_[1] = v[1]; v_[2] = v[2]; v_[3] = v[3];
  return *this;
}

// MATH OPERATORS

/**
 * Addition/assignment operator.
 * @param v The vector to add to the current vector.
 * @return The vector += v.
 */
inline gmVector4& gmVector4::operator +=(const gmVector4& v)
{
  v_[0] += v[0]; v_[1] += v[1]; v_[2] += v[2]; v_[3] += v[3];
  return *this;
}

/**
 * Subtraction/assignment operator.
 * @param v The vector to subtraction from the current vector.
 * @return The vector -= v.
 */
inline gmVector4& gmVector4::operator -=(const gmVector4& v)
{
  v_[0] -= v[0]; v_[1] -= v[1]; v_[2] -= v[2]; v_[3] -= v[3];
  return *this;
}

/**
 * Multiplication/assignment operator.
 * @param v The vector to multiply with the current vector.
 * @return The vector *= v.
 */
inline gmVector4& gmVector4::operator *=(double c)
{
  v_[0] *= c; v_[1] *= c; v_[2] *= c; v_[3] *= c;
  return *this;
}

/**
 * Division/assignment operator.
 * @param v The vector to divide by the current vector.
 * @return The vector /= v.
 * @note Modified code to be more efficient by inverting c once and
 * multiplying each term.
 */
inline gmVector4& gmVector4::operator /=(double c)
{
  assert(!gmIsZero(c));
  
  v_[0] *= c; v_[1] *= c; v_[2] *= c; v_[3] *= c;
  return *this;
}

/**
 * Addition operator.  Adds two vectors together.
 * @param v The vector to add to the current vector.
 * @return A new vector resulting from the addition of two vectors.
 */
inline gmVector4 gmVector4::operator +(const gmVector4& v) const
{
  return gmVector4(v_[0] + v[0], v_[1] + v[1], v_[2] + v[2], v_[3] + v[3]);
}

/**
 * Subtraction operator.  Subtracts one vector from another.
 * @param v The vector to subtract from the current vector.
 * @return A new vector resulting from the difference of two vectors.
 */
inline gmVector4 gmVector4::operator -(const gmVector4& v) const
{
  return gmVector4(v_[0] - v[0], v_[1] - v[1], v_[2] - v[2], v_[3] - v[3]);
}

/**
 * Unary negation operator. 
 * @return A new vector with all elements negated.
 */
inline gmVector4 gmVector4::operator -() const
{
  return gmVector4(-v_[0], -v_[1], -v_[2], -v_[3]);
}

/**
 * Scalar multiplication operator.  Multiplies a vector with a double.
 * @param c The double to multiply with the vector.
 * @return A new vector resulting from the multiplication of c with each
 * element of the vector.
 */
inline gmVector4 gmVector4::operator *(double c) const
{
  return gmVector4(v_[0] * c, v_[1] * c, v_[2] * c, v_[3] * c);
}

/**
 * Scalar division operator.  Divides a vector by a double.
 * @param c The double to divide through the vector.
 * @return A new vector resulting from the division of c with each
 * element of the vector.
 * @note Modified code to be more efficient by inverting c once and
 * multiplying each term.
 */
inline gmVector4 gmVector4::operator /(double c) const
{
  assert(!gmIsZero(c));
  double ci = 1.0/c;
  return gmVector4(v_[0] * ci, v_[1] * ci, v_[2] * ci, v_[3] * ci);
}

/**
 * Scalar multiplication operator.  Multiplies a vector with a double.
 * @param c The double to multiply with the vector.
 * @param v The vector to multiply with the double.
 * @return A new vector resulting from the multiplication of c with each
 * element of the vector v.
 */
inline gmVector4 operator *(double c, const gmVector4& v)
{
  return gmVector4(c * v[0], c * v[1], c * v[2], c * v[3]);
}

/**
 * Equivalence operator.  Checks to see if two vectors are equal.
 * @param v The vector to check for equivalence.
 * @return True if all elements are fuzzy equal, false otherwise.
 */
inline bool gmVector4::operator ==(const gmVector4& v) const
{
  return (gmFuzEQ(v_[0], v[0]) && gmFuzEQ(v_[1], v[1]) &&
          gmFuzEQ(v_[2], v[2]) && gmFuzEQ(v_[3], v[3]));
}

/**
 * Non-equivalence operator.  Checks to see if two vectors are not equal.
 * @param v The vector to check for non-equivalence.
 * @return False if all elements are fuzzy equal, true otherwise.
 */
inline bool gmVector4::operator !=(const gmVector4& v) const
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
inline gmVector4& gmVector4::clamp(double lo, double hi)
{
  gmClamp(v_[0], lo, hi); gmClamp(v_[1], lo, hi);
  gmClamp(v_[2], lo, hi); gmClamp(v_[3], lo, hi);
  return *this;
}

/**
 * Returns the length of the vector.
 * @return The length of the vector.
 */
inline double gmVector4::length() const
{
  return sqrt(gmSqr(v_[0]) + gmSqr(v_[1]) + gmSqr(v_[2]) + gmSqr(v_[3]));
}

/**
 * Returns the length of the vector squared.
 * @return The length of the vector squared.
 */
inline double gmVector4::lengthSquared() const
{
  return gmSqr(v_[0]) + gmSqr(v_[1]) + gmSqr(v_[2]) + gmSqr(v_[3]);
}

/**
 * Returns a normalized copy of the vector.  This method tests to make sure
 * the length of the vector is not zero.  If the length of the vector is
 * zero, it simply returns the vector without doing the normalization.
 * @return A new vector with length = 1.0 unless the length of the vector
 * was zero.
 */
inline gmVector4& gmVector4::normalize()
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
inline gmVector4& gmVector4::normalizeNonZero()
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
inline void gmVector4::copyTo(float f[4]) const
{
  f[0] = (float)v_[0]; f[1] = (float)v_[1]; f[2] = (float)v_[2]; f[3] = (float)v_[3];
}

/**
 * Puts the elements of the vector in an array of doubles.  This is useful
 * if you need to pass the vector as a standard C array.
 * @param f An array of doubles to copy the vector to.
 */
inline void gmVector4::copyTo(double f[4]) const
{
  f[0] = v_[0]; f[1] = v_[1]; f[2] = v_[2]; f[3] = v_[3];
}

/**
 * Returns the distance between two vectors.  
 * @param v1 The first vector (point).
 * @param v2 The second vector (point).
 * @return The distance between the two vectors (points).
 */
inline double distance(const gmVector4& v1, const gmVector4& v2)
{
  return sqrt(gmSqr(v1[0] - v2[0]) + gmSqr(v1[1] - v2[1]) +
              gmSqr(v1[2] - v2[2]) + gmSqr(v1[3] - v2[3]));
}

/**
 * Returns the distance between two vectors squared.  
 * @param v1 The first vector (point).
 * @param v2 The second vector (point).
 * @return The distance squared between the two vectors (points).
 */
inline double distanceSquared(const gmVector4& v1, const gmVector4& v2)
{
  return gmSqr(v1[0] - v2[0]) + gmSqr(v1[1] - v2[1]) +
         gmSqr(v1[2] - v2[2]) + gmSqr(v1[3] - v2[3]);
}

/** 
 * Returns the dot product of two vectors.
 * @param v1 The first vector (point).
 * @param v2 The second vector (point).
 * @return v1 dot v2.
 */
inline double dot(const gmVector4& v1, const gmVector4& v2)
{
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3];
}

/** 
 * Returns a linear interpolation between two vectors. 
 * @param f The number at which to evaluate the linear interpolation.
 * @param v1 The vector at which f = 0.
 * @param v2 The vector at which f = 1.
 * @return The vector representing the linear interpolation between v1 and
 * v2 evaluated at f.
 */
inline gmVector4 lerp(double f, const gmVector4& v1, const gmVector4& v2)
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
inline std::ostream & operator << ( std::ostream& os, const gmVector4& v)
{
  os << "< " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << " >";
  return os;
}

#endif 

