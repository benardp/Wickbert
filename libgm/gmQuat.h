/**
 * This file contains functions for Quaternions.
 * @file  gmQuat.h
 * @date  4 May 2003
 * @author Jerry O. Talton III
 */

#ifndef GMQUAT_H
#define GMQUAT_H

#include "gmUtils.h"

class gmVector3;
class gmMat4;

/**
 * The gmQuat class provides methods and operators for creating and
 * working with Quaternions.
 * @nosubgrouping
 */
class gmQuat {

protected:
  gmVector3 v;  ///< Private storage for the vector elements
  double s;  ///< Private storage for the scalar

public:
  /** @name Constructors
   * @{
   */
  /// Default constructor.
  gmQuat();
  /// Copy constructor.
  gmQuat(const gmQuat&);
  /// Constructs a new quaternion object using a vector and a double
  gmQuat(const gmVector3&, double); 
  /// Constructs a new quaternion object using four doubles.
  gmQuat(double, double, double, double);
  /** @} end of constructors */

  /** @name Access methods
   * @{
   */
  /// Returns the vector element.
  gmVector3& vector();
  /// Returns a const vector.
  const gmVector3& vector() const;
  /// Returns the scalar element.
  double scalar();
  /// Returns a const scalar.
  const double scalar() const;
  /** @} end of Access methods*/

  /** @name Assignment
   * @{
   */
  /// Sets the quaternion's elements to the given values.
  gmQuat& assign(double, double, double, double);
  /// Sets the value of the quaternion to another quaternion.
  gmQuat& operator =(const gmQuat&);
  /** @} end of assignment */

  /** @name Math Operators
   * @{
   */
  /// Addition operator.
  gmQuat operator +(const gmQuat&) const;
  /// Subtraction operator.
  gmQuat operator -(const gmQuat&) const;
  /// Multiplication operator.
  gmQuat operator *(const gmQuat&) const;
  /// Scalar multiplication operator.
  gmQuat operator *(const double) const;
  /// Scalar division operator.
  gmQuat operator /(const double) const;

  /// Addition/assignment operator.
  gmQuat& operator +=(const gmQuat&);
  /// Subtraction/assignment operator.
  gmQuat& operator -=(const gmQuat&);
  /// Multiplication/assignment operator.
  gmQuat& operator *=(const gmQuat&);
  /// Multiplication by scalar/assignment operator.
  gmQuat& operator *=(const double);
  /// Division by scalar/assignment operator.
  gmQuat& operator /=(const double);

  /// Scalar multiplication operators.
  friend gmQuat operator *(double, const gmQuat&); 

  /// Equivalence operator.
  bool operator ==(const gmQuat&) const;
  /// Non-equivalence operator.
  bool operator !=(const gmQuat&) const;
  /** @} end of math operations */

  /** @name Operations
   * @{
   */
  /// Returns the norm of the quaternion.
  double norm() const;
  /// Returns a normalized copy of the quaternion.
  gmQuat& normalize();
  /// Returns the conjugate of the quaternion.
  gmQuat conjugate() const;
  /// Returns the inverse of the quaternion.
  gmQuat inverse() const;
  /// Returns the exponentiation of the quaternion
  gmQuat exp() const;
  /// Takes the log of a quaternion.
  gmQuat log() const;

  /// Returns a rotation about the vector
  friend gmQuat axis_rot_quat(const gmVector3&, double);

  /// Returns the quaternion that maps the first vector onto the second
  friend gmQuat vec_to_vec_quat(const gmVector3&, const gmVector3&);

  /// Returns the result of rotating a vector by a quaternion
  friend gmVector3 rotate_by_quat(const gmVector3&, const gmQuat&);

  /// Returns the rotation matrix corresponding to the current quaternion
  friend gmMatrix4 quat_to_matrix(const gmQuat& q);

  /// Puts the elements of the quaternion in an array of floats.
  void copyTo(float [4]) const;
  /// Puts the elements of the quaternion in an array of doubles.
  void copyTo(double [4]) const;

  /** @} end operations */

  /** @name Output
   * @{
   */
  /// Output stream concatenation operator.
  friend std::ostream & operator << ( std::ostream&, const gmQuat&);
  /** @} end output */
};

// CONSTRUCTORS

/**
 * Default constructor. 
 */
inline gmQuat::gmQuat()
{
	v[0] = v[1] = v[2] = 0.0; 
	s = 1.0;
}

/**
 * Copy constructor.  Creates a new quaternion by using the values of another
 * quaternion.
 * @param q The quaternion to copy.
 */
inline gmQuat::gmQuat(const gmQuat& q)
{
  v = q.v;
  s = q.s;
}

/**
 * Creates a new quaternion by using the values of a vector
 * and a scalar.
 * @param vec The vector to copy.
 * @param sc The scalar to copy.
 */
inline gmQuat::gmQuat(const gmVector3& vec, double sc)
{
  v = vec;
  s = sc;
}

/**
 * Constructs a new vector object using three doubles.
 * @param x The value to set for v[0].
 * @param y The value to set for v[1].
 * @param z The value to set for v[2].
 * @param w The value to be set for s.
 */
inline gmQuat::gmQuat(double x, double y, double z, double w)
{
  v[0] = x; v[1] = y; v[2] = z;
  s = w;
}


// ACCESS METHODS

/**
 * Returns the vector element of the quaternion.
 * @return The value of the vector element.
 */
inline gmVector3& gmQuat::vector() 
{
	return v;
}

/**
 * Returns the const vector element of the quaternion.  
 * @return The const vector element of the quaternion element.
 */
inline const gmVector3& gmQuat::vector() const
{
  return v;
}

/**
 * Returns the scalar element of the quaternion.
 * @return The value of the scalar element of the quaternion.
 */
inline double gmQuat::scalar() 
{
	return s;
}

/**
 * Returns the const scalar element of the quaternion.
 * @return The const value of the scalar element of the quaternion.
 */
inline const double gmQuat::scalar() const 
{
	return s;
}


// ASSIGNMENT

/**
 * Sets the quaternion's elements to the given values.
 * @param x The value to set for vector[0].
 * @param y The value to set for vector[1].
 * @param z The value to set for vector[2].
 * @param w The value to set for s.
 * @return The vector with new values.
 */
inline gmQuat& gmQuat::assign(double x, double y, double z, double w)
{
  v[0] = x; v[1] = y; v[2] = z;
  s = w;
  return *this;
}

/**
 * Sets the value of the quaternion to another quaternion.
 * @param v The other quaternion to use for setting values.
 * @return The quaternion with new values.
 */
inline gmQuat& gmQuat::operator =(const gmQuat& q)
{
  v = q.v; 
  s = q.s;
  return *this;
}


// MATH OPERATORS

/**
 * Addition operator.  Adds two quaternions together.
 * @param q The quaternions to add to the current vector.
 * @return A new quaternion resulting from the addition of two quaternions.
 */
inline gmQuat gmQuat::operator +(const gmQuat& q) const
{
  return gmQuat(v + q.v, s + q.s);
}

/**
 * Subtraction operator.  Subtracts one quaternion from another.
 * @param q The quaternion to subtract from the current quaternion.
 * @return A new quaternion resulting from the difference of two quaternions.
 */
inline gmQuat gmQuat::operator -(const gmQuat& q) const
{
  return gmQuat(v - q.v, s - q.s);
}

/**
 * Multiplication operator.  Multiplies a quaternion with a quaternion.
 * @param q The double to multiply with the quaternion.
 * @return A new quaternion resulting from the multiplication of c with each
 * element of the quaternion.
 */
inline gmQuat gmQuat::operator *(const gmQuat& q) const
{
  return gmQuat(cross(v, q.v) + q.s * v + s * q.v, s * q.s - dot(v, q.v));
}

/**
 * Scalar multiplication operator.  Multiplies a quaternion with a double.
 * @param c The double to multiply with the quaternion.
 * @return A new quaternion resulting from the multiplication of c with each
 * element of the quaternion.
 */
inline gmQuat gmQuat::operator *(const double c) const
{
  return gmQuat(v * c, s * c);
}

/**
 * Scalar division operator.  Divides a quaternion by a double.
 * @param c The double to divide through the quaternion.
 * @return A new quaternion resulting from the division of c with each
 * element of the quaternion.
 */
inline gmQuat gmQuat::operator /(double c) const
{
  assert(!gmIsZero(c));
  double ci = 1.0 / c;
  return gmQuat(v * ci, s * ci);
}


/**
 * Addition/assignment operator.
 * @param q The quaternion to add to the current quaternion.
 * @return The vector += q.
 */
inline gmQuat& gmQuat::operator +=(const gmQuat& q)
{
  v += q.v;  s += q.s;
  return *this;
}

/**
 * Subtraction/assignment operator.
 * @param q The quaternion to subtract from the current quaternion.
 * @return The quaternion -= q.
 */
inline gmQuat& gmQuat::operator -=(const gmQuat& q)
{
  v -= q.v;  s -= q.s;
  return *this;
}

/**
 * Multiplication/assignment operator.
 * @param q The quaternion to multiply with the current quaternion.
 * @return The quaternion *= q.
 */
inline gmQuat& gmQuat::operator *=(const gmQuat& q)
{
  v = cross(v, q.v) + q.s * v + s * q.v;
  s = s * q.s - dot(v, q.v);
  return *this;
}

/**
 * Multiplication/assignment operator with a scalar.
 * @param c The scalar to multiply with the current quaternion.
 * @return The quaternion *= c.
 */
inline gmQuat& gmQuat::operator *=(const double c)
{
  v *= c;  s *= c;
  return *this;
}

/**
 * Division/assignment operator with a scalar.
 * @param c The scalar to divide through the current quaternion.
 * @return The quaternion /= c.
 */
inline gmQuat& gmQuat::operator /=(double c)
{
  assert(!gmIsZero(c));
  double ci = 1.0 / c;
  v *= ci; s *= ci;
  return *this;
}


/**
 * Scalar multiplication operator.  Multiplies a quaternion with a double.
 * @param c The double to multiply with the quaternion.
 * @param q The quaternion to multiply with the double.
 * @return A new quaternion resulting from the multiplication of c with each
 * element of the quaternion q.
 */
inline gmQuat operator *(double c, const gmQuat& q)
{
  return gmQuat(c * q.v, c * q.s);
}


/**
 * Equivalence operator.  Checks to see if two quaternions are equal.
 * @param q The quaternion to check for equivalence.
 * @return True if all elements are fuzzy equal, false otherwise.
 */
inline bool gmQuat::operator ==(const gmQuat& q) const
{
  return
    (gmFuzEQ(v[0], q.v[0]) && gmFuzEQ(v[1], q.v[1]) && gmFuzEQ(v[2], q.v[2]) && gmFuzEQ(s, q.s));
}

/**
 * Non-equivalence operator.  Checks to see if two quaternions are not equal.
 * @param q The quaternion to check for non-equivalence.
 * @return False if all elements are fuzzy equal, true otherwise.
 */
inline bool gmQuat::operator !=(const gmQuat& q) const
{
  return (!(*this == q));
}


// OPERATIONS

/**
 * Returns the norm of the quaternion.
 * @return The norm of the quaternion.
 */
inline double gmQuat::norm() const
{
  return v.lengthSquared() + gmSqr(s);
}

/**
 * Conjugate operation.
 * @return The conjugate of the quaternion.
 */
inline gmQuat gmQuat::conjugate() const
{
  return gmQuat(-v, s);
}

/**
 * Returns a normalized copy of the quaternion.  This method tests to make sure
 * the length of the vector is not zero.  If the length of the vector is
 * zero, it simply returns the vector without doing the normalization.
 * @return A new vector with length = 1.0 unless the length of the vector
 * was zero.
 */
inline gmQuat& gmQuat::normalize()
{
  double n = norm();
  if (!gmIsZero(n))
    *this = conjugate() / n;
  return *this;
}

/**
 * Inverse operation.
 * @return The inverse of the quaternion.
 */
inline gmQuat gmQuat::inverse() const
{
  return conjugate() / norm();
}

/**
 * Exponential operation, only well defined if the scalar part is 0.
 * @return The exponential of the quaternion.
 */
inline gmQuat gmQuat::exp() const
{
   assert(gmIsZero(s));

   double theta = v.length();
   double c = cos(theta);

   if( theta > gmEPSILON )
   {
     double si = sin(theta) / theta;
     return gmQuat(si * v, c);
   }
   else return gmQuat(v, c);
}

/**
 * Logarithm operation, only well defined if the quaternion is a UNIT
 * quaternion.
 * @return The logarithm of the quaternion.
 */
inline gmQuat gmQuat::log() const
{
	double scale = v.length();

	assert(gmFuzEQ(scale, 1.0));

    double theta = atan2(scale, s);

    if(scale > 0.0) scale = theta / scale;

    return gmQuat(scale * v, 0.0);
}


/**
 * Returns a quaternion representing a rotation of phi degrees about
 * the vector v.
 * param v The vector to rotate about.
 * param phi The degrees to rotate.
 * @return A quaternion representing the rotation.
 */
inline gmQuat axis_rot_quat(const gmVector3& v, double phi) 
{
	gmVector3 u = v;
    u = u.normalize();

    double s = sin(phi / 2.0);
    return gmQuat(s * u[0], s * u[1], s * u[2], cos(phi / 2.0));
}

/**
 * Returns a quaternion representing a rotation that will map
 * the vector u onto the vector v through the shortest rotation.
 * param u The vector to rotate.
 * param v The vector to rotate onto.
 * @return A quaternion representing the rotation.
 */
inline gmQuat vec_to_vec_quat(const gmVector3& u, const gmVector3& v) 
{
	gmVector3 w = cross(u, v);

	// If the vectors are colinear, blow up
	if (w == gmVector3(0.0, 0.0, 0.0)) return gmQuat();

	return axis_rot_quat(w, acos(dot(u, v)));   
}

/**
 * Returns a quaternion representing a rotation that will map
 * the vector u onto the vector v through the shortest rotation.
 * param v The vector to rotate.
 * param q The quaternion representing the rotation.
 * @return The rotated vector.
 */
inline gmVector3 rotate_by_quat(const gmVector3& v, const gmQuat& q) 
{
	gmQuat q_vec = gmQuat(v, 0.0);
	gmQuat q_inv = q.inverse();

	return (q * q_vec * q_inv).vector();
	
}

/**
 * Returns the rotation matrix represented by the quaternion.
 * param q The quaternion.
 * @return A rotation matrix representing the quaternion.
 */
inline gmMatrix4 quat_to_matrix(const gmQuat& q)
{
    gmMatrix4 M;

    const double x = q.v[0];
    const double y = q.v[1];
    const double z = q.v[2];
    const double w = q.s;
    const double s = 2 / q.norm();

    M[0][0]=1-s*(y*y+z*z); M[0][1]=s*(x*y-w*z);   M[0][2]=s*(x*z+w*y);   M[0][3]=0;
    M[1][0]=s*(x*y+w*z);   M[1][1]=1-s*(x*x+z*z); M[1][2]=s*(y*z-w*x);   M[1][3]=0;
    M[2][0]=s*(x*z-w*y);   M[2][1]=s*(y*z+w*x);   M[2][2]=1-s*(x*x+y*y); M[2][3]=0;
    M[3][0]=0;             M[3][1]=0;             M[3][2]=0;             M[3][3]=1;

    return M;
}

/**
 * Puts the elements of the quaternion in an array of floats.  This is useful
 * if you need to pass the quaternion as a standard C array.
 * @param f An array of floats to copy the quaternion to.
 */
inline void gmQuat::copyTo(float f[4]) const
{
  f[0] = (float)v[0]; f[1] = (float)v[1]; f[2] = (float)v[2]; f[3] = (float)s;
}

/**
 * Puts the elements of the quaternion in an array of doubles.  This is useful
 * if you need to pass the quaternion as a standard C array.
 * @param f An array of doubles to copy the quaternion to.
 */
inline void gmQuat::copyTo(double f[4]) const
{
  f[0] = v[0]; f[1] = v[1]; f[2] = v[2]; f[3] = s;
}


// OUTPUT

/**
 * Output stream concatenation operator.  This allows the text output of a
 * quaternion object using the << operator.
 * @param os The output stream to print to.
 * @param q The quaternion to output.
 * @return The resulting output stream.
 */
inline std::ostream & operator << ( std::ostream& os, const gmQuat& q)
{
  os << "[ <" << q.v[0] << " " << q.v[1] << " " << q.v[2] << ">, " << q.s << " ]";
  return os;
}


#endif 
