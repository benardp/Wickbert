/**
 * This file contains the class declaration and implementation for an
 * Infinite object, which adds the values negative infinity and positive
 * infinity to the basic data types int, float, double, and long double.
 * Most of the basic arithmetic and logical comparison operators are
 * provided.  The only thing missing at this point is checking for overflow
 * errors and handling them appropriately.
 *
 * @file Infinite.h
 * @author Terry Fleury
 * @date November 9, 2001
 */

#ifndef _INFINITE_H
#define _INFINITE_H

#ifdef _MSC_VER
# pragma warning( disable: 4661)
#endif

#include <limits.h>
#include <float.h>
#include <math.h>
#include <iostream>
#include "libgm/gm.h"

/**
 * Infinite is a class that allows for values to be stored as negative or
 * positive infinity in addition to standard integer/floating point 
 * values.  This class is needed because there is currently no standard 
 * method of handling math operations when one or both operands of a 
 * math operator (such as +,-,*,/) is either negative infinity or 
 * positive infinity.  This class attempts to handle most of the basic
 * math and comparion operators.  One notable exception: If both of the
 * operands are valid non-infinite numbers and you attempt to do a math
 * operation which would result in overflow, this class does NOT currently
 * catch the condition and set the appropriate infinity value!
 *
 * So, how do you use this class?  Currently, this class will handle
 * only ints, floats, doubles, and long doubles.  You declare a new
 * Infinite instance by passing in the type to the template.  For
 * example, if you want to use "standard" double math, but want to
 * be able to keep track of infinity conditions, you would do:
 *
 * Infinite<double> a;
 *
 * The default constructor would set the value of a to be 0.0.  You can
 * then do many of the typical math operations as if it was a simple
 * double type.  Also, this class has adapters to automatically promote
 * base types to Infinite types and to demote Infinite types to base
 * types.  This allows you to combine base types and Infinite types
 * with the standard math/comparison operators. 
 *
 * There are two class (static) member constants to represent positive
 * infinity (I_INF) and negative infinity (I_NEGINF).  So if you wanted to
 * create a couple of new Infinite instances with positive and negative
 * infinity values, you would do something like this:
 *
 * Infinite<double> my_infinity(Infinite<double>::I_INF);
 * Infinite<double> my_negative_infinity(Infinite<double>::I_NEGINF);
 *
 * Notice that you need to actually declare the type when referencing I_INF
 * and I_NEGINF.  This is so that the interal constants can be set to the
 * largest available number depending on the type.
 *
 * Some of the various operations on one or two infinite valued arguments
 * are not well defined.  In the chart below, you will find the results of
 * the various math and comparion operators.  NOTICE that I have redefined
 * the bitwise operator (^) to act as a power function.  You can still do
 * a.pow(b), and you can also do a^b.  Cool!
 * 
 * <HR>
 * <PRE>
 * +I = positive infinity
 * -I = negative infinity
 *  0 = the value zero 
 *  x = a value that is neither infinite nor zero
 *  # = results that should probably be 'undefined' - maybe change later?
 *
 *                                Mathematical Operations
 *  a   b  -a  e^a a+b  a-b       a*b         a/b               a^b          
 * --- --- --- --- --- ------ ----------- ----------- ----------------------
 * +I  +I  -I  +I  +I  a+(-b)     +I           0                +I          
 * +I  -I           0#            -I           0                 0     
 * -I  -I  +I   0  -I             +I           0                 0
 * -I  +I           0#            -I           0                +I#
 *  0  +I   0   1  +I              0           0                 0
 *  0  -I          -I              0           0                 0
 * +I   0          +I              0          +I#                1
 * -I   0          -I              0          -I#                1
 *  x  +I  -x  e^x +I         (x>0)?+I:-I      0                +I
 *  x  -I          -I         (x>0)?-I:+I      0                 0
 * +I   x          +I         (x>0)?+I:-I (x>0)?+I:-I       (x<0)?0:+I
 * -I   x          -I         (x>0)?-I:+I (x>0)?-I:+I (x<0)?0:(odd(b)?-I:+I)
 *
 *           Logical Operations
 *  a   b  a <= b a >= b  a < b   a > b 
 * --- --- ------ ------ ------- -------
 * +I  +I     T      T   !(a>=b) !(a<=b)
 * +I  -I     F      T  
 * -I  -I     T      T 
 * -I  +I     T      F  
 *  0  +I     T      F  
 *  0  -I     F      T  
 * +I   0     F      T  
 * -I   0     T      F  
 *  x  +I     T      F  
 *  x  -I     F      T  
 * +I   x     F      T  
 * -I   x     T      F  
 * </PRE>
 * <HR>
 *
 * @note Currently, this class does NOT handle overflow conditions when
 *       you do a valid math operation on two non-infinite values
 *       but would end up with a value that is out of range of the
 *       underlying datatype.
 * @todo Check for overflow on +,-,*,/ when the two values are not
 *       I_INF or I_NEGINF (ie. would use a "standard" math operator).
 */
template<class Type>
class Infinite
{
  private:
    Type _val;   ///< The value of the number.

  public:

	static const Type I_INF;
	static const Type I_NEGINF;

    /// @name Constructors @{
    Infinite(Type x = 0);  ///< Default constructor
    /// @} end Constructors 

    /// @name Interrogation @{
    Type val() const;
    bool isInfinite() const;
    bool isPositiveInfinite() const;
    bool isNegativeInfinite() const;
    bool hasInfiniteValue() const;
    operator Type() const;
    /// @} end Interrogation 

    /// @name Function Methods @{
    Infinite fabs() const;
    Infinite operator-() const;
    Infinite exp() const;
    Infinite sqrt() const;
    Infinite pow(const Infinite& x) const;
    Infinite pow(const Type& x) const;
    static Infinite min(const Type& x, const Type& y);
    static Infinite max(const Type& x, const Type& y);
    /// @} end Function Methods 

    /// @name Mathematical Operators @{
    Infinite& operator+=(const Infinite& x);
    Infinite& operator+=(const Type& x);
    Infinite& operator-=(const Infinite& x);
    Infinite& operator-=(const Type& x);
    Infinite& operator*=(const Infinite& x);
    Infinite& operator*=(const Type& x);
    Infinite& operator/=(const Infinite& x);
    Infinite& operator/=(const Type& x);
    Infinite& operator^=(const Infinite& x);
    Infinite& operator^=(const Type& x);
    /// @} end Mathematical Operators
}; 




/// The class constant for positive infinity.
template<class Type>
const Type Infinite<Type>::I_INF = (((Type)0.5 == 0.5) ? 
      ((sizeof(Type) == sizeof(float)) ? FLT_MAX : DBL_MAX) : INT_MAX);

/// The class constant for negative infinity.
template<class Type>
const Type Infinite<Type>::I_NEGINF = (-Infinite<Type>::I_INF);

/// Binary addition for two Infinite objects.
template<class Type>
Infinite<Type> operator+(const Infinite<Type>& x,const Infinite<Type>& y)
{
  Infinite<Type> z = x;
  return (z += y);
}

/// Binary addition for one Infinite type object and one basic type object.
template<class Type>
Infinite<Type> operator+(const Infinite<Type>& x, const Type& y)
{
  Infinite<Type> z = x;
  return (z += y);
}

/// Binary addition for one basic type object and one Infinite type object.
template<class Type>
Infinite<Type> operator+(const Type& x, const Infinite<Type>& y)
{
  Infinite<Type> z = y;
  return (z += x);
}

/// Binary subtraction for two Infinite objects.
template<class Type>
Infinite<Type> operator-(const Infinite<Type>& x,const Infinite<Type>& y)
{
  Infinite<Type> z = x;
  return (z -= y);
}

/// Binary subtraction for one Infinite type and one basic type object.
template<class Type>
Infinite<Type> operator-(const Infinite<Type>& x, const Type& y)
{
  Infinite<Type> z = x;
  return (z -= y);
}

/// Binary subtraction for one basic type and one Infinite type object.
template<class Type>
Infinite<Type> operator-(const Type& x, const Infinite<Type>& y)
{
  Infinite<Type> z = Infinite<Type>(x);
  return (z -= y);
}

/// Binary multiply for two Infinite objects.
template<class Type>
Infinite<Type> operator*(const Infinite<Type>& x,const Infinite<Type>& y)
{
  Infinite<Type> z = x;
  return (z *= y);
}

/// Binary multiply for one Infinite type object and one basic type object.
template<class Type>
Infinite<Type> operator*(const Infinite<Type>& x, const Type& y)
{
  Infinite<Type> z = x;
  return (z *= y);
}
     
/// Binary multiply for one basic type object and one Infinite type object.
template<class Type>
Infinite<Type> operator*(const Type& x, const Infinite<Type>& y)
{
  Infinite<Type> z = y;
  return (z *= x);
}

/// Binary divide for two Infinite objects.
template<class Type>
Infinite<Type> operator/(const Infinite<Type>& x,const Infinite<Type>& y)
{
  Infinite<Type> z = x;
  return (z /= y);
}

/// Binary divide for one Infinite type object and one basic type object.
template<class Type>
Infinite<Type> operator/(const Infinite<Type>& x, const Type& y)
{
  Infinite<Type> z = x;
  return (z /= y);
}

/// Binary divide for one basic type object and one Infinite type object.
template<class Type>
Infinite<Type> operator/(const Type& x, const Infinite<Type>& y)
{
  Infinite<Type> z = Infinite<Type>(x);
  return (z /= y);
}

/// Binary power op for two Infinite objects.
template<class Type>
Infinite<Type> operator^(const Infinite<Type>& x,const Infinite<Type>& y)
{
  Infinite<Type> z = x;
  return (z ^= y);
}

/// Binary power op for one Infinite type object and one basic type object.
template<class Type>
Infinite<Type> operator^(const Infinite<Type>& x, const Type& y)
{
  Infinite<Type> z = x;
  return (z ^= y);
}

/// Binary power op for one basic type object and one Infinite type object.
template<class Type>
Infinite<Type> operator^(const Type& x, const Infinite<Type>& y)
{
  Infinite<Type> z = Infinite<Type>(x);
  return (z ^= y);
}

/// Binary Less-than-equal operator for two Infinite types.
template<class Type>
bool operator<=(const Infinite<Type>& x, const Infinite<Type>& y) 
{
  if (x.isNegativeInfinite() || y.isInfinite())
    return true;
  else if (x.hasInfiniteValue() || y.hasInfiniteValue())
    return false;
  else
    return gmFuzLEQ(x.val(),y.val());
}

/// Binary Less-than-equal operator for one Infinite and one basic object.
template<class Type>
bool operator<=(const Infinite<Type>& x, const Type& y) 
{ 
  return (x <= Infinite<Type>(y));
}

/// Binary Less-than-equal operator for one basic and one Infinite object.
template<class Type>
bool operator<=(const Type& x, const Infinite<Type>& y) 
{ 
  return (Infinite<Type>(x) <= y);
}

/// Binary Greater-than-equal operator for two Infinite types.
template<class Type>
bool operator>=(const Infinite<Type>& x, const Infinite<Type>& y) 
{ 
  if (x.isInfinite() || y.isNegativeInfinite())
    return true;
  else if (x.hasInfiniteValue() || y.hasInfiniteValue())
    return false;
  else 
    return gmFuzGEQ(x.val(),y.val());
}

/// Binary Greater-than-equal operator for one Infinite, one basic object.
template<class Type>
bool operator>=(const Infinite<Type>& x, const Type& y) 
{ 
  return (x >= Infinite<Type>(y));
}

/// Binary Greater-than-equal operator for one basic, one Infinite object.
template<class Type>
bool operator>=(const Type& x, const Infinite<Type>& y) 
{ 
  return (Infinite<Type>(x) >= y);
}

/// Binary Less-than operator for two Infinite types.
template<class Type>
bool operator<(const Infinite<Type>& x, const Infinite<Type>& y)
{
  return !(x >= y);
}

/// Binary Less-than operator for one Infinite and one basic type object.
template<class Type>
bool operator<(const Infinite<Type>& x, const Type& y) 
{ 
  return (x < Infinite<Type>(y)); 
}

/// Binary Less-than operator for one basic and one Infinite type object.
template<class Type>
bool operator<(const Type& x, const Infinite<Type>& y) 
{
  return (Infinite<Type>(x) < y); 
}

/// Binary Greater-than operator for two Infinite type objects.
template<class Type>
bool operator>(const Infinite<Type>& x, const Infinite<Type>& y) 
{ 
  return !(x <= y);
}

/// Binary Greater-than operator for one Infinite and one basic type object.
template<class Type>
bool operator>(const Infinite<Type>& x, const Type& y) 
{ 
  return (x > Infinite<Type>(y)); 
}

/// Binary Greater-than operator for one basic and one Infinite type object.
template<class Type>
bool operator>(const Type& x, const Infinite<Type>& y) 
{ 
  return (Infinite<Type>(x) > y);
}

/// Binary Equality comparison for two Infinite type objects.
template<class Type>
bool operator==(const Infinite<Type>& x, const Infinite<Type>& y) 
{
  return gmFuzEQ(x.val(),y.val());
}

/// Binary Equality comparison for one Infinite and one basic type object.
template<class Type>
bool operator==(const Infinite<Type>& x, const Type& y) 
{ 
  return (x == Infinite<Type>(y)); 
}

/// Binary Equality comparison for one basic and one Infinite type object.
template<class Type>
bool operator==(const Type& x, const Infinite<Type>& y) 
{ 
  return (Infinite<Type>(x) == y); 
}

/// Binary Inequality comparison for two Infinite type objects.
template<class Type>
bool operator!=(const Infinite<Type>& x, const Infinite<Type>& y) 
{ 
  return (!(x == y)); 
}

/// Binary Inequality comparison for one Infinite and one basic type object.
template<class Type>
bool operator!=(const Infinite<Type>& x, const Type& y) 
{ 
  return (x != Infinite<Type>(y)); 
}

/// Binary Inequality comparison for one basic and one Infinite type object.
template<class Type>
bool operator!=(const Type& x, const Infinite<Type>& y) 
{ 
  return (Infinite<Type>(x) != y); 
}

/// Basic output stream operator for printing out the value of the Infinite.
template<class Type>
std::ostream& operator<<(std::ostream& s, const Infinite<Type>& x)
{
  if (x.isInfinite())
    return s << "+INF";
  else if (x.isNegativeInfinite())
    return s << "-INF";
  else
    return s << x.val();
}

/// Basic input stream operator to read in an Infinite value.
template<class Type>
std::istream& operator>>(std::istream& s, Infinite<Type>& x)
{
  Type temp = 0;
  bool neg = false;
  bool inf = false;
  char c = 0;

  // First, see if there is a leading +/- sign.
  s >> c;
  if (c == '-')      // Leading '-' sign?
    neg = true;
  else if (c == '+') // Leading '+' sign?
    neg = false;
  else               // Neither, put the char back on the stream.
    s.putback(c);

  // Now, try to read in a number or "INF".
  s >> c;
  if ((c == 'I') || (c == 'i'))  // Start of "INF"?
    {
      s >> c;
      if ((c == 'N') || (c == 'n'))
        {
          s >> c;
          if ((c == 'F') || (c == 'f'))
            inf = true;
          else 
            s.clear(std::ios::badbit); // Set the state to bad.
        }
      else 
        s.clear(std::ios::badbit); // Set the state to bad.
    }
  else  // Try to read in the rest as a number.
    {
      s.putback(c);
      s >> temp;
    }

  if (s)
    {
      if (inf)
        x = Infinite<Type>(Infinite<Type>::I_INF);
      else
        x = Infinite<Type>(temp);
      if (neg)
        x = (-x);
    }
  return s;
}

template<class Type>
Infinite<Type>::Infinite(Type x)
{
	// First, make sure that this Type has at least as much 
	// precision as an "int" and that this Type is not "unsigned".
	assert((sizeof(Type) >= sizeof(int)) && (Type)(-1) < 0);
	// Set the internal value.
	_val = x;
}

/// Returns the value of the object as a basic type (eg. int, double).
template<class Type>
Type Infinite<Type>::val() const 
{ 
	return _val; 
}

/// Returns true if the value is I_INF (positive infinity).
template<class Type>
bool Infinite<Type>::isInfinite() const 
{ 
	return gmFuzEQ(_val,I_INF); 
}

/// Returns true if the value is I_INF (positive infinity).
template<class Type>
bool Infinite<Type>::isPositiveInfinite() const 
{ 
	return isInfinite(); 
}

/// Returns true if the value is I_NEGINF (negative infinity).
template<class Type>
bool Infinite<Type>::isNegativeInfinite() const 
{ 
	return gmFuzEQ(_val,I_NEGINF); 
}

/// Returns true if the value is either I_INF or I_NEGINF.
template<class Type>
bool Infinite<Type>::hasInfiniteValue() const 
{
	return (isInfinite()||isNegativeInfinite());
}

/// A type conversion operator.
template<class Type>
Infinite<Type>::operator Type() const 
{ 
	return (Type)_val; 
}

/// Returns the absolute value of the object.
template<class Type>
Infinite<Type> Infinite<Type>::fabs() const 
{ 
	return (_val < 0) ? -(*this) : (*this); 
}

/// Unary negation.
template<class Type>
Infinite<Type> Infinite<Type>::operator-() const 
{ 
	return Infinite(-_val); 
}

/// Returns e to the power of the Infinite object
template<class Type>
Infinite<Type> Infinite<Type>::exp() const
{
	if (isInfinite())
		return (*this);
	else if (isNegativeInfinite())
		return Infinite(0);
	else if (gmIsZero(_val))
		return Infinite<Type>(1);
	else
		return Infinite<Type>(::exp(_val));
}

/**
* Returns the square root of the Infinite object.
* @bug Should raise an error flag when called on -INF, but currently
*      returns 0.
*/
template<class Type>
Infinite<Type> Infinite<Type>::sqrt() const
{
	if (isInfinite())
		return (*this);
	else if (isNegativeInfinite())  // Should be an error!
		return Infinite<Type>(0);
	else 
		return Infinite<Type>(::sqrt(_val));
}

/// Power function for Infinite types, returns object^x.
template<class Type>
Infinite<Type> Infinite<Type>::pow(const Infinite<Type>& x) const 
{ 
	Infinite<Type> y = *this;
	return (y ^= x);
}

/// Power function for basic types, returns object^x.
template<class Type>
Infinite<Type> Infinite<Type>::pow(const Type& x) const 
{
	Infinite<Type> y = *this; 
	return (y ^= Infinite<Type>(x)); 
}

/// Minimum of two basic types.
template<class Type>
Infinite<Type> Infinite<Type>::min(const Type& x, const Type& y)
{
	return ((x < y) ? x : y);
}

/// Maximum of two basic types.
template<class Type>
Infinite<Type> Infinite<Type>::max(const Type& x, const Type& y)
{
	return ((x > y) ? x : y);
}

/// Addition/assignment operator for Infinite types.
template<class Type>
Infinite<Type>& Infinite<Type>::operator+=(const Infinite<Type>& x)
{
	if ((isInfinite() && x.isNegativeInfinite()) ||
		(isNegativeInfinite() && x.isInfinite()))
		_val = 0;
	if (isInfinite() || x.isInfinite())
		_val = I_INF;
	else if (isNegativeInfinite() || x.isNegativeInfinite())
		_val = I_NEGINF;
	else 
		_val += x;
	return (*this);
}

/// Addition/assignment operator for basic types.
template<class Type>
Infinite<Type>& Infinite<Type>::operator+=(const Type& x) 
{ 
	return ((*this) += Infinite<Type>(x)); 
}

/// Subtraction/assignment operator for Infinite types.
template<class Type>
Infinite<Type>& Infinite<Type>::operator-=(const Infinite<Type>& x) 
{ 
	return ((*this) += (-x)); 
}

/// Subtraction/assignment operator for basic types.
template<class Type>
Infinite<Type>& Infinite<Type>::operator-=(const Type& x) 
{ 
	return ((*this) -= Infinite<Type>(x)); 
}

/// Multiplication/assignment operator for Infinite types.
template<class Type>
Infinite<Type>& Infinite<Type>::operator*=(const Infinite<Type>& x)
{
	if (gmIsZero(_val) || gmIsZero(x.val()))
		_val = 0;
	else if (hasInfiniteValue() && x.hasInfiniteValue())
    {
		if ((isInfinite() && x.isInfinite()) ||
			(isNegativeInfinite() && x.isNegativeInfinite()))
			_val = I_INF;
		else
			_val = I_NEGINF;
    }
	else if (isInfinite())
		_val = (x.val() > 0) ? I_INF : I_NEGINF;
	else if (isNegativeInfinite())
		_val = (x.val() > 0) ? I_NEGINF : I_INF;
	else if (x.isInfinite())
		_val = (_val > 0) ? I_INF : I_NEGINF;
	else if (x.isNegativeInfinite())
		_val = (_val > 0) ? I_NEGINF : I_INF;
	else
		_val *= x;
	
	return (*this);
}

/// Multiplication/assignment operator for basic types.
template<class Type>
Infinite<Type>& Infinite<Type>::operator*=(const Type& x) 
{ 
	return ((*this) *= Infinite<Type>(x)); 
}

/// Division/assignment operator for Infinite types.
template<class Type>
Infinite<Type>& Infinite<Type>::operator/=(const Infinite<Type>& x)
{
	if (x.hasInfiniteValue())
		_val = 0;
	else if (gmIsZero(x.val()))
		_val = (_val >= 0) ? I_INF : I_NEGINF;
	else if (isInfinite())
		_val = (x.val() >= 0) ? I_INF : I_NEGINF;
	else if (isNegativeInfinite())
		_val = (x.val() >= 0) ? I_NEGINF : I_INF;
	else
		_val /= x;
	return (*this);
}

/// Division/assignment operator for basic types.
template<class Type>
Infinite<Type>& Infinite<Type>::operator/=(const Type& x) 
{ 
	return ((*this) /= Infinite<Type>(x)); 
}

/// Power/assignment operator for Infinite types.
template<class Type>
Infinite<Type>& Infinite<Type>::operator^=(const Infinite<Type>& x)
{
	if (gmIsZero(x.val()))
		_val = 1;
	else if (gmIsZero(_val) || x.isNegativeInfinite())
		_val = 0;
	else if (x.isInfinite())
		_val = I_INF;
	else if (isInfinite())
		_val = (x.val() < 0) ? 0 : I_INF;
	else if (isNegativeInfinite())
    {
		if (x.val() < 0)
			_val = 0;
		else
			_val = (((int)(x.val()) % 2) == 0) ? I_INF : I_NEGINF;
    }
	else
		_val = ::pow(_val,x.val());
	
	return (*this);
}

/// Power/assignment operator for basic types.
template<class Type>
Infinite<Type>& Infinite<Type>::operator^=(const Type& x) 
{ 
	return ((*this) ^= Infinite<Type>(x)); 
}




#endif

