/**
 * This file contains the graphics math utility functions.
 * @file  gmUtils.h
 * @date  15 June 1994
 * @author Ferdi Scheepers
 * @author Stephen F May
 */

#ifndef GMUTILS_H
#define GMUTILS_H

#include <iostream>
#include <math.h>
#include <assert.h>
#include "bool.h"
#include "gmConst.h"

// ONE-ARGUMENT UTILITY FUNCTIONS

/**
 * Returns the absolute value of a floating point number.
 * Usage: f1 = gmAbs(f2);
 * @param f The input number.
 * @return The absolute value of the input number.
 */
inline double gmAbs(double f)
{
  return (f >= 0) ? f : -f;
}

/**
 * Returns the least integer greater than or equal to a floating point
 * number.  
 * Usage: f1 = gmCeil(f2);
 * @param f The input number.
 * @return The "ceiling" integer of the input number. 
 */
inline double gmCeil(double f)
{
  return (f == int(f)) ? f : (f > 0) ? double(int(f) + 1) : double(int(f));
}

/**
 * Returns the cube of a floating point number.
 * Usage: f1 = gmCube(f2);
 * @param f The input number.
 * @return f * f * f
 */
inline double gmCube(double f)
{
  return f * f * f;
}

/**
 * Converts the input angle (as given in radians) to an angle in degrees.
 * Usage: f1 = gmDegrees(f2);
 * @param f The input angle in radians.
 * @return An angle in degrees.
 */
inline double gmDegrees(double f)
{
  return f * gmRADTODEG;
}

/**
 * Returns the greatest integer less than or equal to a floating point
 * number.
 * Usage: f1 = gmFloor(f2);
 * @param f The input number.
 * @return The "floor" integer of the input number.
 */
inline double gmFloor(double f)
{
  return (f == int(f)) ? f : (f > 0) ? double(int(f)) : double(int(f) - 1);
}

/**
 * Returns the inverse of f assuming f does not equal zero.
 * Usage: f1 = gmInv(f2);
 * @param f The input number.
 * @return The inverse of the number.  Note that there is no error checking
 *         to enforce that the number > 0.
 */
inline double gmInv(double f)
{
  return 1.0 / f;
}

/**
 * Returns a boolean as to whether a floating point number is fuzzy equal to
 * zero.
 * Usage: if (gmIsZero(f)) ...
 * @param f The input number.
 * @return True if the input number is fuzzy equal to 0, false otherwise.
 */
inline bool gmIsZero(double f)
{
  return (gmAbs(f) < gmEPSILON);
}

/**
 * Converts the input angle (as given in degrees) to an angle in radians.
 * Usage: f1 = gmRadians(f2);
 * @param f The input angle in degrees.
 * @return An angle in radians.
 */
inline double gmRadians(double f)
{
  return f * gmDEGTORAD;
}

/**
 * Returns a number corresponding to a floating point number rounded to
 * the nearest integer.
 * Usage: f1 = gmRound(f2);
 * @param f The input number.
 * @return f rounded to the nearest integer.  Note that the return value is
 * still of type double, NOT integer.
 */
inline double gmRound(double f)
{
  return (f >= 0) ? double(int(f + 0.5)) : double(- int(0.5 - f));
}

/**
 * Returns the sign of a floating point number.
 * Usage: f1 = gmSign(f2);
 * @param f The input number.
 * @return -1.0 if f < 0, 1.0 otherwise.
 */
inline double gmSign(double f)
{
  return (f < 0) ? -1.0 : 1.0;
}

/**
 * Returns a smooth hermite interpolate of a floating point number.
 * Usage: f1 = gmSmooth(f2);
 * @param f The input number.
 * @return f^2 * (3 - 2f)
 */
inline double gmSmooth(double f)
{
  return (3.0 - 2.0 * f) * f * f;
}

/**
 * Returns the square of a floating point number.
 * Usage: f1 = gmSqr(f2);
 * @param f The input number.
 * @return f * f
 */
inline double gmSqr(double f)
{
  return f * f;
}

/**
 * Returns a number corresponding to a floating point number truncated to
 * an integer.
 * Usage: f1 = gmTrunc(f2);
 * @param f The input number.
 * @return f truncated to an integer.  Note that the return value is
 * still of type double, NOT integer.
 */
inline double gmTrunc(double f)
{
  return double(int(f));
}

/**
 * Returns zero or the sign of a floating point number.
 * Usage: f1 = gmZSign(f2);
 * @param f The input number.
 * @return If the input number equals 0, then return 0.0, otherwise return
 * the sign of the input number (-1.0 or 1.0).
 */
inline double gmZSign(double f)
{
  return (f > 0) ? 1.0 : (f < 0) ? -1.0 : 0.0;
}

// TWO-ARGUMENT UTILITY FUNCTIONS

/**
 * Returns a boolean as to whether two floating point numbers are fuzzy
 * equal.
 * Usage: if (gmFuzEQ(f1,f2)) ...
 * @param f The first input number.
 * @param g The second input number.
 * @return True if the two input numbers are fuzzy equal, false otherwise.
 */
inline bool gmFuzEQ(double f, double g)
{
  return (f <= g) ? (f >= g - gmEPSILON) : (f <= g + gmEPSILON);
}

/**
 * Returns a boolean as to whether one floating point number is fuzzy
 * greater than or equal to another floating point number.
 * Usage: if (gmFuzGEQ(f1,f2)) ...
 * @param f The first input number.
 * @param g The second input number.
 * @return True if the first input number is fuzzy >= the second input
 *         number, false otherwise.
 */
inline bool gmFuzGEQ(double f, double g)
{
  return (f >= g - gmEPSILON);
}

/**
 * Returns a boolean as to whether one floating point number is fuzzy
 * less than or equal to another floating point number.
 * Usage: if (gmFuzLEQ(f1,f2)) ...
 * @param f The first input number.
 * @param g The second input number.
 * @return True if the first input number is fuzzy <= the second input
 *         number, false otherwise.
 */
inline bool gmFuzLEQ(double f, double g)
{
  return (f <= g + gmEPSILON);
}

/**
 * Returns the maximum of two floating point numbers.
 * Usage: f = gmMax(f1,f2);
 * @param f The first input number.
 * @param g The second input number.
 * @return The maximum value of the two input numbers.
 */
inline double gmMax(double f, double g)
{
  return (f > g) ? f : g;
}

/**
 * Returns the minimum of two floating point numbers.
 * Usage: f = gmMin(f1,f2);
 * @param f The first input number.
 * @param g The second input number.
 * @return The minimum value of the two input numbers.
 */
inline double gmMin(double f, double g)
{
  return (f < g) ? f : g;
}

/**
 * Swaps two floating point numbers. 
 * Usage: gmSwap(f1,f2);
 * @param f The first input number.
 * @param g The second input number.
 */
inline void gmSwap(double& f, double& g)
{
  double gmTmp = f; f = g; g = gmTmp;
}

/**
 * Swaps two integer numbers. 
 * Usage: gmSwap(i1,i2);
 * @param i The first input number.
 * @param j The second input number.
 */
inline void gmSwap(int& i, int& j)
{
  int gmTmp = i; i = j; j = gmTmp;
}

// MULTI-ARGUMENT UTILITY FUNCTIONS

/**
 * Clamps one floating number to be in the range of two other floating point
 * numbers.
 * Usage: gmClamp(f,f1,f2);
 * @param f The number to be clamped to the range [l,h].
 * @param l The low value of the range.
 * @param h The high value of the range.
 * @remarks This function does not check to enforce l <= h.
 */
inline void gmClamp(double &f, double l, double h)
{
  if(f < l) f = l;
  if(f > h) f = h;
}

/**
 * Returns a linear interpolation bewteen two floating point numbers.
 * Usage: f3 = gmLerp(f,f1,f2);
 * @param f The number at which to evaluate the linear interpolation.
 * @param l The value at which f = 0.
 * @param h The value at which f = 1.
 * @return The linear interpolation between l and h evaluated at f.
 */
inline double gmLerp(double f, double l, double h)
{
  return l + ((h - l) * f );
}

/**
 * Returns the maximum of three floating point numbers.
 * Usage: f = gmMax(f1,f2,f3);
 * @param f The first input number.
 * @param g The second input number.
 * @param h The third input number.
 * @return The maximum value of the three input numbers.
 */
inline double gmMax(double f, double g, double h)
{
  return (f > g) ? gmMax(f, h) : gmMax(g, h);
}

/**
 * Returns the minimum of three floating point numbers.
 * Usage: f = gmMin(f1,f2,f3);
 * @param f The first input number.
 * @param g The second input number.
 * @param h The third input number.
 * @return The minimum value of the three input numbers.
 */
inline double gmMin(double f, double g, double h)
{
  return (f < g) ? gmMin(f, h) : gmMin(g, h);
}

/**
 * Returns a hermite interpolation bewteen two floating point numbers.
 * Usage: f3 = gmSlide(f,f1,f2);
 * @param f The number at which to evaluate the hermite interpolation.
 * @param l The value at which f = 0.
 * @param h The value at which f = 1.
 * @return The hermite interpolation between l and h evaluated at f.
 */
inline double gmSlide(double f, double l, double h)
{
  return (f < 0) ? l : (f > 1) ? h : gmLerp(gmSmooth(f), l, h);
}

#endif


