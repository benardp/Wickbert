/**
 * This file contains the method definitions for a RandomStream.  To use the
 * RandomStream, simply create a new RandomStream object (either with or
 * without a given seed value) and then continually call next() (or one of
 * its relatives) to get a new random number.  The stream is fairly long
 * (2^32-1 numbers in the stream) so it takes quite a while before you wrap
 * around and generate the same numbers again.
 * @file RandomStream.cpp
 * @date 3 August 2001
 * @author Terry Fleury (tfleury@uiuc.edu)
 */

#include <time.h>
#include <algorithm>
#include "RandomStream.h"

// Define the value for class constant modulus
const long RandomStream::modulus = 2147483647;   ///< 2^31 - 1 

/**
 * Default constructor.  This method uses the current machine time as the
 * seed for the random number stream.
 */
RandomStream::RandomStream()
{
  unsigned long currTime = time(NULL);
  if (currTime > modulus)
    lastZ = currTime - modulus;
  else
    lastZ = currTime;
}

/** 
 * Constructor with a seed number.  This method takes in a particular value
 * for the stream's beginning value.  Thus you can create the same stream of
 * numbers using a particular seed value.
 */
RandomStream::RandomStream(long seed)
{
  lastZ = seed;
  if (lastZ == -modulus-1)  // The largest negative integer doesn't negate!
    lastZ = modulus;
  else if (lastZ < 0)       // Don't allow negative seed values
    lastZ = -lastZ;
}

/**
 * Returns a random integer number in the range [1,2^32-2].  Use this method
 * if you want to get a random number as a long integer value.
 * @return A random number in the range [1,2^32-2].
 */
long RandomStream::nextLong()
{
  lastZ = 48271 * (lastZ % 44488) - 3399 * (lastZ / 44488);
  if (lastZ < 0)
    lastZ += modulus;

  return lastZ;
}

/**
 * Returns a random floating point number in the range (0,1).  Notice that
 * neither 0.0 nor 1.0 can be returned by this method.
 * @return A random number in the range (0,1).
 */
double RandomStream::next()
{
  return (nextLong() / (double)modulus);
}

/**
 * Returns a random floating point number in the range [0,max).  Notice that
 * this method can return 0.0 but cannot return any number >= max (assuming
 * that max is a positive value).
 * @param max The upper limit for the returned random number.
 * @return A random number in the range [0,max).
 */
double RandomStream::next(double max)
{
  return ((nextLong()-1) / (double)(modulus-1)) * max;
}

/**
 * Returns a random floating point number in the range [low,high].  Notice
 * that this method can return values in the full range of low to high,
 * including the endpoint values.  If you attempt to pass in values where
 * low > high, the numbers are automatically swapped to enforce low <= high.
 * @param low The minimum value of the random number.
 * @param high The maximum value of the random number.
 * @return A random number in the range [low, high] inclusive.
 */
double RandomStream::next(double low, double high)
{
  if (low > high)
    std::swap(low,high);
  return ((nextLong()-1) / (double)(modulus-2)) * (high - low) + low;
}

