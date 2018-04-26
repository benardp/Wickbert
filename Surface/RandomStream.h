/**
 * This file contains the class declaration of a RandomStream.  Basically,
 * this is a never-ending stream of random numbers.  The nice things about
 * this particular random number generator are (1) it has a very long period
 * (2^31 - 1 numbers in the stream); (2) the numbers are uniformly random
 * over the interval; and (3) the stream generates numbers which do well for
 * when pairs of random numbers are needed (ie. minimal lattice effect).
 * @file RandomStream.h
 * @date 3 August 2001
 * @author Terry Fleury (tfleury@uiuc.edu)
 */

#ifndef RANDOMSTREAM_H
#define RANDOMSTREAM_H

class RandomStream 
{
  private:
    static const long modulus;   ///< Constant = 2^31-1
    long lastZ;                  ///< Previous random number

  public:
    /// Default constructor.
    RandomStream();
    /// Constructor with a seed number.
    RandomStream(long);
    /// Returns a random integer number in the range [1,2^32-2].
    long nextLong();
    ///  Returns a random floating point number in the range (0,1).
    double next();
    /// Returns a random floating point number in the range [0,max).
    double next(double);
    /// Returns a random floating point number in the range [low,high].
    double next(double,double);
};

#endif

