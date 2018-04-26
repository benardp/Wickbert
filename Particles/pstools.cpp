/**
 * Implementation of particle system tools
 * @file pstools.cpp
 * @date June 12, 2001
 * @author Hans Pedersen
 * @remarks Adapted for use with the Advanced Surface Library by Ed Bachta
 */

#include <math.h>
#include <stdlib.h>

#include "pstools.h"

/**
 * Generates a random number.
 * @returns The random number.
 */
double psRandom(void) {

	return (double)rand()/RAND_MAX;

} // end psRandom

int posIntRand(int maxNb)
{
	return rand()%maxNb;
}

/**
 * Clamps x to [0, max-1].
 * @param x   Value to clamp.
 * @param max Maximum value to clamp to.
 */
void clampInt(int &x, int max) {

  if (x < 0)
    x = 0;
  else if (x >= max)
    x = max - 1;

} // end clampInt

/**
 * Generates a random vector orthogonal to x.
 * @param   x A vector.
 * @returns A vector orthogonal to x.
 */
gmVector3 randOrtho(gmVector3 x) {

  gmVector3 y(0.0,0.0,0.0);

  if (x[0] != 0.0) {

    y[0] = -x[1];
    y[1] = x[0];

  } else if (x[1] != 0.0) {

    y[1] = -x[2];
    y[2] = x[1];

  } else if (x[2] != 0.0) {

    y[0] = -x[2];

  }

  return y;

} // end randOrtho

/**
 * Does a fast evaluation of exp.
 * @param   arg power to raise e to.
 * @returns e^arg.
 */
double fastExp(double arg) {

#define MY_EXP_RES 10000
#define MY_EXP_RANGE 6.
#define MY_EXP_DELTA .0012
#define MY_EXP_TABLE_RES 10001   // Number of entries in exp function table

  static double exp_table[MY_EXP_TABLE_RES];
  static bool first_time = true;

  // First call initializes the table

  if (first_time) {
    for ( int i = 0; i < MY_EXP_TABLE_RES; i++)
      exp_table[i] = exp(MY_EXP_DELTA*(i-MY_EXP_TABLE_RES/2));
    first_time = false;
  }

  // Look up in precomputed table

  if (arg < -MY_EXP_RANGE) {

    return 0.0;
  
  } else {
    
    short tmp = (short)((arg+MY_EXP_RANGE)/MY_EXP_DELTA);

    if (tmp >= MY_EXP_RES)
      return exp_table[MY_EXP_RES-1];
    else if (tmp < 0)
      return exp_table[0];
    else
      return exp_table[tmp];
  }

} // end fastExp



