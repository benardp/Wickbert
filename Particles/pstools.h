/**
 * Declaration of particle system helper functions
 * @file pstools.h
 * @date June 12, 2001
 * @author Hans Pedersen
 * @remarks Adapted for use with the Advanced Surface Library by Ed Bachta
 */

#ifndef __ASL_PSTOOLS_H__
#define __ASL_PSTOOLS_H__

#include "libgm/gm.h"

#define SET_FLAG 1
#define NOT_SET_FLAG 0
//this seed can be used if pseudo random behavior is needed.
//e.g. placing only some particles on faces of a mesh.
//on the other hand, orientation coming from the mesh is a different attribute 
//but we would like the same random order
//after its use, srand should be reinitialized by the time() function - EE
static const int STANDARD_RANDOM_SEED=1;

int posIntRand(int maxNb);//<gives a random pos int using rand()%maxNb

double fastExp(double); ///< Fast approximation to exp.
double psRandom(void); ///< Random number generator.

void clampInt(int &x, int max); ///< Clamps x to [0,max-1].

gmVector3 randOrtho(gmVector3 x); ///< Generates a radom orthogonal vector.

#endif

