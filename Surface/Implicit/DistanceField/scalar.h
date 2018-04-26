#ifndef SCALAR_H
#define SCALAR_H

#include <float.h>
#include <assert.h>
#include <cmath>


typedef double scalar;

#define SCALAR_MAX     FLT_MAX
#define SCALAR_MIN     FLT_MIN
#define SCALAR_EPSILON FLT_EPSILON
#define FABS(x)        fabs(x)

#define SQRT(x)        sqrt(x)

#define COS(x)         cos(x)
#define SIN(x)         sin(x)
#define TAN(x)         tan(x)

#define ACOS(x)         acos(x)
#define ASIN(x)         asin(x)
#define ATAN(x)         atan(x)

#define SIGN(x)         ((x < 0) ? (scalar)-1.0 : ((x > 0) ? (scalar) 1.0 : (scalar) 0.0))

#define rand_scalar(x) (x * rand() / (scalar) RAND_MAX)
#define rand_int(x)    (int) ((scalar) (x+1) * (rand() / (scalar) (RAND_MAX + 1.0)))

/*
 * Performing floating point to integer conversion is *very* slow on x86 systems.
 * Here FLOOR and CEIL are implemented using lrint on linux platforms and in asm on windows.
 *
 * Code and techniques are from:
 *
 *      "Fast Rounding of Floating Point Numbers in C/C++ on Wintel Platform"
 *       Laurent de Soras 2004
 *       http://ldesoras.free.fr
 *
 */

#ifdef WIN32 
__inline long int lrintf(float x){
  int i; 
  __asm { fld x 
	  fistp i 
	} 
  return (i); 
}

__inline long int lrint(double x){
  int i; 
  __asm { fld x 
	  fistp i 
	} 
  return (i); 
}
#endif

inline int ROUND(float x){ return lrintf(x); }
inline int FLOOR(float x){ return  (lrintf(2*x-0.5F) >> 1); }
inline int CEIL(float x) { return -(lrintf(-0.5F-2*x) >> 1); }


inline int ROUND(double x){ return lrint(x); }
inline int FLOOR(double x){ return  (lrint(2*x-0.5) >> 1); }
inline int CEIL(double x) { return -(lrint(-0.5-2*x) >> 1); }



//#define	FLOAT_TO_INT(in,out)  __asm__ __volatile__ ("fistpl %0" : "=m" (out) : "t" (in) : "st") ;

#ifndef M_PI //for WIN32
#define M_PI 3.14159265358979323846
#endif

#endif
