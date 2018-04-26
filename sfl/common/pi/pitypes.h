/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   pitypes.h
 * Created: Spring 2000
 * Authors: Gilbert Baumann, Markus Noga, Tim Weyrich
 *     $Id: pitypes.h,v 1.6 2003/11/19 10:37:30 weyrich Exp $
 *
 * Copyright (c) 2000--2003
 * Gilbert Baumann <unk6@stud.uni-karlsruhe.de>
 * Markus Noga <markus@noga.de>
 * Tim Weyrich <weyrich@inf.ethz.ch>
 *
 * This file is part of the Pointshop3D system.
 * See http: *www.pointshop3d.com/ for more information.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330,
 * Boston, MA 02111-1307 USA
 *
 * Contact info@pointshop3d.com if any conditions of this
 * licensing are not clear to you.
 *
 * ---------------------------------------------------------------
 *
 * $Log: pitypes.h,v $
 * Revision 1.6  2003/11/19 10:37:30  weyrich
 * added Pointshop3D headers
 *
 *
 */

/*! \file   pitypes_h
    \brief  Platform-independent system routines.
    \author Gilbert Baumann, Markus L. Noga, Tim Weyrich
*/

#ifndef __pitypes_h__
#define __pitypes_h__

#include "pidefs.h"

#include <limits.h>
#include <float.h>

/* Note:
 *
 * It is important, that pitypes.h doesn't include other files that
 * may depend from code-generation, since some code-gererators include
 * this file!
 * */

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __STRICT_ANSI__
#  undef  NO_INT64
#  define NO_INT64
#endif

/* the following definitions may differ on other machines... */

typedef void*          voidp;

#ifndef _TIFF_DATA_TYPEDEFS_
  typedef signed char      int8;
  typedef unsigned char   uint8;
  typedef short           int16;
  typedef unsigned short uint16;

#  ifdef __alpha
    typedef int            int32;
    typedef unsigned int  uint32;
#  else
    typedef long           int32;
    typedef unsigned long uint32;
#  endif
#endif

typedef unsigned char   uchar;
typedef signed char     sint8;
typedef short          sint16;
#ifdef __alpha
  typedef int           sint32;
#else
  typedef long          sint32;
#endif
  
#ifndef NO_INT64
#  ifdef WIN32
     typedef __int64           int64;
     typedef __int64          sint64;
     typedef unsigned __int64 uint64;
#  else
#    ifndef __STRICT_ANSI__
       typedef long long           int64;
       typedef long long          sint64;
       typedef unsigned long long uint64;
#    endif
#  endif
#endif

#ifdef WIN32
  typedef int32 ssize_t;
#endif

/* other basic types used be several world-related modules: */
typedef int32 bool_t;

#ifndef __cplusplus
#  define true 1
#  define false 0
#endif

#if 0
  /* Deprecated, as far as I'm concerned, since float.h is ANSI (?)
   *
   *                              [ Tim, Fri May 10 11:53:36 2002 ]
   */

   /* Radix of exponent representation */
#  ifndef FLT_RADIX
#    define FLT_RADIX 2
#  endif
     /* Number of base-FLT_RADIX digits in the significand of a float */
#  ifndef FLT_MANT_DIG
#    define FLT_MANT_DIG 24
#  endif
     /* Number of decimal digits of precision in a float */
#  ifndef FLT_DIG
#    define FLT_DIG 6
#  endif
     /* Addition rounds to 0: zero, 1: nearest, 2: +inf, 3: -inf, -1: unknown */
#  ifndef FLT_ROUNDS
#    define FLT_ROUNDS 1
#  endif
     /* Difference between 1.0 and the minimum float greater than 1.0 */
#  ifndef FLT_EPSILON
#    define FLT_EPSILON 1.19209290e-07F
#  endif
     /* Minimum int x such that FLT_RADIX**(x-1) is a normalised float */
#  ifndef FLT_MIN_EXP
#    define FLT_MIN_EXP (-125)
#  endif
     /* Minimum normalised float */
#  ifndef FLT_MIN
#    define FLT_MIN 1.17549435e-38F
#  endif
     /* Minimum int x such that 10**x is a normalised float */
#  ifndef FLT_MIN_10_EXP
#    define FLT_MIN_10_EXP (-37)
#  endif
     /* Maximum int x such that FLT_RADIX**(x-1) is a representable float */
#  ifndef FLT_MAX_EXP
#    define FLT_MAX_EXP 128
#  endif
     /* Maximum float */
#  ifndef FLT_MAX
#    define FLT_MAX 3.40282347e+38F
#  endif
     /* Maximum int x such that 10**x is a representable float */
#  ifndef FLT_MAX_10_EXP
#    define FLT_MAX_10_EXP 38
#  endif

     /* Number of base-FLT_RADIX digits in the significand of a double */
#  ifndef DBL_MANT_DIG
#    define DBL_MANT_DIG 53
#  endif
     /* Number of decimal digits of precision in a double */
#  ifndef DBL_DIG
#    define DBL_DIG 15
#  endif
     /* Difference between 1.0 and the minimum double greater than 1.0 */
#  ifndef DBL_EPSILON
#    define DBL_EPSILON 2.2204460492503131e-16
#  endif
     /* Minimum int x such that FLT_RADIX**(x-1) is a normalised double */
#  ifndef DBL_MIN_EXP
#    define DBL_MIN_EXP (-1021)
#  endif
     /* Minimum normalised double */
#  ifndef DBL_MIN
#    define DBL_MIN 2.2250738585072014e-308
#  endif
     /* Minimum int x such that 10**x is a normalised double */
#  ifndef DBL_MIN_10_EXP
#    define DBL_MIN_10_EXP (-307)
#  endif
     /* Maximum int x such that FLT_RADIX**(x-1) is a representable double */
#  ifndef DBL_MAX_EXP
#    define DBL_MAX_EXP 1024
#  endif
     /* Maximum double */
#  ifndef DBL_MAX
#    define DBL_MAX 1.7976931348623157e+308
#  endif
     /* Maximum int x such that 10**x is a representable double */
#  ifndef DBL_MAX_10_EXP
#    define DBL_MAX_10_EXP 308
#  endif

#endif

#ifdef __cplusplus
}
#endif

#endif
