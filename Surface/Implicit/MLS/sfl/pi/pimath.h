/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   pimath.h
 * Created: Wed May 17 23:54:26 2000
 * Authors: Gilbert Baumann, Markus Noga, Tim Weyrich
 *     $Id: pimath.h,v 1.8 2003/11/21 09:45:30 weyrich Exp $
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
 * $Log: pimath.h,v $
 * Revision 1.8  2003/11/21 09:45:30  weyrich
 * nuked g++-3 warnings
 *
 * Revision 1.7  2003/11/19 10:37:30  weyrich
 * added Pointshop3D headers
 *
 * Revision 1.6  2003/11/14 17:15:09  weyrich
 * another Win32 fix
 *
 * Revision 1.5  2003/09/02 09:47:14  wickem
 * removed prototype for isinf
 *
 * Revision 1.4  2003/06/23 14:14:39  weyrich
 * more Win32 port...
 *
 * Revision 1.3  2002/05/23 14:02:33  weyrich
 * First official libsfl release.
 *
 * Revision 1.2  2002/01/15 14:24:54  weyrich
 * some IEEE-float cosmetics...
 *
 * Revision 1.1.1.2  2001/11/18 20:14:12  weyrich
 * first check-in
 *
 * Revision 1.3  2000/06/29 00:15:08  mlnoga
 * misc win32 porting changes
 *
 * Revision 1.2  2000/06/22 01:14:27  mlnoga
 * first version of low-level sound driver; portability adaptions
 *
 * Revision 1.1  2000/05/17 22:11:28  tim
 * added pimath.[ch]
 *
 */

#ifndef __PIMATH__
#define __PIMATH__

#include "pidefs.h"

#if defined( _WIN32 ) && !defined( __pisocket_h__ )
#  include "pisocket.h"  /* to prevent winsock/winsock2 botch */
#endif

#include <math.h>

#ifndef M_E
#  define       M_E             2.7182818284590452354   /* e */
#endif
#ifndef M_LOG2E
#  define       M_LOG2E         1.4426950408889634074   /* log 2e */
#endif
#ifndef M_LOG10E
#  define       M_LOG10E        0.43429448190325182765  /* log 10e */
#endif
#ifndef M_LN2
#  define       M_LN2           0.69314718055994530942  /* log e2 */
#endif
#ifndef M_LN10
#  define       M_LN10          2.30258509299404568402  /* log e10 */
#endif
#ifndef M_PI
#  define       M_PI            3.14159265358979323846  /* pi */
#endif
#ifndef M_PI_2
#  define       M_PI_2          1.57079632679489661923  /* pi/2 */
#endif
#ifndef M_PI_4
#  define       M_PI_4          0.78539816339744830962  /* pi/4 */
#endif
#ifndef M_1_PI
#  define       M_1_PI          0.31830988618379067154  /* 1/pi */
#endif
#ifndef M_2_PI
#  define       M_2_PI          0.63661977236758134308  /* 2/pi */
#endif
#ifndef M_2_SQRTPI
#  define       M_2_SQRTPI      1.12837916709551257390  /* 2/sqrt(pi) */
#endif
#ifndef M_SQRT2
#  define       M_SQRT2         1.41421356237309504880  /* sqrt(2) */
#endif
#ifndef M_SQRT1_2
#  define       M_SQRT1_2       0.70710678118654752440  /* 1/sqrt(2) */
#endif

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32)
# define isnan _isnan
  extern PI_API int isinf(double);
  extern PI_API double rint(double x);

  _CRTIMP int  __cdecl rand(void);
  static __inline double drand48() {
    return (double) rand();
  }

# define  ilogb(x)  ((int)_logb(x))

#elif defined(__STRICT_ANSI__)

#  ifndef rint
#    define  rint   pi_rint
#  endif
#  ifndef isnan
#    define  isnan  pi_isnan
#  endif

  extern PI_API int pi_isnan(double) __PI_THROW;
/*  extern PI_API int isinf(double) __PI_THROW; */

  extern PI_API double pi_rint(double) __PI_THROW;
#endif

#ifdef __cplusplus
}
#endif

#endif
