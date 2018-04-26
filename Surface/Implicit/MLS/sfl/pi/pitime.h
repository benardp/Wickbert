/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   pitime.h
 * Created: Spring 2000
 * Authors: Markus Noga, Tim Weyrich
 *     $Id: pitime.h,v 1.5 2003/11/19 10:37:30 weyrich Exp $
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
 * $Log: pitime.h,v $
 * Revision 1.5  2003/11/19 10:37:30  weyrich
 * added Pointshop3D headers
 *
 *
 */

/*! \file   pitime_h
    \brief  Platform-independent time routines.
    \author Markus L. Noga, Tim Weyrich
*/

#ifndef __pitime_h__
#define __pitime_h__

#include "pidefs.h"
#include "pitypes.h"

#include <time.h>
#ifndef _WIN32
#  include <sys/time.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
  /* For potential previous definitions of timercmp, etc... */
#  include "win32windows.h"

#  if 0
  struct timeval
  {
    uint32  tv_sec;     /* Seconds.  */
    uint32  tv_usec;    /* Microseconds.  */
  };
#  endif

  struct timezone
  {
    int tz_minuteswest;         /* Minutes west of GMT.  */
    int tz_dsttime;             /* Nonzero if DST is ever in effect.  */
  };

#  if defined __USE_GNU || defined __USE_BSD
  typedef struct timezone *__timezone_ptr_t;
#  else
  typedef void *__timezone_ptr_t;
#  endif

/* Get the current time of day and timezone information,
   putting it into *TV and *TZ.  If TZ is NULL, *TZ is not filled.
   Returns 0 on success, -1 on errors.
   NOTE: This form of timezone information is obsolete.
   Use the functions and variables declared in <time.h> instead.  */
  extern PI_API int gettimeofday (struct timeval * /*__restrict*/ __tv,
                                  __timezone_ptr_t /*__restrict*/ __tz) __PI_THROW;
#endif

/* Convenience macros for operations on timevals.
   NOTE: `timercmp' does not work for >= or <=.  */
#ifndef timerisset
#define	timerisset(tvp)		((tvp)->tv_sec || (tvp)->tv_usec)
#endif
#ifndef timerclear
#define	timerclear(tvp)		((tvp)->tv_sec = (tvp)->tv_usec = 0)
#endif
#ifndef timercmp
#define	timercmp(a, b, CMP) 						      \
  (((a)->tv_sec == (b)->tv_sec) ? 					      \
   ((a)->tv_usec CMP (b)->tv_usec) : 					      \
   ((a)->tv_sec CMP (b)->tv_sec))
#endif
#ifndef timeradd
#define	timeradd(a, b, result)						      \
  do {									      \
    (result)->tv_sec = (a)->tv_sec + (b)->tv_sec;			      \
    (result)->tv_usec = (a)->tv_usec + (b)->tv_usec;			      \
    if ((result)->tv_usec >= 1000000)					      \
      {									      \
	++(result)->tv_sec;						      \
	(result)->tv_usec -= 1000000;					      \
      }									      \
  } while (0)
#endif
#ifndef timersub
#define	timersub(a, b, result)						      \
  do {									      \
    (result)->tv_sec = (a)->tv_sec - (b)->tv_sec;			      \
    (result)->tv_usec = (a)->tv_usec - (b)->tv_usec;			      \
    if ((result)->tv_usec < 0) {					      \
      --(result)->tv_sec;						      \
      (result)->tv_usec += 1000000;					      \
    }									      \
  } while (0)
#endif

/*! Return system time in ms.*/
extern PI_API uint32 piSystemTime(void);

#ifdef __cplusplus
}
#endif

#endif /* __pitime_h__ */
