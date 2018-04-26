/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   pistring.h
 * Created: Wed May 17 23:08:04 2000
 * Authors: Gilbert Baumann, Markus Noga, Tim Weyrich
 *     $Id: pistring.h,v 1.9 2003/12/30 20:24:09 weyrich Exp $
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
 * $Log: pistring.h,v $
 * Revision 1.9  2003/12/30 20:24:09  weyrich
 * last commit before continuing at MERL
 *
 * Revision 1.8  2003/12/09 13:31:40  weyrich
 * Win32: stricmp/strnicmp
 *
 * Revision 1.7  2003/11/19 10:37:30  weyrich
 * added Pointshop3D headers
 *
 * Revision 1.6  2003/11/14 17:15:09  weyrich
 * another Win32 fix
 *
 * Revision 1.5  2003/06/23 14:14:40  weyrich
 * more Win32 port...
 *
 * Revision 1.4  2002/05/23 14:02:33  weyrich
 * First official libsfl release.
 *
 * Revision 1.3  2002/05/14 13:27:00  weyrich
 * Win32 up again: fix
 *
 * Revision 1.2  2002/05/07 09:03:12  weyrich
 * _ip; [v]snprintf
 *
 * Revision 1.1.1.2  2001/11/18 20:14:13  weyrich
 * first check-in
 *
 * Revision 1.3  2000/08/04 12:33:20  tim
 * some SGI-port stuff...
 *
 * Revision 1.2  2000/06/29 00:15:08  mlnoga
 * misc win32 porting changes
 *
 * Revision 1.1  2000/05/17 20:36:19  tim
 * added pistring.[ch]
 *
 */

#ifndef __PISTRING__
#define __PISTRING__

#if defined( _WIN32 ) && !defined( __pisocket_h__ )
#  include <pisocket.h>  /* to prevent winsock/winsock2 botch */
#endif

#include <string.h>
#include "pidefs.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __STRICT_ANSI__
#  ifdef _WIN32
  _CRTIMP char * __cdecl strdup(const char *);
  _CRTIMP char * __cdecl strndup(const char *, size_t size);
#  else
  extern PI_API char *strdup(const char *s) __PI_THROW;
  extern PI_API char *strndup(const char *s, size_t size) __PI_THROW;
#  endif
#endif

#ifdef _WIN32
#  define strcasecmp  _stricmp
#  define strncasecmp _strnicmp
#  define snprintf    _snprintf
#  define vsnprintf   _vsnprintf
#elif defined( __STRICT_ANSI__ ) && defined( __GNUC__ )

#  if 1
#    include <stdarg.h>
extern int snprintf (char * __s, size_t __maxlen,
                     const char *__format, ...)
                    __PI_THROW __attribute__ ((__format__ (__printf__, 3, 4)));

extern int vsnprintf (char *__restrict __s, size_t __maxlen,
                      const char *__format, va_list __arg)
                     __PI_THROW __attribute__ ((__format__ (__printf__, 3, 0)));

#  else
extern int snprintf (char *__restrict __s, size_t __maxlen,
                     __const char *__restrict __format, ...)
                    __THROW __attribute__ ((__format__ (__printf__, 3, 4)));

extern int vsnprintf (char *__restrict __s, size_t __maxlen,
                      __const char *__restrict __format, _G_va_list __arg)
                     __THROW __attribute__ ((__format__ (__printf__, 3, 0)));
#  endif
#endif

#if !defined( _WIN32 ) && defined( __STRICT_ANSI__ )

extern PI_API int strcasecmp (const char *__s1, const char *__s2) __PI_THROW;

/* Compare no more than N chars of S1 and S2, ignoring case.  */
extern PI_API int strncasecmp (const char *__s1, const char *__s2,
                               size_t __n) __PI_THROW;
#endif

#ifdef __cplusplus
}
#endif

#endif
