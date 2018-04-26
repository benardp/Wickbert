/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   pidefs.h
 * Created: Wed Mar  8 18:26:25 2000
 * Authors: Gilbert Baumann, Markus Noga, Tim Weyrich
 *     $Id: pidefs.h,v 1.5 2003/11/19 10:37:30 weyrich Exp $
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
 * $Log: pidefs.h,v $
 * Revision 1.5  2003/11/19 10:37:30  weyrich
 * added Pointshop3D headers
 *
 * Revision 1.4  2002/05/14 13:06:14  weyrich
 * Win32 up again
 *
 * Revision 1.3  2001/12/12 10:30:56  mnef
 * Fixed warnings about c++ comments in ANSI
 *
 * Revision 1.2  2001/11/27 12:48:14  weyrich
 * piendian.h: definition of PI_LOCAL_ENDIANESS is now preprocessor safe
 *
 * Revision 1.1.1.2  2001/11/18 20:14:12  weyrich
 * first check-in
 *
 * Revision 1.5  2000/08/01 09:02:46  tim
 * USEx()-macros to avoid non-usage warnings
 *
 * Revision 1.4  2000/06/27 09:14:02  tim
 * started doxygen-setup
 *
 * Revision 1.3  2000/06/22 01:14:27  mlnoga
 * first version of low-level sound driver; portability adaptions
 *
 * Revision 1.2  2000/05/17 22:11:11  tim
 * moved some stuff to pimath.h
 *
 * Revision 1.1  2000/05/17 14:37:09  mlnoga
 * moved remaining plugin/pi to system/, patched sources accordingly.
 *
 * Revision 1.4  2000/05/16 08:01:40  gilbert
 * C++ style comments are not ANSI
 *
 * Revision 1.3  2000/05/09 13:58:15  gilbert
 * isnan anders.
 *
 * Revision 1.2  2000/03/18 04:09:52  tim
 * started Win32-adaption
 *
 * Revision 1.1  2000/03/18 02:31:44  tim
 * host objects (and wrapper generation) should work the first time; most of the code is pure ANSI again
 *
 */

#ifndef __PIDEFS__
#define __PIDEFS__

#if defined(_WIN32) && !defined(WIN32)
#  define WIN32
#endif

#ifdef _WIN32
#  ifndef __STRICT_ANSI__
#    define __STRICT_ANSI__
#  endif
#endif

#if defined(__STRICT_ANSI__) && !defined(__STDC__)
#  define __STDC__ 1
#elif !defined(__STRICT_ANSI__) && defined(__STDC__)
#  define __STRICT_ANSI__ 1
#endif

#ifndef __PI_THROW
#  if defined __cplusplus && (__GNUC__ >= 3 || __GNUC_MINOR__ >= 8)
#    define __PI_THROW  throw ()
#  else
#    define __PI_THROW
#  endif
#endif

/*! \brief Indicates a function won't return.
    Applicable only to function prototypes.
 */
#if defined(_WIN32)
#  define PI_NORETURN    __declspec(noreturn)
#elif defined (GNUC)
#  define PI_NORETURN    __attribute__ ((__noreturn__))
#else
/* on SGIs (I still don't know, how to determine, whether we are on an
   SGI, so it is the else-branch :-(   ): */
#  define PI_NORETURN
#endif

/*! \brief Directive to process function inline.
 */
#if defined(_WIN32)
# define PI_INLINE      __inline
#elif defined(__GNUC__)
# define PI_INLINE      __inline__
#else
/* on SGIs (I still don't know, how to determine, whether we are on an
   SGI, so it is the else-branch :-(   ): */
# define PI_INLINE      __inline
#endif
#if defined(PROFILE) /* || defined(__STRICT_ANSI__) */
#  undef  PI_INLINE
#  define PI_INLINE  /* no inlining while profiling/in ANSI-C */
#endif

/*! \brief Struct attribute to force packing for GNU compilers.
 *
 *  Use "#pragma pack(n)" for VC++.
 */
#ifdef  __GNUC__
#  define PI_PACKED __attribute__ ((packed))
#else
#  define PI_PACKED
#endif

#define  USE(v1)                       do{(void)(v1);}while(0)
#define  USE1(v1)                      do{(void)(v1);}while(0)
#define  USE2(v1,v2)                   do{(void)(v1);(void)(v2);}while(0)
#define  USE3(v1,v2,v3)                do{(void)(v1);(void)(v2);(void)(v3);}while(0)
#define  USE4(v1,v2,v3,v4)             do{(void)(v1);(void)(v2);(void)(v3);(void)(v4);}while(0)
#define  USE5(v1,v2,v3,v4,v5)          do{USE4(v1,v2,v3,v4);USE1(v5);}while(0)
#define  USE6(v1,v2,v3,v4,v5,v6)       do{USE4(v1,v2,v3,v4);USE2(v5,v6);}while(0)
#define  USE7(v1,v2,v3,v4,v5,v6,v7)    do{USE4(v1,v2,v3,v4);USE3(v5,v6,v7);}while(0)
#define  USE8(v1,v2,v3,v4,v5,v6,v7,v8) do{USE4(v1,v2,v3,v4);USE4(v5,v6,v7,v8);}while(0)

#ifdef WIN32
/* The following ifdef block is the standard way of creating macros which make exporting
 * from a DLL simpler. All files within this DLL are compiled with the PI_EXPORTS
 * symbol defined on the command line. this symbol should not be defined on any project
 * that uses this DLL. This way any other project whose source files include this file see
 * PI_API functions as being imported from a DLL, wheras this DLL sees symbols
 * defined with this macro as being exported.
 */
#  ifdef PI_EXPORTS
#    define PI_API __declspec(dllexport)
#  else
#    define PI_API __declspec(dllimport)
#  endif
#else
#  define PI_API
#endif

#endif
