/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   piendian.h
 * Created: Sun Jan 30 22:01:04 2000
 * Authors: Tim Weyrich, Marco Nef
 *     $Id: piendian.h,v 1.12 2003/11/19 10:37:30 weyrich Exp $
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
 * $Log: piendian.h,v $
 * Revision 1.12  2003/11/19 10:37:30  weyrich
 * added Pointshop3D headers
 *
 * Revision 1.11  2003/11/14 17:15:09  weyrich
 * another Win32 fix
 *
 * Revision 1.10  2002/11/26 16:16:49  weyrich
 * nuked pi_NOP_ip bug
 *
 * Revision 1.9  2002/05/14 13:06:14  weyrich
 * Win32 up again
 *
 * Revision 1.8  2002/05/07 09:03:11  weyrich
 * _ip; [v]snprintf
 *
 * Revision 1.7  2002/02/28 15:33:32  weyrich
 * messed around with win32 compatibility
 *
 * Revision 1.6  2002/02/05 10:53:19  mnef
 * no C++ comments anymore
 *
 * Revision 1.5  2001/11/27 15:56:47  mnef
 * renamed pi_endian-macros
 *
 * Revision 1.4  2001/11/27 13:40:25  mnef
 * PI_NETWORK_ENDIAN, pi_to_networkendian, pi_form_networkendian
 *
 * Revision 1.3  2001/11/27 13:25:11  mnef
 * pi_to_bigendian, pi_to_littleendian, pi_from_bigendian, pi_from_littleendian
 *
 * Revision 1.2  2001/11/27 12:48:14  weyrich
 * piendian.h: definition of PI_LOCAL_ENDIANESS is now preprocessor safe
 *
 * Revision 1.1  2001/11/20 16:36:50  weyrich
 * new files: piendian.[ch]
 *
 */

#ifndef __PIENDIAN_H__
#define __PIENDIAN_H__

/* FIXME: the pitypes.h-include must become better -- waiting
 * for Markus to sort his pi-headers
 * */
#include "pidefs.h"
#include "pitypes.h"

/*****************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

#define  PI_BIG_ENDIAN      0x04030201UL
#define  PI_LITTLE_ENDIAN   0x01020304UL

#define  PI_NETWORK_ENDIAN  PI_BIG_ENDIAN
/* Reference: www.cs.rpi.edu/courses/sysprog/sockets/byteorder.html */

#if defined(__PI_LITTLE_ENDIAN) || defined(_WIN32)
#  define  PI_LOCAL_ENDIANESS PI_LITTLE_ENDIAN
#  define  PI_HAVE_BIG_ENDIAN      0
#  define  PI_HAVE_LITTLE_ENDIAN   1
#  define  PI_HAVE_NETWORK_ENDIAN  0
#elif defined(__PI_BIG_ENDIAN)
#  define  PI_LOCAL_ENDIANESS PI_BIG_ENDIAN
#  define  PI_HAVE_BIG_ENDIAN      1
#  define  PI_HAVE_LITTLE_ENDIAN   0
#  define  PI_HAVE_NETWORK_ENDIAN  1
#else
#  define  PI_LOCAL_ENDIANESS  (*(uint32 *)"\4\3\2\1")
#  if  PI_LOCAL_ENDIANESS == PI_BIG_ENDIAN
#    define  PI_HAVE_BIG_ENDIAN      1
#    define  PI_HAVE_LITTLE_ENDIAN   0
#    define  PI_HAVE_NETWORK_ENDIAN  1
#  else
#    define  PI_HAVE_BIG_ENDIAN      0
#    define  PI_HAVE_LITTLE_ENDIAN   1
#    define  PI_HAVE_NETWORK_ENDIAN  0
#  endif
#endif

/* IMPORTANT NOTE: pi_swap_xxx swap in-place! */

extern PI_API uint16   *pi_swap_uint16 (uint16 *ui);
extern PI_API uint32   *pi_swap_uint32 (uint32 *ui);
#ifndef NO_INT64
extern PI_API uint64   *pi_swap_uint64 (uint64 *ui);
#endif

#if !defined(_WIN32) && defined(__STRICT_ANSI__)
#  ifdef _WIN32
#    pragma warning(disable : 4505) /* disable stupid, stupid warnings */
#  endif
#  define  DECL_STYLE  static PI_INLINE
#else
#  define  DECL_STYLE  extern PI_API PI_INLINE
#endif

DECL_STYLE int16   *pi_swap_int16(int16 *i)   { return ((int16 *)pi_swap_uint16(((uint16 *)(i)))); }
DECL_STYLE int32   *pi_swap_int32(int32 *i)   { return ((int32 *)pi_swap_uint32(((uint32 *)(i)))); }
#ifndef NO_INT64
DECL_STYLE sint64  *pi_swap_int64(sint64 *i)  { return ((sint64 *)pi_swap_uint64(((uint64 *)(i)))); }
#endif

DECL_STYLE float   *pi_swap_float(float *f)   { return ((float *)pi_swap_uint32(((uint32 *)(f)))); }
#ifndef NO_INT64
DECL_STYLE double  *pi_swap_double(double *d) { return ((double *)pi_swap_uint64(((uint64 *)(d)))); }
#else
DECL_STYLE double  *pi_swap_double(double *d)
{
  uint32  *a=(uint32*)d, h;
  pi_swap_uint32 (&a[0]);
  pi_swap_uint32 (&a[1]);
  h=a[0]; a[0]=a[1]; a[1]=h;
  return (double*)a;
}
#endif

/* And finally, to make code generation easier: */

DECL_STYLE char   *pi_swap_char(char *c)   { return c; }
DECL_STYLE uchar  *pi_swap_uchar(uchar *c) { return c; }
DECL_STYLE int8   *pi_swap_int8(int8 *i)   { return i; }
DECL_STYLE uint8  *pi_swap_uint8(uint8 *i) { return i; }

#if defined(__STRICT_ANSI__)
#  define FAKE_PIENDIAN_USAGE do{(void)pi_swap_int8;(void)pi_swap_int16;(void)pi_swap_int32;    \
                                 (void)pi_swap_float;(void)pi_swap_double;(void)pi_swap_uint8;  \
                                 (void)pi_swap_char;(void)pi_swap_uchar;}while(0)
#else
#  if 0
#    define FAKE_PIENDIAN_USAGE do{(void)pi_swap_int8;(void)pi_swap_int16;(void)pi_swap_int32;(void)pi_swap_int64;      \
                                  (void)pi_swap_float;(void)pi_swap_double;(void)pi_swap_uint8;(void)pi_swap_char;      \
                                  (void)pi_swap_uchar;}while(0)
#  else
#    define FAKE_PIENDIAN_USAGE
#  endif
#endif

/*****************************************************************************/

#if PI_HAVE_LITTLE_ENDIAN

#  define pi_char_to_be(x)    (*pi_swap_char(&(x)))
#  define pi_uchar_to_be(x)   (*pi_swap_uchar(&(x)))
#  define pi_int8_to_be(x)    (*pi_swap_int8(&(x)))
#  define pi_uint8_to_be(x)   (*pi_swap_uint8(&(x)))
#  define pi_int16_to_be(x)   (*pi_swap_int16(&(x)))
#  define pi_uint16_to_be(x)  (*pi_swap_uint16(&(x)))
#  define pi_int32_to_be(x)   (*pi_swap_int32(&(x)))
#  define pi_uint32_to_be(x)  (*pi_swap_uint32(&(x)))
#  define pi_float_to_be(x)   (*pi_swap_float(&(x)))
#  define pi_double_to_be(x)  (*pi_swap_double(&(x)))

#  define pi_char_to_le(x)    (x)
#  define pi_uchar_to_le(x)   (x)
#  define pi_int8_to_le(x)    (x)
#  define pi_uint8_to_le(x)   (x)
#  define pi_int16_to_le(x)   (x)
#  define pi_uint16_to_le(x)  (x)
#  define pi_int32_to_le(x)   (x)
#  define pi_uint32_to_le(x)  (x)
#  define pi_float_to_le(x)   (x)
#  define pi_double_to_le(x)  (x)

#  define pi_char_from_be(x)   (*pi_swap_char(&(x)))
#  define pi_uchar_from_be(x)  (*pi_swap_uchar(&(x)))
#  define pi_int8_from_be(x)   (*pi_swap_int8(&(x)))
#  define pi_uint8_from_be(x)  (*pi_swap_uint8(&(x)))
#  define pi_int16_from_be(x)  (*pi_swap_int16(&(x)))
#  define pi_uint16_from_be(x) (*pi_swap_uint16(&(x)))
#  define pi_int32_from_be(x)  (*pi_swap_int32(&(x)))
#  define pi_uint32_from_be(x) (*pi_swap_uint32(&(x)))
#  define pi_float_from_be(x)  (*pi_swap_float(&(x)))
#  define pi_double_from_be(x) (*pi_swap_double(&(x)))

#  define pi_char_from_le(x)   (x)
#  define pi_uchar_from_le(x)  (x)
#  define pi_int8_from_le(x)   (x)
#  define pi_uint8_from_le(x)  (x)
#  define pi_int16_from_le(x)  (x)
#  define pi_uint16_from_le(x) (x)
#  define pi_int32_from_le(x)  (x)
#  define pi_uint32_from_le(x) (x)
#  define pi_float_from_le(x)  (x)
#  define pi_double_from_le(x) (x)

#else  /* PI_LOCAL_ENDIANESS == PI_BIGENDIAN */

#  define pi_char_to_be(x)    (x)
#  define pi_uchar_to_be(x)   (x)
#  define pi_int8_to_be(x)    (x)
#  define pi_uint8_to_be(x)   (x)
#  define pi_int16_to_be(x)   (x)
#  define pi_uint16_to_be(x)  (x)
#  define pi_int32_to_be(x)   (x)
#  define pi_uint32_to_be(x)  (x)
#  define pi_float_to_be(x)   (x)
#  define pi_double_to_be(x)  (x)

#  define pi_char_to_le(x)    (*pi_swap_char(&(x)))
#  define pi_uchar_to_le(x)   (*pi_swap_uchar(&(x)))
#  define pi_int8_to_le(x)    (*pi_swap_int8(&(x)))
#  define pi_uint8_to_le(x)   (*pi_swap_uint8(&(x)))
#  define pi_int16_to_le(x)   (*pi_swap_int16(&(x)))
#  define pi_uint16_to_le(x)  (*pi_swap_uint16(&(x)))
#  define pi_int32_to_le(x)   (*pi_swap_int32(&(x)))
#  define pi_uint32_to_le(x)  (*pi_swap_uint32(&(x)))
#  define pi_float_to_le(x)   (*pi_swap_float(&(x)))
#  define pi_double_to_le(x)  (*pi_swap_double(&(x)))

#  define pi_char_from_be(x)          (x)
#  define pi_uchar_from_be(x)         (x)
#  define pi_int8_from_be(x)          (x)
#  define pi_uint8_from_be(x)         (x)
#  define pi_int16_from_be(x)         (x)
#  define pi_uint16_from_be(x)        (x)
#  define pi_int32_from_be(x)         (x)
#  define pi_uint32_from_be(x)        (x)
#  define pi_float_from_be(x)         (x)
#  define pi_double_from_be(x)        (x)

#  define pi_char_from_le(x)          (*pi_swap_char(&(x)))
#  define pi_uchar_from_le(x)         (*pi_swap_uchar(&(x)))
#  define pi_int8_from_le(x)          (*pi_swap_int8(&(x)))
#  define pi_uint8_from_le(x)         (*pi_swap_uint8(&(x)))
#  define pi_int16_from_le(x)         (*pi_swap_int16(&(x)))
#  define pi_uint16_from_le(x)        (*pi_swap_uint16(&(x)))
#  define pi_int32_from_le(x)         (*pi_swap_int32(&(x)))
#  define pi_uint32_from_le(x)        (*pi_swap_uint32(&(x)))
#  define pi_float_from_le(x)         (*pi_swap_float(&(x)))
#  define pi_double_from_le(x)        (*pi_swap_double(&(x)))

#endif

#define pi_char_to_net     pi_char_to_be
#define pi_uchar_to_net    pi_uchar_to_be
#define pi_int8_to_net     pi_int8_to_be
#define pi_uint8_to_net    pi_uint8_to_be
#define pi_int16_to_net    pi_int16_to_be
#define pi_uint16_to_net   pi_uint16_to_be
#define pi_int32_to_net    pi_int32_to_be
#define pi_uint32_to_net   pi_uint32_to_be
#define pi_float_to_net    pi_float_to_be
#define pi_double_to_net   pi_double_to_be

#define pi_char_from_net   pi_char_from_be
#define pi_uchar_from_net  pi_uchar_from_be
#define pi_int8_from_net   pi_int8_from_be
#define pi_uint8_from_net  pi_uint8_from_be
#define pi_int16_from_net  pi_int16_from_be
#define pi_uint16_from_net pi_uint16_from_be
#define pi_int32_from_net  pi_int32_from_be
#define pi_uint32_from_net pi_uint32_from_be
#define pi_float_from_net  pi_float_from_be
#define pi_double_from_net pi_double_from_be

  /** in-place operations ******************************************************/
  /*
   * FIXME: The other macros are working in-place, too! Only the
   * result is different. I'll have to scan all my code to come aroud
   * this problem...
   */

#define  pi_NOP_ip(x)  do{}while(0)

#if PI_HAVE_LITTLE_ENDIAN

#  define pi_char_to_be_ip    pi_swap_char
#  define pi_uchar_to_be_ip   pi_swap_uchar
#  define pi_int8_to_be_ip    pi_swap_int8
#  define pi_uint8_to_be_ip   pi_swap_uint8
#  define pi_int16_to_be_ip   pi_swap_int16
#  define pi_uint16_to_be_ip  pi_swap_uint16
#  define pi_int32_to_be_ip   pi_swap_int32
#  define pi_uint32_to_be_ip  pi_swap_uint32
#  define pi_float_to_be_ip   pi_swap_float
#  define pi_double_to_be_ip  pi_swap_double

#  define pi_char_to_le_ip(x)    pi_NOP_ip(x)
#  define pi_uchar_to_le_ip(x)   pi_NOP_ip(x)
#  define pi_int8_to_le_ip(x)    pi_NOP_ip(x)
#  define pi_uint8_to_le_ip(x)   pi_NOP_ip(x)
#  define pi_int16_to_le_ip(x)   pi_NOP_ip(x)
#  define pi_uint16_to_le_ip(x)  pi_NOP_ip(x)
#  define pi_int32_to_le_ip(x)   pi_NOP_ip(x)
#  define pi_uint32_to_le_ip(x)  pi_NOP_ip(x)
#  define pi_float_to_le_ip(x)   pi_NOP_ip(x)
#  define pi_double_to_le_ip(x)  pi_NOP_ip(x)

#  define pi_char_from_be_ip   pi_swap_char
#  define pi_uchar_from_be_ip  pi_swap_uchar
#  define pi_int8_from_be_ip   pi_swap_int8
#  define pi_uint8_from_be_ip  pi_swap_uint8
#  define pi_int16_from_be_ip  pi_swap_int16
#  define pi_uint16_from_be_ip pi_swap_uint16
#  define pi_int32_from_be_ip  pi_swap_int32
#  define pi_uint32_from_be_ip pi_swap_uint32
#  define pi_float_from_be_ip  pi_swap_float
#  define pi_double_from_be_ip pi_swap_double

#  define pi_char_from_le_ip(x)   pi_NOP_ip(x)
#  define pi_uchar_from_le_ip(x)  pi_NOP_ip(x)
#  define pi_int8_from_le_ip(x)   pi_NOP_ip(x)
#  define pi_uint8_from_le_ip(x)  pi_NOP_ip(x)
#  define pi_int16_from_le_ip(x)  pi_NOP_ip(x)
#  define pi_uint16_from_le_ip(x) pi_NOP_ip(x)
#  define pi_int32_from_le_ip(x)  pi_NOP_ip(x)
#  define pi_uint32_from_le_ip(x) pi_NOP_ip(x)
#  define pi_float_from_le_ip(x)  pi_NOP_ip(x)
#  define pi_double_from_le_ip(x) pi_NOP_ip(x)

#else  /* PI_LOCAL_ENDIANESS == PI_BIGENDIAN */

#  define pi_char_to_be_ip(x)    pi_NOP_ip(x)
#  define pi_uchar_to_be_ip(x)   pi_NOP_ip(x)
#  define pi_int8_to_be_ip(x)    pi_NOP_ip(x)
#  define pi_uint8_to_be_ip(x)   pi_NOP_ip(x)
#  define pi_int16_to_be_ip(x)   pi_NOP_ip(x)
#  define pi_uint16_to_be_ip(x)  pi_NOP_ip(x)
#  define pi_int32_to_be_ip(x)   pi_NOP_ip(x)
#  define pi_uint32_to_be_ip(x)  pi_NOP_ip(x)
#  define pi_float_to_be_ip(x)   pi_NOP_ip(x)
#  define pi_double_to_be_ip(x)  pi_NOP_ip(x)

#  define pi_char_to_le_ip    pi_swap_char
#  define pi_uchar_to_le_ip   pi_swap_uchar
#  define pi_int8_to_le_ip    pi_swap_int8
#  define pi_uint8_to_le_ip   pi_swap_uint8
#  define pi_int16_to_le_ip   pi_swap_int16
#  define pi_uint16_to_le_ip  pi_swap_uint16
#  define pi_int32_to_le_ip   pi_swap_int32
#  define pi_uint32_to_le_ip  pi_swap_uint32
#  define pi_float_to_le_ip   pi_swap_float
#  define pi_double_to_le_ip  pi_swap_double

#  define pi_char_from_be_ip(x)          pi_NOP_ip(x)
#  define pi_uchar_from_be_ip(x)         pi_NOP_ip(x)
#  define pi_int8_from_be_ip(x)          pi_NOP_ip(x)
#  define pi_uint8_from_be_ip(x)         pi_NOP_ip(x)
#  define pi_int16_from_be_ip(x)         pi_NOP_ip(x)
#  define pi_uint16_from_be_ip(x)        pi_NOP_ip(x)
#  define pi_int32_from_be_ip(x)         pi_NOP_ip(x)
#  define pi_uint32_from_be_ip(x)        pi_NOP_ip(x)
#  define pi_float_from_be_ip(x)         pi_NOP_ip(x)
#  define pi_double_from_be_ip(x)        pi_NOP_ip(x)

#  define pi_char_from_le_ip          pi_swap_char
#  define pi_uchar_from_le_ip         pi_swap_uchar
#  define pi_int8_from_le_ip          pi_swap_int8
#  define pi_uint8_from_le_ip         pi_swap_uint8
#  define pi_int16_from_le_ip         pi_swap_int16
#  define pi_uint16_from_le_ip        pi_swap_uint16
#  define pi_int32_from_le_ip         pi_swap_int32
#  define pi_uint32_from_le_ip        pi_swap_uint32
#  define pi_float_from_le_ip         pi_swap_float
#  define pi_double_from_le_ip        pi_swap_double

#endif

#define pi_char_to_net_ip      pi_char_to_be_ip
#define pi_uchar_to_net_ip     pi_uchar_to_be_ip
#define pi_int8_to_net_ip      pi_int8_to_be_ip
#define pi_uint8_to_net_ip     pi_uint8_to_be_ip
#define pi_int16_to_net_ip     pi_int16_to_be_ip
#define pi_uint16_to_net_ip    pi_uint16_to_be_ip
#define pi_int32_to_net_ip     pi_int32_to_be_ip
#define pi_uint32_to_net_ip    pi_uint32_to_be_ip
#define pi_float_to_net_ip     pi_float_to_be_ip
#define pi_double_to_net_ip    pi_double_to_be_ip

#define pi_char_from_net_ip    pi_char_from_be_ip
#define pi_uchar_from_net_ip   pi_uchar_from_be_ip
#define pi_int8_from_net_ip    pi_int8_from_be_ip
#define pi_uint8_from_net_ip   pi_uint8_from_be_ip
#define pi_int16_from_net_ip   pi_int16_from_be_ip
#define pi_uint16_from_net_ip  pi_uint16_from_be_ip
#define pi_int32_from_net_ip   pi_int32_from_be_ip
#define pi_uint32_from_net_ip  pi_uint32_from_be_ip
#define pi_float_from_net_ip   pi_float_from_be_ip
#define pi_double_from_net_ip  pi_double_from_be_ip

/*****************************************************************************/

#ifdef __cplusplus
}
#endif

#endif
