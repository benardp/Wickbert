/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   normalenc.h
 * Created: Sun Apr 28 14:21:09 2002
 * Authors: Tim Weyrich <weyrich@inf.ethz.ch>
 *
 * Copyright (c) 2001--2003, Computer Graphics Lab, ETH Zurich
 *
 * This file is part of the Pointshop3D system.
 * See http://www.pointshop3d.com/ for more information.
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
 * $Log: normalenc.h,v $
 * Revision 1.3  2003/11/19 13:04:02  weyrich
 * added Pointshop3D headers
 *
 * Revision 1.2  2003/10/09 14:58:42  wickem
 * added isovalue + cell fields
 *
 * Revision 1.1.1.1  2002/05/08 08:58:47  weyrich
 * surfel file format
 *
 */

#ifndef __NORMALENC_H__
#define __NORMALENC_H__

#ifdef __cplusplus
extern "C" {
#endif

#if 1
#  include <pitypes.h>
#else
  typedef long int32;
  typedef unsigned long uint32;
#endif

  typedef enum
    {
      NONE,
      ICOSAHEDRAL_8,
      CUBIC_16,
      CUBIC_32
    }
  normalenc_compression_enum;

  /* Normal compression expects the vector to be normalized. The only
   * exception is the null vector, which can be compressed, as well.
   *
   * CUBIC_16 guarantees an error below 0.032 (0.023 in monte carlo
   * experiment).
   *
   * CUBIC_32 guarantees an error below 0.00013 (0.000085 in monte
   * carlo experiment)
   *
   * These error measures are the length of the difference vector
   * between an encoded and the decoded vector.
   */

#define  NORMALENC_CODE_INVALID  0xffffffffUL

  uint32  normalenc_encode_normal_ico8(double x, double y, double z);
  uint32  normalenc_encode_normal_cub16(double x, double y, double z);
  uint32  normalenc_encode_normal_cub32(double x, double y, double z);

  int     normalenc_decode_normal_ico8_3f(uint32 code, float *x, float *y, float *z);
  int     normalenc_decode_normal_cub16_3f(uint32 code, float *x, float *y, float *z);
  int     normalenc_decode_normal_cub32_3f(uint32 code, float *x, float *y, float *z);

  int     normalenc_decode_normal_ico8_3d(uint32 code, double *x, double *y, double *z);
  int     normalenc_decode_normal_cub16_3d(uint32 code, double *x, double *y, double *z);
  int     normalenc_decode_normal_cub32_3d(uint32 code, double *x, double *y, double *z);

  uint32  normalenc_encode_normal(normalenc_compression_enum mode,
                                  double x, double y, double z);

  int     normalenc_decode_normal_3f(normalenc_compression_enum mode,
                                     uint32 code, float *x, float *y, float *z);
  int     normalenc_decode_normal_3d(normalenc_compression_enum mode,
                                     uint32 code, double *x, double *y, double *z);

  /* TODO: normalenc_encode_scaled_normal(), etc... */

#ifdef __cplusplus
}
#endif

#endif

/* Some Emacs-Hints -- please don't remove:
 *
 *  Local Variables:
 *  tab-width:4
 *  End:
 */
