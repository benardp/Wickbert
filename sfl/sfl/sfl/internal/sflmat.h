// ### CHECK-ME-BEFORE-RELEASE ###
//
// Title:   sflmat.h
// Created: Thu May 16 18:41:59 2002
// Authors: Tim Weyrich <weyrich@inf.ethz.ch>
//
// Copyright (c) 2001--2003, Computer Graphics Lab, ETH Zurich
//
// This file is part of the Pointshop3D system.
// See http://www.pointshop3d.com/ for more information.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General
// Public License along with this library; if not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330,
// Boston, MA 02111-1307 USA
//
// Contact info@pointshop3d.com if any conditions of this
// licensing are not clear to you.
//
// ---------------------------------------------------------------
//
// $Log: sflmat.h,v $
// Revision 1.3  2003/11/19 13:04:03  weyrich
// added Pointshop3D headers
//
// Revision 1.2  2003/03/20 14:52:52  weyrich
// smaller fixes of implicit coordinate system transformations
//
// Revision 1.1  2002/05/17 13:58:16  weyrich
// PointShop3D integration nearly finished
//

#ifndef __SFLMAT_H__
#define __SFLMAT_H__

namespace sflMat
{
  void  setIdentity16(float mat[16]);
  void  setIdentity16(double mat[16]);
  void  setIdentity32(float mat[32]);
  void  setIdentity32(double mat[32]);

  void  setZero16(double mat[16]);

  bool  invMat16(double dest[16], const double mat[16]);

  void  trpMat16(float mat[16]);
  void  trpMat16(double mat[16]);
  void  trpMat16(float dest[16], const float mat[16]);
  void  trpMat16(double dest[16], const double mat[16]);

  bool  updateInverse(double mat[32]);

  bool  isIdentity(double mat[16]);
  bool  is3x3(double mat[16]);
  bool  is3x4(double mat[16]);
  bool  hasNegDet(double mat[16]);

  void  mulMat16Mat16(float dest[16], const float a[16], const float b[16]);
  void  mulMat16Mat16(double dest[16], const double a[16], const double b[16]);
  void  mulMat32Mat32(double dest[32], const double a[32], const double b[32]);
  void  mulMat32Inv32(double dest[32], const double a[32], const double bInv[32]);

  void  cpyMat16(float dest[16], const float a[16]);
  void  cpyMat16(float dest[16], const double a[16]);
  void  cpyMat16(double dest[16], const double a[16]);
  void  cpyMat32(float dest[16], const float a[16]);
  void  cpyMat32(double dest[16], const double a[16]);
  void  cpyTrp16(float dest[16], const float a[16]);
  void  cpyTrp16(double dest[16], const double a[16]);

  void  dumpMat32(FILE *out, const double mat[32]);
  void  dumpMat16(FILE *out, const double mat[16]);

  void  debugInversion(void);
};

#endif

// Some Emacs-Hints -- please don't remove:
//
//  Local Variables:
//  mode:C++
//  tab-width:4
//  End:
