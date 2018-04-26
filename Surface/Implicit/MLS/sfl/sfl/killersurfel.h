// ### CHECK-ME-BEFORE-RELEASE ###
//
// Title:   killersurfel.h
// Created: Sun Nov 16 15:18:24 2003
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
// $Log: killersurfel.h,v $
// Revision 1.5  2004/04/27 15:34:06  weyrich
// experimental: fixing propmask/layout inconsistencies caused by the Pointshop3D writer
//
// Revision 1.4  2003/11/19 13:04:02  weyrich
// added Pointshop3D headers
//
// Revision 1.3  2003/11/19 09:50:49  weyrich
// fix: g++ optimization bug
//
// Revision 1.2  2003/11/18 09:44:45  weyrich
// fix: shininess output; new: clamping of double->uintXX on output
//
// Revision 1.1  2003/11/17 15:05:52  weyrich
// new API allowing for more than 32 surfel properties
//

#ifndef __KILLERSURFEL_H__
#define __KILLERSURFEL_H__

#include "propmasktype.h"

namespace sfl
{
#include "normalenc.h"

  class SFL_API KillerSurfel
  {
    double  readTransform[32];

    double  positionTransform[16];
    bool    positionTransformIsId;
    bool    positionTransformIs3x3;
    bool    positionTransformIs3x4;
    bool    positionTransformHasNegDet;

    double  normalTransform[16]; // Only preserving normalization, if normalTransformIs3x3 is true!! (or id, of course)
    bool    normalTransformIsId;
    bool    normalTransformIs3x3;
    bool    normalTransformIs3x4;
    bool    normalTransformHasNegDet;

  public:
    uint32  surfelPadding; //!< Constantly 1 for now; will be read from tag SURFEL_SET_SURFEL_PADDING.

  private:
    PropMaskType  propertyMask;
    PropMaskType  propertyHints;

    uint32  propertyOffset[MAX_PROPERTY_BITS];
    uint32  propertyPadding[MAX_PROPERTY_BITS];
    uint32  propertySize[MAX_PROPERTY_BITS];

    int16   layoutWriteUInt8ToUInt8Lookup[3*8*MAX_PROPERTY_BITS][2]; // FIXME: Allocate array dynamically!
    int16   layoutWriteDoubleToFloatLookup[3*MAX_PROPERTY_BITS][2]; // FIXME: Allocate array dynamically!
    int16   layoutWriteDoubleToUInt16Lookup[MAX_PROPERTY_BITS][2]; // FIXME: Allocate array dynamically!
    int16   layoutWriteDoubleToUInt8Lookup[MAX_PROPERTY_BITS][2]; // FIXME: Allocate array dynamically!

    typedef void  (KillerSurfel::*ConvFuncType)(void);
    ConvFuncType  conversions[MAX_PROPERTY_BITS];
    int   numConversions;

    uint32  colorSize;  //!< Disk size of color information, according to the current color model.
    uint32  diskSize;           //!< Disk size of one surfel on disk, including surfel padding.

    void  updatePropertyLayout(int diskSizeOverride=-1, int numProps=-1, const uint16 *propIndcs=NULL, const uint32 *propOffs=NULL);
    void  updateReadConversions(uint32 surfelSetFormatVersion);
    void  updateWriteConversions(void);

    void  convNormalToIco8(void);
    void  convNormalToCub16(void);
    void  convNormalToCub32(void);

    void  convIco8ToNormal(void);
    void  convCub16ToNormal(void);
    void  convCub32ToNormal(void);
    void  convRadiusToSamplingDensity(void);
    void  convRadiusToSamplingDensityWithoutNormal(void);
    void  convSamplingDensityToRadius(void);
    void  convColorModelRGBToL(void);
    void  convColorModelLToRGB(void);
    void  convColorModelRGBToYUV(void);
    void  convColorModelYUVToRGB(void);
    void  convColorModelRGBToYIQ(void);
    void  convColorModelYIQToRGB(void);
    void  convColorModelYUVToYIQ(void);
    void  convColorModelYIQToYUV(void);
    void  convCopyColorAmbientToAny(void);
    void  convCopyColorDiffuseToAny(void);
    void  convCopyColorSpecularToAny(void);

    void  convRadiusNormalToTangentAxes(void);
    void  convTangentAxesToRadius(void);
    void  convTangentAxesToNormal(void);

    void  convPositionTransform3x3(void);
    void  convPositionTransform3x4(void);
    void  convPositionTransform(void);

    void  convNormalTransform3x3(void);
    void  convNormalTransform3x4(void);
    void  convNormalTransform(void);
    void  convTangentTransform3x3(void);
    void  convTangentTransform3x4(void);
    void  convTangentTransform(void);

    void  convTangentTransform3x3NegDet(void);
    void  convTangentTransform3x4NegDet(void);
    void  convTangentTransformNegDet(void);

    uint32  spillGeneric(uint8 *dest) const;
    uint32  slurpGeneric(uint8 *src) const;

  public:
    KillerSurfel();

    double         position[3];
    double         normal[3];
    union {
      uint8   ico8;
      uint16  cub16;
      uint32  cub32;
    } normalCmpr; //!< Compressed normal.
    double         radius;
    double         weight;
    double         samplingDensity[3];
    double         textureUV[2];
    unsigned char  colors[4][4]; //!< RGBA, LxxA, or YUVA for Amb./Diff./Spec./Any.
    unsigned char  coeffs[4];   //!< Coefficients for Amb./Diff./Spec./Any.
    double         shininess;

    double         tangentAxes[2][3];

    unsigned       flags;

    double         isovalue;
    unsigned       origin_cell;

    double         cameraDepth; // currently written as float
    int            cameraId;    // currently written as uint8
    double         cameraXY[2]; // currently written as uint16s (integer coordinates!)

    double         normalInfoRays[2][3]; // currently written as floats
    double         normalInfoThetaPhi[2]; // currently written as floats

    KillerSurfel  &setPropertyMask(const PropMaskType &mask);
    KillerSurfel  &setPropertyMask(const PropMaskType &mask, int diskSizeOverride, int numProps, const uint16 *propIndcs, const uint32 *propOffs);
    const PropMaskType  &getPropertyMask(void) const { return propertyMask; }

    KillerSurfel  &setPropertyHints(unsigned long hints, uint32 surfelSetFormatVersion);
    KillerSurfel  &setPropertyHints(const PropMaskType &hints, uint32 surfelSetFormatVersion);
    const PropMaskType  &getPropertyHints(void) const { return propertyHints; }

    KillerSurfel  &setReadTransform(double readTransform[32]);

    uint32  getDiskSize(void) const { return diskSize; }

    uint32  getPropertyOffset(PropertyBitEnum idx) const { return propertyOffset[(int)idx]; }
    uint32  getPropertyPadding(PropertyBitEnum idx) const { return propertyPadding[(int)idx]; }
    uint32  getPropertySize(PropertyBitEnum idx) const { return propertySize[(int)idx]; }

    PI_INLINE normalenc_compression_enum  getNormalEncoding(void) const
    { return (normalenc_compression_enum)((propertyMask.getUInt32(0) & SFLPROP_NORMAL_COMPR_MASK) >> SFLPROP_NORMAL_COMPR_SHIFT); }

    uint32  (KillerSurfel::*spill)(uint8 *dest) const;
    uint32  (KillerSurfel::*slurp)(uint8 *src) const;

    PI_INLINE void  applyConversions(void)
    {
      int  i;
      for (i=0; i<numConversions; i++)
        (this->*conversions[i])();
    }
  };

}; // namespace sfl


#endif

// Some Emacs-Hints -- please don't remove:
//
//  Local Variables:
//  mode:C++
//  tab-width:4
//  End:
