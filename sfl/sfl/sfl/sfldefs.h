// ### CHECK-ME-BEFORE-RELEASE ###
//
// Title:   sfldefs.h
// Created: Sat Apr 28 19:47:12 2001
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
// $Log: sfldefs.h,v $
// Revision 1.13  2004/04/27 15:34:06  weyrich
// experimental: fixing propmask/layout inconsistencies caused by the Pointshop3D writer
//
// Revision 1.12  2003/11/19 13:04:02  weyrich
// added Pointshop3D headers
//
// Revision 1.11  2003/11/18 09:44:45  weyrich
// fix: shininess output; new: clamping of double->uintXX on output
//
// Revision 1.10  2003/11/17 15:05:52  weyrich
// new API allowing for more than 32 surfel properties
//
// Revision 1.9  2003/10/17 11:59:34  weyrich
// added property ids for wuermlin again(!)
//
// Revision 1.8  2003/10/13 12:26:50  wickem
// added flag definitions
//
// Revision 1.6  2003/10/09 14:58:42  wickem
// added isovalue + cell fields
//
// Revision 1.5  2003/03/20 14:52:52  weyrich
// smaller fixes of implicit coordinate system transformations
//
// Revision 1.4  2002/12/18 17:05:45  weyrich
// added feature relations support (experimental)
//
// Revision 1.3  2002/11/18 13:09:12  weyrich
// added markers and tangential vectors
//
// Revision 1.2  2002/05/17 13:58:16  weyrich
// PointShop3D integration nearly finished
//
// Revision 1.1.1.1  2002/05/08 08:58:47  weyrich
// surfel file format
//
// Revision 1.1.1.1  2002/03/26 17:32:13  weyrich
// first check-in
//

#ifndef __SFLDEFS__
#define __SFLDEFS__

/*
 *  Win32 Stuff
 */

#ifdef _WIN32
/* The following ifdef block is the standard way of creating macros which make exporting
 * from a DLL simpler. All files within this DLL are compiled with the SFL_EXPORTS
 * symbol defined on the command line. this symbol should not be defined on any project
 * that uses this DLL. This way any other project whose source files include this file see
 * SFL_API functions as being imported from a DLL, wheras this DLL sees symbols
 * defined with this macro as being exported.
 */
#  ifdef SFL_EXPORTS
#    define SFL_API __declspec(dllexport)
#    include "../stdafx.h"
#  else
#    define SFL_API __declspec(dllimport)
#  endif
#else
#  define SFL_API
#endif

/*
 *  Supported Surfel Properties
 */

namespace sfl
{
  typedef enum
    {
      POSITION_BIT = 0,
      NORMAL_BIT = 1,

      RADIUS_BIT = 4,
      WEIGHT_BIT = 5,
      SAMPLING_DENSITY_BIT = 6,
      TEXTURE_UV_BIT = 7,

      AMBIENT_COLOR_BIT = 8,
      AMBIENT_COEFF_BIT = 9,
      DIFFUSE_COLOR_BIT = 10,
      DIFFUSE_COEFF_BIT = 11,
      SPECULAR_COLOR_BIT = 12,
      SPECULAR_COEFF_BIT = 13,
      SHININESS_BIT = 14,

      TANGENT_AXES_BIT = 24,         /*!< \brief If set, radius and normal are always set, but not stored. */

      FLAGS_BIT = 27,                   /* for possible flags see below */

      ISOVALUE_BIT = 28,                /* (1 float) addded by wickem */
      ORIGIN_CELL_BIT = 29,             /* (1 uint32) added by wickem */

      CAMERA_DEPTH_BIT = 30,            /* (1 float) added by weyrich for wuermlin */
      CAMERA_ID_BIT = 31,               /* (1 uint8) added by weyrich for wuermlin */
      CAMERA_XY_BIT = 32,               /* (2 uint16) added by weyrich for wuermlin */

      NORMAL_INFO_RAY_ONE_BIT = 33,     /* (3 floats) added by weyrich for wuermlin */
      NORMAL_INFO_RAY_TWO_BIT = 34,     /* (3 floats) added by weyrich for wuermlin */
      NORMAL_INFO_THETA_PHI_BIT = 35,   /* (2 floats) added by weyrich for wuermlin */

      MAX_PROPERTY_BITS  /* Highest property index representable in a KillerSurfel */
    }
  PropertyBitEnum;
};

/*
 *  Supported Surfel Properties (Old API; supported for downwards compatibility)
 */

#define  SFLPROP_POSITION             (1 << sfl::POSITION_BIT)
#define  SFLPROP_NORMAL               (1 << sfl::NORMAL_BIT)

#define  SFLPROP_NORMAL_COMPR_MASK    0x0000000c /*!< \brief These two bits encode the normal compression. */
#define  SFLPROP_NORMAL_COMPR_SHIFT            2 /*!< \brief Right-shift, that is needed to right-align the compression mask. */

#define  SFLPROP_RADIUS               (1 << sfl::RADIUS_BIT)
#define  SFLPROP_WEIGHT               (1 << sfl::WEIGHT_BIT)
#define  SFLPROP_SAMPLING_DENSITY     (1 << sfl::SAMPLING_DENSITY_BIT)
#define  SFLPROP_TEXTURE_UV           (1 << sfl::TEXTURE_UV_BIT)

#define  SFLPROP_AMBIENT_COLOR        (1 << sfl::AMBIENT_COLOR_BIT)
#define  SFLPROP_AMBIENT_COEFF        (1 << sfl::AMBIENT_COEFF_BIT)
#define  SFLPROP_DIFFUSE_COLOR        (1 << sfl::DIFFUSE_COLOR_BIT)
#define  SFLPROP_DIFFUSE_COEFF        (1 << sfl::DIFFUSE_COEFF_BIT)
#define  SFLPROP_SPECULAR_COLOR       (1 << sfl::SPECULAR_COLOR_BIT)
#define  SFLPROP_SPECULAR_COEFF       (1 << sfl::SPECULAR_COEFF_BIT)
#define  SFLPROP_SHININESS            (1 << sfl::SHININESS_BIT)

#define  SFLPROP_COLOR_MODEL_L        0x00010000
#define  SFLPROP_COLOR_MODEL_LA       0x00090000
#define  SFLPROP_COLOR_MODEL_RGB      0x00070000
#define  SFLPROP_COLOR_MODEL_RGBA     0x000f0000
#define  SFLPROP_COLOR_MODEL_YUV      0x00100000
#define  SFLPROP_COLOR_MODEL_YUVA     0x00180000

#define  SFLPROP_COLOR_MODEL_MASK     0x001f0000 /*!< \brief masks all color model information */

#define  SFLPROP_TANGENT_AXES         (1 << sfl::TANGENT_AXES_BIT)

#define  SFLPROP_FLAGS                (1 << sfl::FLAGS_BIT)
#define  SFLPROP_USER_FLAGS           SFLPROP_FLAGS

#define  SFLPROP_ISOVALUE             (1 << sfl::ISOVALUE_BIT)
#define  SFLPROP_ORIGIN_CELL          (1 << sfl::ORIGIN_CELL_BIT)

#define  SFLPROP_CAMERA_DEPTH         (1 << sfl::CAMERA_DEPTH_BIT)
#define  SFLPROP_CAMERA_ID            (1 << sfl::CAMERA_ID_BIT)


#define  SFL_MAX_PROPERTY_FLAGS  128  /* internal buffer size (### still needed?) */

/*
 * [User] Flags
 */

//
// PointShop Flags
//

#define  SFLPROP_FLAG_CLEAR            0x00000000   // no flags set

// selections. these do not render the surfel differently
#define  SFLPROP_FLAG_SELECTED1        0x00000001 // selection 1
#define  SFLPROP_FLAG_SELECTED2        0x00000002 // selection 2
#define  SFLPROP_FLAG_SELECTED3        0x00000004 // selection 3

#define  SFLPROP_FLAG_EMPHASISE        0x00000008 // render the surfel differently, eg when selected

#define  SFLPROP_FLAG_COVERED          0x00000010 // used internally (only for temporary computations!)
#define  SFLPROP_FLAG_INSIDE           0x00000020 // used for CSG calculations, true if point is inside other model
#define  SFLPROP_FLAG_FEATURE          0x00000040 // determines if surfel is feature surfel
#define  SFLPROP_FLAG_CLIP_ORIENTATION 0x00000080 // set if union csg (positive clipping) see renderer
#define SFLPROP_FLAG_CLIP_BOUNDARY_A   0x00000100 // surfel close to boundary A that needs to be clipped
#define SFLPROP_FLAG_CLIP_BOUNDARY_B   0x00000200 // surfel close to boundary B that needs to be clipped
#define SFLPROP_FLAG_INVISIBLE         0x00000400 // set the flag if the surfel should not be rendered

//
// Martin Wicke's flags
//

#define SFLPROP_FLAG_EDGESURFEL        0x00000800 // surfel is part of a patch edge (Martin Wicke)
#define SFLPROP_FLAG_EDGEOUT           0x00001000 // surfel center is outside of actual surface (Martin Wicke)

#define SFLPROP_FLAG_ALL_FLAGS         0xffffffff // all flags set



/*
 *  File Format
 */

#include "gff.h"

#define  SFL_MAGIC    "SFLF"
#define  SFL_VERSION  0x0000

#define  SFL_SURFEL_FORMAT_VERSION  0x00000001UL

/* Tags for the main header
 */
#define  SFLTAG_TITLE               0x0100 /* optional */
#define  SFLTAG_AUTHOR              0x0101 /* optional */
#define  SFLTAG_APPLICATION         0x0102 /* optional */
#define  SFLTAG_DESCRIPTION         0x0103 /* optional */
#define  SFLTAG_VERSION             0x0104 /* optional */
#define  SFLTAG_DATE                0x0105 /* optional */

#define  SFLTAG_RIGHT_VEC           0x0110 /* default: 1 0 0 */
#define  SFLTAG_UP_VEC              0x0111 /* default: 0 1 0 */
#define  SFLTAG_BACK_VEC            0x0112 /* default: 0 0 1 */

#define  SFLTAG_DEFAULT_UNITS       0x0113 /* default: 1.0 */

#define  SFLTAG_SCENE_CAMERA_TRANSFORM 0x0115 /* default: id */

#define  SFLTAG_BOUNDING_BOX        0x0118 /* optional */
/* The overall bounding box is given as (x1, y1, z1, x2, y2, z2) in
   world coordinates, using DEFAULT_UNITS. */

#define  SFLTAG_SURFEL_SET_LIST_IDX 0x0120 /* Child index of the surfel-set directory. Mandatory. */

/* Tags for surfel set list header
 */
#define  SFLTAG_NUM_SURFEL_SETS     0x0130 /* mandatory */

/* Tags for surfel set header
 *
 *   Besides the tags below, the following tags may be repeated:
 *
 *     SFLTAG_BOUNDING_BOX
 *
 *   The bounding box is given in object coordinates (before the
 *   transform)
 *
 */
#define  SFLTAG_SURFEL_SET_NUM_RESOLUTIONS  0x0140 /* default: 1 */
#define  SFLTAG_SURFEL_SET_RES0_IDX         0x0141 /* default: numResolutions-1 */
#define  SFLTAG_SURFEL_SET_FORMAT_VERSION   0x0150 /* mandatory (at the moment: (uint32)0x0) */
#define  SFLTAG_SURFEL_SET_SURFEL_SIZE      0x0151 /* mandatory (uint32) */
#define  SFLTAG_SURFEL_SET_PROPERTIES       0x0152 /* bit-vector(!); mandatory (uint32s) */
#define  SFLTAG_SURFEL_SET_PROP_INDICES     0x0153 /* mandatory: field indices corresponding to offsets (uint16s) */
#define  SFLTAG_SURFEL_SET_PROP_OFFSETS     0x0154 /* mandatory: field offsets inside the surfel (uint32s) */
#define  SFLTAG_SURFEL_SET_PROP_TYPES       0x0155 /* not used at the moment (uint32s) */
#define  SFLTAG_SURFEL_SET_UNITS            0x0160 /* not written, at the moment (float) */
#define  SFLTAG_SURFEL_SET_TRANSFORM        0x0161 /* default: id (16 floats)*/
#define  SFLTAG_SURFEL_SET_IDENTIFIER       0x0162 /* default: "surfelSet%02d" */
#define  SFLTAG_SURFEL_SET_MARKER_IDXVEC    0x0163 /* optional; length defines number of markers */
#define  SFLTAG_SURFEL_SET_MARKER_UVEC      0x0164 /* optional; if present, length must be number of markers */
#define  SFLTAG_SURFEL_SET_MARKER_VVEC      0x0165 /* optional; if present, length must be number of markers */
#define  SFLTAG_SURFEL_SET_EXP_FEATURE_RELS 0x0166 /* optional; feature relations; length=2*numRels. Format: directed index tuples */

/* Tags for resolution header
 */
#define  SFLTAG_SURFEL_RSL_NUM_SURFELS      0x0170 /* mandatory */
#define  SFLTAG_SURFEL_RSL_DATA             0x0171 /* LDB follows */

/* Surfel Set Tags used for MERL's Surface Reflectance Fields
 */

#define  SFLTAG_SURFEL_SET_SRF_CALIBFLOATS  0x00050000

#endif

// Some Emacs-Hints -- please don't remove:
//
//  Local Variables:
//  mode:C++
//  tab-width:4
//  End:
