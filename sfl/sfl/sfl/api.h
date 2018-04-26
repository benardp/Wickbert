// ### CHECK-ME-BEFORE-RELEASE ###
//
// Title:   api.h
// Created: Wed Apr 10 18:09:16 2002
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
// $Log: api.h,v $
// Revision 1.18  2004/04/27 15:34:05  weyrich
// experimental: fixing propmask/layout inconsistencies caused by the Pointshop3D writer
//
// Revision 1.17  2003/11/20 15:24:51  weyrich
// fix: argument list termination PropMaskType::setBits()
//
// Revision 1.16  2003/11/20 09:13:12  weyrich
// win32 fix: dump/readSurfelUserFlags linkage
//
// Revision 1.15  2003/11/19 13:57:37  weyrich
// removed misleading normal/densitity 4-tuple write functions
//
// Revision 1.14  2003/11/19 13:04:00  weyrich
// added Pointshop3D headers
//
// Revision 1.13  2003/11/19 09:50:48  weyrich
// fix: g++ optimization bug
//
// Revision 1.12  2003/11/18 09:44:45  weyrich
// fix: shininess output; new: clamping of double->uintXX on output
//
// Revision 1.11  2003/11/17 15:05:52  weyrich
// new API allowing for more than 32 surfel properties
//
// Revision 1.10  2003/10/09 15:40:57  weyrich
// added properties 'camera_id' and 'camera_weight' for Stephan Wuermlin
//
// Revision 1.9  2003/10/09 14:58:41  wickem
// added isovalue + cell fields
//
// Revision 1.8  2003/03/20 14:52:51  weyrich
// smaller fixes of implicit coordinate system transformations
//
// Revision 1.7  2002/12/18 17:05:45  weyrich
// added feature relations support (experimental)
//
// Revision 1.6  2002/12/12 13:02:31  weyrich
// new: default value for user flags: 0
//
// Revision 1.5  2002/11/18 13:09:12  weyrich
// added markers and tangential vectors
//
// Revision 1.4  2002/06/13 14:24:19  weyrich
// implemented YUV conversions
//
// Revision 1.3  2002/05/23 13:59:56  weyrich
// First official release.
//
// Revision 1.2  2002/05/17 13:58:15  weyrich
// PointShop3D integration nearly finished
//
// Revision 1.1.1.1  2002/05/08 08:58:47  weyrich
// surfel file format
//

#ifndef __SFL_API_H__
#define __SFL_API_H__

#include <math.h>

//#define  HAVE_SFL_UINT64

#ifdef HAVE_SFL_UINT64
namespace sfl
{
#  ifdef _WIN32
  typedef unsigned __int64 uint64;
#  else
  typedef long long uint64;
#  endif
};
#endif

// ### TODO
//
//   - More than 32 property bits
//
//   - Optional: Automatische Anwendung Transformations/Units/CoordSys
//     (Bisher nur CoordSys)
//
//   - Internen Fehlerzustand speichern -- Rueckgabe bei close() oder
//     getError(). (Evtl. auch getErrString() vorsehen.)
//
//   - Switch for error level; redirection of error messages (e.g.,
//     qDebug()).
//
//   - FIXME: SEGV if no properties are dumped (empty property mask)
//

#include <pisystem.h>

#include "sfldefs.h"
#include "killersurfel.h"

namespace sfl
{
  typedef enum
    {
      AMBIENT,
      DIFFUSE,
      SPECULAR,
      ANY                       // Doesn't work for coeffs!!
    }
  colorEnum;

  template<class T> uint8  clampToUInt8(T t) { return t < 0 ? 0 : (t > 255 ? 255 : (uint8)t); }
  template<class T> uint8  clampUnsignedToUInt8(T t) { return t > 255 ? 255 : (uint8)t; }

  template<class T> uint16  clampToUInt16(T t) { return t < 0 ? 0 : (t > 65535 ? 65535 : (uint16)t); }
  template<class T> uint16  clampUnsignedToUInt16(T t) { return t > 65535 ? 65535 : (uint16)t; }

  SFL_API void  memFree(void *mem);     //! Use this to free malloc'ed strings returned by the sfl-API.

  class SFL_API Stream
  {
  protected:
    gff_file       *gff;

    KillerSurfel  surfelBuffer;

    size_t        rwBufferSize;
    uint8        *rwBuffer[2];  //!< Two buffers of size rwBufferSize. rwBuffer[1] follows \e directly rwBuffer[0].
    uint8        *rwBufferHead; //!< Buffer write position.
    uint8        *rwBufferTail; //!< Buffer read position.

    uint32  fillRwBuffer(void); //!< Reads at max rwBufferSize bytes from the gff stream.
    uint32  flushRwBuffer0(void); //!< Writes rwBuffer[0] to the gff stream.
    uint32  flushRwBuffer(void); //!< Writes both buffers to the gff stream.

    bool  apiRowMajorMatrices;

    void  apiSwapMat4x4(float mat[16]) const;
    void  apiSwapMat4x4(double mat[16]) const;
    void  apiSwapMat4x4ToTemp4x4(float temp[16], const float mat[16]) const;
    void  apiSwapMat4x4ToTemp4x4(double temp[16], const double mat[16]) const;

    bool  apiGetTransform4x4(const double *src, float *mat, bool internal=false) const;
    bool  apiGetTransform4x4(const double *src, double *mat, bool internal=false) const;
    bool  apiGetInvTransform4x4(const double *src, float *mat, bool internal=false) const;
    bool  apiGetInvTransform4x4(const double *src, double *mat, bool internal=false) const;
    bool  apiGetRotation4x4(const double *src, float *mat, bool internal=false) const;
    bool  apiGetRotation4x4(const double *src, double *mat, bool internal=false) const;
    bool  apiGetScalationRotation4x4(const double *src, float *mat, bool internal=false) const;
    bool  apiGetScalationRotation4x4(const double *src, double *mat, bool internal=false) const;
    bool  apiGetTranslation4x4(const double *src, float *mat, bool internal=false) const;
    bool  apiGetTranslation4x4(const double *src, double *mat, bool internal=false) const;
    bool  apiGetScalationTranslation4x4(const double *src, float *mat, bool internal=false) const;
    bool  apiGetScalationTranslation4x4(const double *src, double *mat, bool internal=false) const;

    Stream();
    ~Stream();

  public:
    void  setMatrixApiRowMajor(void) { apiRowMajorMatrices = true; }
    void  setMatrixApiColumnMajor(void) { apiRowMajorMatrices = false; }
  };

  class SFL_API InStream : public sfl::Stream
  {
  private:
    gff_directory  *header;
    double          defaultUnits;
    double          fileCoordSystem[32];
    double          appCoordSystem[32];

    double          sceneCameraTransform[32];

    bool            readWorldFlag;

    gff_directory  *surfelSetList;

    gff_directory  *surfelSet;
    int             surfelSetIdx;
    uint32          surfelSetFormatVersion;
    uint32          surfelSetNumResolutions;
    PropMaskType    surfelSetProperties;
    double          surfelSetUnits;
    double          surfelSetTransform[32];
    int             surfelSetNumMarkers;
    int             surfelSetExpNumFeatureRels;
    int             surfelSetSRFNumCalibrations;

    gff_directory  *resolution;
    int             resolutionIdx;
    uint32          resolutionNumSurfels;

    uint32  hdVersion;

    static bool    queryField(gff_directory *dir, int tagId) { return !gff_query_field(dir, tagId, NULL, NULL); }
    static int     queryFieldNumElements(gff_directory *dir, int tagId)
    { int n; if (!gff_query_field(dir, tagId, NULL, &n)) return n; else return -1; }

    static char   *getFieldString(gff_directory *dir, int tagId, const char *deflt=NULL);
    static bool    getFieldString(gff_directory *dir, int tagId, char *dest, size_t buflen, const char *deflt=NULL);
    static bool    getFieldUInt8(gff_directory *dir, int tagId, uint8 &val, bool haveDefault=false, uint8 deflt=0);
    static uint8   getFieldUInt8(gff_directory *dir, int tagId, bool haveDefault=false, uint8 deflt=0);
    static bool    getFieldUInt32(gff_directory *dir, int tagId, uint32 &val, bool haveDefault=false, uint32 deflt=0);
    static uint32  getFieldUInt32(gff_directory *dir, int tagId, bool haveDefault=false, uint32 deflt=0);
    static bool    getFieldDouble(gff_directory *dir, int tagId, double &val, bool haveDefault=false, double deflt=0.0);
    static double  getFieldDouble(gff_directory *dir, int tagId, bool haveDefault=false, double deflt=0.0);

    static uint32  getFieldUInt32ArrIdx(gff_directory *dir, int tagId, int idx, uint32 deflt);
    static double  getFieldDoubleArrIdx(gff_directory *dir, int tagId, int idx, double deflt);

    InStream();
    ~InStream();

  public:
    static InStream  *open(const char *fName, const char *mode="");
    static int        close(InStream *in);

    // Query methods return true, if the corresponding field is defined.
    //

    // Header Field Queries
    bool   queryTitle(void) const       { return queryField(header, SFLTAG_TITLE); }
    bool   queryAuthor(void) const      { return queryField(header, SFLTAG_AUTHOR); }
    bool   queryApplication(void) const { return queryField(header, SFLTAG_APPLICATION); }
    bool   queryDescription(void) const { return queryField(header, SFLTAG_DESCRIPTION); }
    bool   queryVersion(void) const     { return queryField(header, SFLTAG_VERSION); }
    bool   queryDate(void) const        { return queryField(header, SFLTAG_DATE); }

    bool   queryRightVector(void) const { return queryField(header, SFLTAG_RIGHT_VEC); }
    bool   queryUpVector(void) const    { return queryField(header, SFLTAG_UP_VEC); }
    bool   queryBackVector(void) const  { return queryField(header, SFLTAG_BACK_VEC); }

    bool   querySceneCameraTransform(void) const { return queryField(header, SFLTAG_SCENE_CAMERA_TRANSFORM); }

    bool   queryDefaultUnits(void) const{ return queryField(header, SFLTAG_DEFAULT_UNITS); }
    bool   queryOverAllBoundingBox(void) const { return queryField(header, SFLTAG_BOUNDING_BOX); }

    // Get methods returning a boolean return \c true, if the
    // corresponding field is defined. (If a default value exists,
    // they \e always return \c true.)
    //
    // The get methods are not coded very efficiently.
    //

    // Header Field Gets
    char  *getTitle(void) const                         { return getFieldString(header, SFLTAG_TITLE); }
    bool   getTitle(char *dest, size_t buflen) const    { return getFieldString(header, SFLTAG_TITLE, dest, buflen); }
    char  *getAuthor(void) const                        { return getFieldString(header, SFLTAG_AUTHOR); }
    bool   getAuthor(char *dest, size_t buflen) const   { return getFieldString(header, SFLTAG_AUTHOR, dest, buflen); }
    char  *getApplication(void) const                   { return getFieldString(header, SFLTAG_APPLICATION); }
    bool   getApplication(char *dest, size_t buflen) const { return getFieldString(header, SFLTAG_APPLICATION, dest, buflen); }
    char  *getDescription(void) const                   { return getFieldString(header, SFLTAG_DESCRIPTION); }
    bool   getDescription(char *dest, size_t buflen) const {return getFieldString(header, SFLTAG_DESCRIPTION, dest, buflen); }
    bool    getVersion(uint32 &version) const           { return getFieldUInt32(header, SFLTAG_VERSION, version); }
    uint32  getVersion(void) const                      { return getFieldUInt32(header, SFLTAG_VERSION); }
    char   *getDate(void) const                         { return getFieldString(header, SFLTAG_DATE); }
    bool    getDate(char *dest, size_t buflen) const    { return getFieldString(header, SFLTAG_DATE, dest, buflen); }

    bool    getRightVector(double &x, double &y, double &z) const
    {
      x = fileCoordSystem[0];
      y = fileCoordSystem[4];
      z = fileCoordSystem[8];
      return true;
    }
    bool    getUpVector(double &x, double &y, double &z) const
    {
      x = fileCoordSystem[1];
      y = fileCoordSystem[5];
      z = fileCoordSystem[9];
      return true;
    }
    bool    getBackVector(double &x, double &y, double &z) const
    {
      x = fileCoordSystem[2];
      y = fileCoordSystem[6];
      z = fileCoordSystem[10];
      return true;
    }

    bool  getSceneCameraTransform4x4(float *mat) const { return apiGetTransform4x4(sceneCameraTransform, mat); }
    bool  getSceneCameraTransform4x4(double *mat) const { return apiGetTransform4x4(sceneCameraTransform, mat); }
    bool  getSceneCameraInvTransform4x4(float *mat) const { return apiGetInvTransform4x4(sceneCameraTransform, mat); }
    bool  getSceneCameraInvTransform4x4(double *mat) const { return apiGetInvTransform4x4(sceneCameraTransform, mat); }

    bool  getSceneCameraRotation4x4(float *mat) const { return apiGetRotation4x4(sceneCameraTransform, mat); }
    bool  getSceneCameraRotation4x4(double *mat) const { return apiGetRotation4x4(sceneCameraTransform, mat); }
    bool  getSceneCameraScalationRotation4x4(float *mat) const { return apiGetScalationRotation4x4(sceneCameraTransform, mat); }
    bool  getSceneCameraScalationRotation4x4(double *mat) const { return apiGetScalationRotation4x4(sceneCameraTransform, mat); }
    bool  getSceneCameraTranslation4x4(float *mat) const { return apiGetTranslation4x4(sceneCameraTransform, mat); }
    bool  getSceneCameraTranslation4x4(double *mat) const { return apiGetTranslation4x4(sceneCameraTransform, mat); }
    bool  getSceneCameraScalationTranslation4x4(float *mat) const { return apiGetScalationTranslation4x4(sceneCameraTransform, mat); }
    bool  getSceneCameraScalationTranslation4x4(double *mat) const { return apiGetScalationTranslation4x4(sceneCameraTransform, mat); }

    bool    getDefaultUnits(double &units) const        { units = defaultUnits; return true; }
    double  getDefaultUnits(void) const                 { return defaultUnits; }

    bool   getOverAllBoundingBox(double &x1, double &y1, double &z1, double &x2, double &y2, double &z2)
    {
      x1 = getFieldDoubleArrIdx(header, SFLTAG_BOUNDING_BOX, 0, 0.0);
      y1 = getFieldDoubleArrIdx(header, SFLTAG_BOUNDING_BOX, 1, 0.0);
      z1 = getFieldDoubleArrIdx(header, SFLTAG_BOUNDING_BOX, 2, 0.0);
      x2 = getFieldDoubleArrIdx(header, SFLTAG_BOUNDING_BOX, 3, 0.0);
      y2 = getFieldDoubleArrIdx(header, SFLTAG_BOUNDING_BOX, 4, 0.0);
      z2 = getFieldDoubleArrIdx(header, SFLTAG_BOUNDING_BOX, 5, 0.0);
      return queryOverAllBoundingBox();
    }

    bool   getOverAllBoundingBox(float &x1, float &y1, float &z1, float &x2, float &y2, float &z2)
    {
      x1 = (float)getFieldDoubleArrIdx(header, SFLTAG_BOUNDING_BOX, 0, 0.0);
      y1 = (float)getFieldDoubleArrIdx(header, SFLTAG_BOUNDING_BOX, 1, 0.0);
      z1 = (float)getFieldDoubleArrIdx(header, SFLTAG_BOUNDING_BOX, 2, 0.0);
      x2 = (float)getFieldDoubleArrIdx(header, SFLTAG_BOUNDING_BOX, 3, 0.0);
      y2 = (float)getFieldDoubleArrIdx(header, SFLTAG_BOUNDING_BOX, 4, 0.0);
      z2 = (float)getFieldDoubleArrIdx(header, SFLTAG_BOUNDING_BOX, 5, 0.0);
      return queryOverAllBoundingBox();
    }

    // Surfel Set List Field Gets
    bool    getNumSurfelSets(uint32 &numSets) const  { return getFieldUInt32(surfelSetList, SFLTAG_NUM_SURFEL_SETS, numSets, true, 0); }
    uint32  getNumSurfelSets(void) const             { return getFieldUInt32(surfelSetList, SFLTAG_NUM_SURFEL_SETS, true, 0); }

    // Seeking A Single Surfel Set
    bool  seekSurfelSet(int idx);

    // Surfel Set Queries/Gets
    uint32  getSurfelSetFormatVersion(void) const { return surfelSetFormatVersion; }
    uint32  getSurfelSetProperties(int idx) const { return surfelSetProperties.getUInt32(idx); }
    const PropMaskType  &getSurfelSetProperties(void) const { return surfelSetProperties; }

    bool  querySurfelSetRes0Index(void) const   { return queryField(surfelSet, SFLTAG_SURFEL_SET_RES0_IDX); }
    bool  querySurfelSetUnits(void) const       { return queryField(surfelSet, SFLTAG_SURFEL_SET_UNITS); }
    bool  querySurfelSetTransform(void) const   { return queryField(surfelSet, SFLTAG_SURFEL_SET_TRANSFORM); }
    bool  querySurfelSetBoundingBox(void) const { return queryField(surfelSet, SFLTAG_BOUNDING_BOX); }

    bool   getSurfelSetRes0Index(int &idx) const { idx = (int)getFieldUInt32(surfelSet,
                                                                             SFLTAG_SURFEL_SET_RES0_IDX,
                                                                             true, surfelSetNumResolutions-1); return true; }
    int    getSurfelSetRes0Index(void) const     { return (int)getFieldUInt32(surfelSet,
                                                                              SFLTAG_SURFEL_SET_RES0_IDX,
                                                                              true, surfelSetNumResolutions-1); }

    bool    getSurfelSetUnits(double &units) const      { units = surfelSetUnits; return true; }
    double  getSurfelSetUnits(void) const               { return surfelSetUnits; }

    bool  getSurfelSetTransform4x4(float *mat) const { return apiGetTransform4x4(surfelSetTransform, mat); }
    bool  getSurfelSetTransform4x4(double *mat) const { return apiGetTransform4x4(surfelSetTransform, mat); }
    bool  getSurfelSetInvTransform4x4(float *mat) const { return apiGetInvTransform4x4(surfelSetTransform, mat); }
    bool  getSurfelSetInvTransform4x4(double *mat) const { return apiGetInvTransform4x4(surfelSetTransform, mat); }

    bool  getSurfelSetRotation4x4(float *mat) const { return apiGetRotation4x4(surfelSetTransform, mat); }
    bool  getSurfelSetRotation4x4(double *mat) const { return apiGetRotation4x4(surfelSetTransform, mat); }
    bool  getSurfelSetScalationRotation4x4(float *mat) const { return apiGetScalationRotation4x4(surfelSetTransform, mat); }
    bool  getSurfelSetScalationRotation4x4(double *mat) const { return apiGetScalationRotation4x4(surfelSetTransform, mat); }
    bool  getSurfelSetTranslation4x4(float *mat) const { return apiGetTranslation4x4(surfelSetTransform, mat); }
    bool  getSurfelSetTranslation4x4(double *mat) const { return apiGetTranslation4x4(surfelSetTransform, mat); }
    bool  getSurfelSetScalationTranslation4x4(float *mat) const { return apiGetScalationTranslation4x4(surfelSetTransform, mat); }
    bool  getSurfelSetScalationTranslation4x4(double *mat) const { return apiGetScalationTranslation4x4(surfelSetTransform, mat); }
    //
    bool   getSurfelSetBoundingBox(double &x1, double &y1, double &z1, double &x2, double &y2, double &z2)
    {
      x1 = getFieldDoubleArrIdx(surfelSet, SFLTAG_BOUNDING_BOX, 0, 0.0);
      y1 = getFieldDoubleArrIdx(surfelSet, SFLTAG_BOUNDING_BOX, 1, 0.0);
      z1 = getFieldDoubleArrIdx(surfelSet, SFLTAG_BOUNDING_BOX, 2, 0.0);
      x2 = getFieldDoubleArrIdx(surfelSet, SFLTAG_BOUNDING_BOX, 3, 0.0);
      y2 = getFieldDoubleArrIdx(surfelSet, SFLTAG_BOUNDING_BOX, 4, 0.0);
      z2 = getFieldDoubleArrIdx(surfelSet, SFLTAG_BOUNDING_BOX, 5, 0.0);
      return querySurfelSetBoundingBox();
    }

    bool   getSurfelSetBoundingBox(float &x1, float &y1, float &z1, float &x2, float &y2, float &z2)
    {
      x1 = (float)getFieldDoubleArrIdx(surfelSet, SFLTAG_BOUNDING_BOX, 0, 0.0);
      y1 = (float)getFieldDoubleArrIdx(surfelSet, SFLTAG_BOUNDING_BOX, 1, 0.0);
      z1 = (float)getFieldDoubleArrIdx(surfelSet, SFLTAG_BOUNDING_BOX, 2, 0.0);
      x2 = (float)getFieldDoubleArrIdx(surfelSet, SFLTAG_BOUNDING_BOX, 3, 0.0);
      y2 = (float)getFieldDoubleArrIdx(surfelSet, SFLTAG_BOUNDING_BOX, 4, 0.0);
      z2 = (float)getFieldDoubleArrIdx(surfelSet, SFLTAG_BOUNDING_BOX, 5, 0.0);
      return querySurfelSetBoundingBox();
    }

    char  *getSurfelSetIdentifier(void) const;
    bool   getSurfelSetIdentifier(char *dest, size_t buflen) const;

    bool    getSurfelSetNumMarkers(int &numMarkers) const { numMarkers = surfelSetNumMarkers; return true; }
    int     getSurfelSetNumMarkers(void) const            { return surfelSetNumMarkers; }

    bool    getSurfelSetMarkerSurfelIndex(int markerIdx, uint32 &surfelIdx) const
    { surfelIdx = getFieldUInt32ArrIdx(surfelSet, SFLTAG_SURFEL_SET_MARKER_IDXVEC, markerIdx, 0); return true; }
    bool    getSurfelSetMarkerU(int markerIdx, float &u) const
    { u = (float)getFieldDoubleArrIdx(surfelSet, SFLTAG_SURFEL_SET_MARKER_UVEC, markerIdx, 0.0); return true; }
    double  getSurfelSetMarkerU(int markerIdx) const
    { return getFieldDoubleArrIdx(surfelSet, SFLTAG_SURFEL_SET_MARKER_UVEC, markerIdx, 0.0); }
    bool    getSurfelSetMarkerV(int markerIdx, float &v) const
    { v = (float)getFieldDoubleArrIdx(surfelSet, SFLTAG_SURFEL_SET_MARKER_VVEC, markerIdx, 0.0); return true; }
    double  getSurfelSetMarkerV(int markerIdx) const
    { return getFieldDoubleArrIdx(surfelSet, SFLTAG_SURFEL_SET_MARKER_VVEC, markerIdx, 0.0); }

    bool    getSurfelSetExpFeatureRels(int &numRels) const { numRels = surfelSetExpNumFeatureRels; return true; }
    int     getSurfelSetExpFeatureRels(void) const         { return surfelSetExpNumFeatureRels; }

    bool    getSurfelSetExpFeatureRel(int relIdx, uint32 &originIdx, uint32 &neighbourIdx) const
    {
      originIdx = getFieldUInt32ArrIdx(surfelSet, SFLTAG_SURFEL_SET_EXP_FEATURE_RELS, 2*relIdx, 0);
      neighbourIdx = getFieldUInt32ArrIdx(surfelSet, SFLTAG_SURFEL_SET_EXP_FEATURE_RELS, 2*relIdx+1, 0);
      return true;
    }

    bool    getSurfelSetSRFNumCalibrations(int &numCalibs) const { numCalibs = surfelSetSRFNumCalibrations; return true; }
    int     getSurfelSetSRFNumCalibrations(void) const         { return surfelSetSRFNumCalibrations; }

    bool    getSurfelSetSRFCalibration(int idx, double cop[3], double pinv[3][3]) const
    {
      int    i, j;
      for (i=0; i<3; i++)
        cop[i] = getFieldDoubleArrIdx(surfelSet, SFLTAG_SURFEL_SET_SRF_CALIBFLOATS, 3*idx+i, 0);
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          pinv[j][i] = getFieldDoubleArrIdx(surfelSet, SFLTAG_SURFEL_SET_SRF_CALIBFLOATS,
                                            12*surfelSetSRFNumCalibrations + 9*idx+i, 0);
      return true;
    }

    // What coordinate system interpretation does the reader have?
    InStream  &setLocalRightVector(double x, double y, double z);
    InStream  &setLocalUpVector(double x, double y, double z);
    InStream  &setLocalBackVector(double x, double y, double z);

    //! Set coordinate system for reads.
    /*! If called \e before seekSurfelSet(), subsequent calls to
     *  readPosition()/readNormal() yield coordinates in world
     *  coordinates, or object coordinates, respectively, according to
     *  \a readWorld.
     *
     *  For a newly created stream, readWorld is \c false by
     *  default. Changing readAsWorldCoordinates() affects all
     *  subsequent surfel sets until changed again by
     *  setReadAsWorldCoordinates().
     *
     *    \param        readWorld       \c true, if world coordinates shall be read.
     */
    InStream  &setReadAsWorldCoordinates(bool readWorld) { readWorldFlag = readWorld; return *this; }

    //! Announcing Properties, that may be read.
    /*! This allows for some internal optimizations during the read
     *  process. Subsequent read attempts of properties that were \e
     *  not listed in \a propMask will lead to undefined results.
     */
    void  setSurfelSetPropertyHints(uint32 propMask=0xffefffffUL) { surfelBuffer.setPropertyHints(propMask, surfelSetFormatVersion); }

    //! Seeking A Single Resolution.
    /*! Makes resolution \a idx current. Subsequent
     *  beginSurfel()/endSurfel() pairs will refer to resolution \a
     *  idx of the current surfel set. seekResolution() sets the
     *  current surfel read position to the first surfel in the
     *  resolution.
     */
    bool  seekResolution(int idx);

    //! Returns the number of surfels in the current resolution.
    uint32  getResolutionNumSurfels(void) const { return resolutionNumSurfels; }

    PI_INLINE int  beginSurfel(void)
    {
      if (rwBufferTail + surfelBuffer.getDiskSize() > rwBufferHead)
        if (fillRwBuffer() < surfelBuffer.getDiskSize())
          return 1;
      rwBufferTail += (surfelBuffer.*surfelBuffer.slurp)(rwBufferTail);
      surfelBuffer.applyConversions();
      return 0;
    }
    PI_INLINE int  endSurfel(void) { return 0; }

    PI_INLINE void  readSurfelPosition3(double &x, double &y, double &z)
    {
      x = surfelBuffer.position[0];
      y = surfelBuffer.position[1];
      z = surfelBuffer.position[2];
    }

    PI_INLINE void  readSurfelPosition3(float &x, float &y, float &z)
    {
      x = (float)surfelBuffer.position[0];
      y = (float)surfelBuffer.position[1];
      z = (float)surfelBuffer.position[2];
    }

    PI_INLINE void  readSurfelPosition4(double &x, double &y, double &z, double &w)
    {
      x = surfelBuffer.position[0];
      y = surfelBuffer.position[1];
      z = surfelBuffer.position[2];
      w = 1.0;
    }

    PI_INLINE void  readSurfelPosition4(float &x, float &y, float &z, float &w)
    {
      x = (float)surfelBuffer.position[0];
      y = (float)surfelBuffer.position[1];
      z = (float)surfelBuffer.position[2];
      w = 1.0f;
    }

    PI_INLINE void  readSurfelNormal3(double &x, double &y, double &z)
    {
      x = surfelBuffer.normal[0];
      y = surfelBuffer.normal[1];
      z = surfelBuffer.normal[2];
    }

    PI_INLINE void  readSurfelNormal3(float &x, float &y, float &z)
    {
      x = (float)surfelBuffer.normal[0];
      y = (float)surfelBuffer.normal[1];
      z = (float)surfelBuffer.normal[2];
    }

    PI_INLINE void  readSurfelNormal4(double &x, double &y, double &z, double &w)
    {
      x = surfelBuffer.normal[0];
      y = surfelBuffer.normal[1];
      z = surfelBuffer.normal[2];
      w = 0.0;
    }

    PI_INLINE void  readSurfelNormal4(float &x, float &y, float &z, float &w)
    {
      x = (float)surfelBuffer.normal[0];
      y = (float)surfelBuffer.normal[1];
      z = (float)surfelBuffer.normal[2];
      w = 0.0f;
    }

    PI_INLINE double  readSurfelRadius(void) { return surfelBuffer.radius; }
    PI_INLINE void    readSurfelRadius(double &r) { r = surfelBuffer.radius; }
    PI_INLINE void    readSurfelRadius(float &r) { r = (float)surfelBuffer.radius; }

    PI_INLINE double  readSurfelWeight(void) { return surfelBuffer.weight; }
    PI_INLINE void    readSurfelWeight(double &w) { w = surfelBuffer.weight; }
    PI_INLINE void    readSurfelWeight(float &w) { w = (float)surfelBuffer.weight; }

    PI_INLINE void  readSurfelSamplingDensity3(double &x, double &y, double &z)
    {
      x = surfelBuffer.samplingDensity[0];
      y = surfelBuffer.samplingDensity[1];
      z = surfelBuffer.samplingDensity[2];
    }

    PI_INLINE void  readSurfelSamplingDensity3(float &x, float &y, float &z)
    {
      x = (float)surfelBuffer.samplingDensity[0];
      y = (float)surfelBuffer.samplingDensity[1];
      z = (float)surfelBuffer.samplingDensity[2];
    }

    PI_INLINE void  readSurfelSamplingDensity4(double &x, double &y, double &z, double &w)
    {
      x = surfelBuffer.samplingDensity[0];
      y = surfelBuffer.samplingDensity[1];
      z = surfelBuffer.samplingDensity[2];
      w = 1.0;
    }

    PI_INLINE void  readSurfelSamplingDensity4(float &x, float &y, float &z, float &w)
    {
      x = (float)surfelBuffer.samplingDensity[0];
      y = (float)surfelBuffer.samplingDensity[1];
      z = (float)surfelBuffer.samplingDensity[2];
      w = 1.0f;
    }

    PI_INLINE void  readSurfelTextureUV(float &u, float &v)
    {
      u = (float)surfelBuffer.textureUV[0];
      v = (float)surfelBuffer.textureUV[1];
    }

    PI_INLINE double  readSurfelColorRf(colorEnum colIdx)
    {
      return (1.0/255.0) * surfelBuffer.colors[(int)colIdx][0];
    }

    PI_INLINE void  readSurfelColorRf(colorEnum colIdx, double &r)
    {
      r = (1.0/255.0) * surfelBuffer.colors[(int)colIdx][0];
    }

    PI_INLINE void  readSurfelColorRf(colorEnum colIdx, float &r)
    {
      r = (float)((1.0/255.0) * surfelBuffer.colors[(int)colIdx][0]);
    }

    PI_INLINE double  readSurfelColorGf(colorEnum colIdx)
    {
      return (1.0/255.0) * surfelBuffer.colors[(int)colIdx][1];
    }

    PI_INLINE void  readSurfelColorGf(colorEnum colIdx, double &g)
    {
      g = (1.0/255.0) * surfelBuffer.colors[(int)colIdx][1];
    }

    PI_INLINE void  readSurfelColorGf(colorEnum colIdx, float &g)
    {
      g = (float)((1.0/255.0) * surfelBuffer.colors[(int)colIdx][1]);
    }

    PI_INLINE double  readSurfelColorBf(colorEnum colIdx)
    {
      return (1.0/255.0) * surfelBuffer.colors[(int)colIdx][2];
    }

    PI_INLINE void  readSurfelColorBf(colorEnum colIdx, double &b)
    {
      b = (1.0/255.0) * surfelBuffer.colors[(int)colIdx][2];
    }

    PI_INLINE void  readSurfelColorBf(colorEnum colIdx, float &b)
    {
      b = (float)((1.0/255.0) * surfelBuffer.colors[(int)colIdx][2]);
    }

    PI_INLINE double  readSurfelColorAf(colorEnum colIdx)
    {
      return (1.0/255.0) * surfelBuffer.colors[(int)colIdx][3];
    }

    PI_INLINE void  readSurfelColorAf(colorEnum colIdx, double &a)
    {
      a = (1.0/255.0) * surfelBuffer.colors[(int)colIdx][3];
    }

    PI_INLINE void  readSurfelColorAf(colorEnum colIdx, float &a)
    {
      a = (float)((1.0/255.0) * surfelBuffer.colors[(int)colIdx][3]);
    }

    PI_INLINE double  readSurfelColorLf(colorEnum colIdx)
    {
      return (1.0/255.0) * surfelBuffer.colors[(int)colIdx][0];
    }

    PI_INLINE void  readSurfelColorLf(colorEnum colIdx, double &l)
    {
      l = (1.0/255.0) * surfelBuffer.colors[(int)colIdx][0];
    }

    PI_INLINE void  readSurfelColorLf(colorEnum colIdx, float &l)
    {
      l = (float)((1.0/255.0) * surfelBuffer.colors[(int)colIdx][0]);
    }

    PI_INLINE void  readSurfelColorLAf(colorEnum colIdx, double &l, double &a)
    {
      l = (1.0/255.0) * surfelBuffer.colors[(int)colIdx][0];
      a = (1.0/255.0) * surfelBuffer.colors[(int)colIdx][3];
    }

    PI_INLINE void  readSurfelColorLAf(colorEnum colIdx, float &l, float &a)
    {
      l = (float)((1.0/255.0) * surfelBuffer.colors[(int)colIdx][0]);
      a = (float)((1.0/255.0) * surfelBuffer.colors[(int)colIdx][3]);
    }

    PI_INLINE void  readSurfelColorRGBf(colorEnum colIdx, double &r, double &g, double &b)
    {
      r = (1.0/255.0) * surfelBuffer.colors[(int)colIdx][0];
      g = (1.0/255.0) * surfelBuffer.colors[(int)colIdx][1];
      b = (1.0/255.0) * surfelBuffer.colors[(int)colIdx][2];
    }

    PI_INLINE void  readSurfelColorRGBf(colorEnum colIdx, float &r, float &g, float &b)
    {
      r = (float)((1.0/255.0) * surfelBuffer.colors[(int)colIdx][0]);
      g = (float)((1.0/255.0) * surfelBuffer.colors[(int)colIdx][1]);
      b = (float)((1.0/255.0) * surfelBuffer.colors[(int)colIdx][2]);
    }

    PI_INLINE void  readSurfelColorRGBAf(colorEnum colIdx, double &r, double &g, double &b, double &a)
    {
      r = (1.0/255.0) * surfelBuffer.colors[(int)colIdx][0];
      g = (1.0/255.0) * surfelBuffer.colors[(int)colIdx][1];
      b = (1.0/255.0) * surfelBuffer.colors[(int)colIdx][2];
      a = (1.0/255.0) * surfelBuffer.colors[(int)colIdx][3];
    }

    PI_INLINE void  readSurfelColorRGBAf(colorEnum colIdx, float &r, float &g, float &b, float &a)
    {
      r = (float)((1.0/255.0) * surfelBuffer.colors[(int)colIdx][0]);
      g = (float)((1.0/255.0) * surfelBuffer.colors[(int)colIdx][1]);
      b = (float)((1.0/255.0) * surfelBuffer.colors[(int)colIdx][2]);
      a = (float)((1.0/255.0) * surfelBuffer.colors[(int)colIdx][3]);
    }

#   define  readSurfelColorYf     readSurfelColorRf
#   define  readSurfelColorUf     readSurfelColorGf
#   define  readSurfelColorVf     readSurfelColorBf
#   define  readSurfelColorYUVf   readSurfelColorRGBf
#   define  readSurfelColorYUVAf  readSurfelColorRGBAf

    PI_INLINE unsigned char  readSurfelColorRuc(colorEnum colIdx)
    {
      return surfelBuffer.colors[(int)colIdx][0];
    }

    PI_INLINE void  readSurfelColorRuc(colorEnum colIdx, unsigned char &r)
    {
      r = surfelBuffer.colors[(int)colIdx][0];
    }

    PI_INLINE unsigned char  readSurfelColorGuc(colorEnum colIdx)
    {
      return surfelBuffer.colors[(int)colIdx][1];
    }

    PI_INLINE void  readSurfelColorGuc(colorEnum colIdx, unsigned char &g)
    {
      g = surfelBuffer.colors[(int)colIdx][1];
    }

    PI_INLINE unsigned char  readSurfelColorBuc(colorEnum colIdx)
    {
      return surfelBuffer.colors[(int)colIdx][2];
    }

    PI_INLINE void  readSurfelColorBuc(colorEnum colIdx, unsigned char &b)
    {
      b = surfelBuffer.colors[(int)colIdx][2];
    }

    PI_INLINE unsigned char  readSurfelColorAuc(colorEnum colIdx)
    {
      return surfelBuffer.colors[(int)colIdx][3];
    }

    PI_INLINE void  readSurfelColorAuc(colorEnum colIdx, unsigned char &a)
    {
      a = surfelBuffer.colors[(int)colIdx][3];
    }

    PI_INLINE unsigned char  readSurfelColorLuc(colorEnum colIdx)
    {
      return surfelBuffer.colors[(int)colIdx][0];
    }

    PI_INLINE void  readSurfelColorLuc(colorEnum colIdx, unsigned char &l)
    {
      l = surfelBuffer.colors[(int)colIdx][0];
    }

    PI_INLINE void  readSurfelColorLAuc(colorEnum colIdx, unsigned char &l, unsigned char &a)
    {
      l = surfelBuffer.colors[(int)colIdx][0];
      a = surfelBuffer.colors[(int)colIdx][3];
    }

    PI_INLINE void  readSurfelColorRGBuc(colorEnum colIdx, unsigned char &r, unsigned char &g, unsigned char &b)
    {
      r = surfelBuffer.colors[(int)colIdx][0];
      g = surfelBuffer.colors[(int)colIdx][1];
      b = surfelBuffer.colors[(int)colIdx][2];
    }

    PI_INLINE void  readSurfelColorRGBAuc(colorEnum colIdx, unsigned char &r, unsigned char &g, unsigned char &b, unsigned char &a)
    {
      r = surfelBuffer.colors[(int)colIdx][0];
      g = surfelBuffer.colors[(int)colIdx][1];
      b = surfelBuffer.colors[(int)colIdx][2];
      a = surfelBuffer.colors[(int)colIdx][3];
    }

#   define  readSurfelColorYuc     readSurfelColorRuc
#   define  readSurfelColorUuc     readSurfelColorGuc
#   define  readSurfelColorVuc     readSurfelColorBuc
#   define  readSurfelColorYUVuc   readSurfelColorRGBuc
#   define  readSurfelColorYUVAuc  readSurfelColorRGBAuc

    PI_INLINE double  readSurfelColorCoefff(colorEnum colIdx)
    {
      // FIXME: Is that Tisu's semantic?
      return (1.0/255.0) * surfelBuffer.coeffs[(int)colIdx];
    }

    PI_INLINE void  readSurfelColorCoefff(colorEnum colIdx, double &c)
    {
      // FIXME: Is that Tisu's semantic?
      c = (1.0/255.0) * surfelBuffer.coeffs[(int)colIdx];
    }

    PI_INLINE void  readSurfelColorCoefff(colorEnum colIdx, float &c)
    {
      // FIXME: Is that Tisu's semantic?
      c = (float)((1.0/255.0) * surfelBuffer.coeffs[(int)colIdx]);
    }

    PI_INLINE unsigned char  readSurfelColorCoeffuc(colorEnum colIdx)
    {
      return surfelBuffer.coeffs[(int)colIdx];
    }

    PI_INLINE void  readSurfelColorCoeffuc(colorEnum colIdx, unsigned char &c)
    {
      c = surfelBuffer.coeffs[(int)colIdx];
    }

    PI_INLINE void  readSurfelShininess(double &s) { s = surfelBuffer.shininess; }
    PI_INLINE void  readSurfelShininess(float &s) { s = (float)surfelBuffer.shininess; }
    PI_INLINE void  readSurfelShininess(uint8 &s) { s = clampToUInt8(surfelBuffer.shininess + 0.5); }

    PI_INLINE void  readSurfelTangentAxisOne3(double &x, double &y, double &z)
    {
      x = surfelBuffer.tangentAxes[0][0];
      y = surfelBuffer.tangentAxes[0][1];
      z = surfelBuffer.tangentAxes[0][2];
    }

    PI_INLINE void  readSurfelTangentAxisOne3(float &x, float &y, float &z)
    {
      x = (float)surfelBuffer.tangentAxes[0][0];
      y = (float)surfelBuffer.tangentAxes[0][1];
      z = (float)surfelBuffer.tangentAxes[0][2];
    }

    PI_INLINE void  readSurfelTangentAxisTwo3(double &x, double &y, double &z)
    {
      x = surfelBuffer.tangentAxes[1][0];
      y = surfelBuffer.tangentAxes[1][1];
      z = surfelBuffer.tangentAxes[1][2];
    }

    PI_INLINE void  readSurfelTangentAxisTwo3(float &x, float &y, float &z)
    {
      x = (float)surfelBuffer.tangentAxes[1][0];
      y = (float)surfelBuffer.tangentAxes[1][1];
      z = (float)surfelBuffer.tangentAxes[1][2];
    }

    PI_INLINE uint32  readSurfelFlags(void) { 
      return (surfelSetProperties.isSet(FLAGS_BIT)) ? surfelBuffer.flags : 0; 
    }
    PI_INLINE void  readSurfelFlags(unsigned &_flags) { 
      _flags = readSurfelFlags(); 
    }
    PI_INLINE void  readSurfelFlags(uint32 &_flags) { 
      _flags = readSurfelFlags(); 
    }

	// Deprecated (user flags are not 'flags'):
    PI_INLINE uint32  readSurfelUserFlags(void) { return readSurfelFlags(); }
    PI_INLINE void  readSurfelUserFlags(unsigned &_flags) { _flags = readSurfelFlags(); }
    PI_INLINE void  readSurfelUserFlags(uint32 &_flags) { _flags = readSurfelFlags(); }


    PI_INLINE double  readSurfelIsovalue(void) { 
      return surfelBuffer.isovalue; 
    }
    PI_INLINE void    readSurfelIsovalue(double &w) { 
      w = surfelBuffer.isovalue; 
    }
    PI_INLINE void readSurfelIsovalue(float &w) { 
      w = (float)surfelBuffer.isovalue; 
    }

    
    PI_INLINE double readSurfelOriginCell(void) { 
      return surfelBuffer.origin_cell; 
    }
    PI_INLINE void readSurfelOriginCell(uint32 &w) { 
      w = surfelBuffer.origin_cell; 
    }
    
    PI_INLINE double  readSurfelCameraDepth(void) { 
      return surfelBuffer.cameraDepth; 
    }
    PI_INLINE void    readSurfelCameraDepth(double &d) { 
      d = surfelBuffer.cameraDepth; 
    }
    PI_INLINE void readSurfelCameraDepth(float &d) { 
      d = (float)surfelBuffer.cameraDepth; 
    }

    PI_INLINE void readSurfelCameraId(uint8 &id) { id = surfelBuffer.cameraId; }
    PI_INLINE uint8  readSurfelCameraId(void) { return surfelBuffer.cameraId; }
    
    PI_INLINE void  readSurfelCameraXY(double &x, double &y) {
      x = surfelBuffer.cameraXY[0];
      y = surfelBuffer.cameraXY[1];
    }
    PI_INLINE void  readSurfelCameraXY(float &x, float &y) {
      x = (float)surfelBuffer.cameraXY[0];
      y = (float)surfelBuffer.cameraXY[1];
    }
    PI_INLINE void  readSurfelCameraXY(int &x, int &y) {
      x = (int)floor(surfelBuffer.cameraXY[0] + 0.5);
      y = (int)floor(surfelBuffer.cameraXY[1] + 0.5);
    }

    PI_INLINE void  readSurfelNormalInfoRayOne3(double &x, double &y, double &z) {
      x = surfelBuffer.normalInfoRays[0][0];
      y = surfelBuffer.normalInfoRays[0][1];
      z = surfelBuffer.normalInfoRays[0][2];
    }
    PI_INLINE void  readSurfelNormalInfoRayOne3(float &x, float &y, float &z) {
      x = (float)surfelBuffer.normalInfoRays[0][0];
      y = (float)surfelBuffer.normalInfoRays[0][1];
      z = (float)surfelBuffer.normalInfoRays[0][2];
    }

    PI_INLINE void  readSurfelNormalInfoRayTwo3(double &x, double &y, double &z) {
      x = surfelBuffer.normalInfoRays[1][0];
      y = surfelBuffer.normalInfoRays[1][1];
      z = surfelBuffer.normalInfoRays[1][2];
    }
    PI_INLINE void  readSurfelNormalInfoRayTwo3(float &x, float &y, float &z) {
      x = (float)surfelBuffer.normalInfoRays[1][0];
      y = (float)surfelBuffer.normalInfoRays[1][1];
      z = (float)surfelBuffer.normalInfoRays[1][2];
    }
    
    PI_INLINE void  readSurfelNormalInfoThetaPhi(double &phi, double &theta) {
      phi   = surfelBuffer.normalInfoThetaPhi[0];
      theta = surfelBuffer.normalInfoThetaPhi[1];
    }
    PI_INLINE void  readSurfelNormalInfoThetaPhi(float &phi, float &theta) {
      phi   = (float)surfelBuffer.normalInfoThetaPhi[0];
      theta = (float)surfelBuffer.normalInfoThetaPhi[1];
    }
  };

  class SFL_API OutStream : public sfl::Stream
  {
    bool  wroteHeader;
    bool  wroteSurfelSetList;

    bool  inHeader;
    bool  inSurfelSetList;
    bool  inSurfelSet;
    bool  inResolution;

    int  numAnnouncedSurfelSets;
    int  numWrittenSurfelSets;

    int  numAnnouncedResolutions;
    int  numWrittenResolutions;

    int  numAnnouncedSurfels;
    int  numWrittenSurfels;

    int  gffDepth;

    void  assertHeader(void);
    void  assertInSurfelSetList(void);

    OutStream();
    ~OutStream();

  public:
    static OutStream  *open(const char *fName, const char *mode="");
    static int         close(OutStream *out);

    int  beginHeader(void);
    OutStream  &setTitle(const char *title);
    OutStream  &setAuthor(const char *author);
    OutStream  &setApplication(const char *appl);
    OutStream  &setDescription(const char *descr);
    OutStream  &setVersion(uint32 version);
    OutStream  &setDate(const char *date);
    OutStream  &setRightVec3(double x, double y, double z);
    OutStream  &setRightVec3(const float *vec);
    OutStream  &setRightVec3(const double *vec);
    OutStream  &setUpVec3(double x, double y, double z);
    OutStream  &setUpVec3(const float *vec);
    OutStream  &setUpVec3(const double *vec);
    OutStream  &setBackVec3(double x, double y, double z);
    OutStream  &setBackVec3(const float *vec);
    OutStream  &setBackVec3(const double *vec);
    OutStream  &setDefaultUnits(double units);
    OutStream  &setOverAllBoundingBox(double x1, double y1, double z1, double x2, double y2, double z2);
    OutStream  &setSceneCameraTransform4x4(const float *mat);
    OutStream  &setSceneCameraTransform4x4(const double *mat);
    OutStream  &setSceneCameraTransformTwice4x4(const float *matSecond, const float *matFirst);
    OutStream  &setSceneCameraTransformTwice4x4(const double *matSecond, const double *matFirst);
    int  endHeader(void);

    int  beginSurfelSetList(uint32 numSurfelSets);
    // call begin/endSurfelSet() pairs inbetween
    int  endSurfelSetList(void);

    int  beginSurfelSet(unsigned propertyMask, ssize_t numResolutions=1);
    int  beginSurfelSet(const PropMaskType &propertyMask, ssize_t numResolutions=1);
    int  setSurfelSetRes0Index(int idx);
    int  setSurfelSetUnits(double units);
    int  setSurfelSetTransform3x3(const float *mat);
    int  setSurfelSetTransform4x4(const float *mat);
    int  setSurfelSetTransform3x3(const double *mat);
    int  setSurfelSetTransform4x4(const double *mat);
    int  setSurfelSetTransformTwice4x4(const float *matSecond, const float *matFirst);
    int  setSurfelSetTransformTwice4x4(const double *matSecond, const double *matFirst);
    int  setSurfelSetBoundingBox(double x1, double y1, double z1, double x2, double y2, double z2);
    // call begin/endResolution() pairs inbetween
    int  setSurfelSetIdentifier(const char *ident);
    int  setSurfelSetMarkerSurfelIndices(int n, const uint32 *idxvec);
    int  setSurfelSetMarkerSurfelUs(int n, const float *uvec); //!< Never write a different number of u's than indices!!
    int  setSurfelSetMarkerSurfelVs(int n, const float *vvec); //!< Never write a different number of v's than indices!!
    int  setSurfelSetExpFeatureRels(int nRels, const uint32 *fromToIdxVec); //!< \a fromToIdxVec must contain 2*\a n indices!
    int  setSurfelSetSRFCalibrations(int nCalibs, const float *cops, const float *pinvs); //!< \a cops must contain 3*\a nCalibs elements, \a pinvs 9*\a nCalibs elements!
    int  setSurfelSetSRFCalibrations(int nCalibs, const double *cops, const double *pinvs); //!< \a cops must contain 3*\a nCalibs elements, \a pinvs 9*\a nCalibs elements!
    int  endSurfelSet(void);

    int  beginResolution(ssize_t numSurfels);
    // call begin/endSurfel() pairs inbetween
    int  endResolution(void);

    PI_INLINE int  beginSurfel(void) { return 0; }
    PI_INLINE int  endSurfel(void)
    {
      surfelBuffer.applyConversions();
      rwBufferHead += (surfelBuffer.*surfelBuffer.spill)(rwBufferHead);
      if (rwBufferHead >= rwBuffer[1])
        if (flushRwBuffer0() != rwBufferSize)
          return 1;
      numWrittenSurfels++;
      return 0;
    }

    PI_INLINE void  dumpSurfelPosition3(double x, double y, double z)
    {
      surfelBuffer.position[0] = x;
      surfelBuffer.position[1] = y;
      surfelBuffer.position[2] = z;
    }

    PI_INLINE void  dumpSurfelPosition4(double x, double y, double z, double w)
    {
      surfelBuffer.position[0] = x/w;
      surfelBuffer.position[1] = y/w;
      surfelBuffer.position[2] = z/w;
    }

    PI_INLINE void  dumpSurfelNormal3(double x, double y, double z)
    {
      surfelBuffer.normal[0] = x;
      surfelBuffer.normal[1] = y;
      surfelBuffer.normal[2] = z;
    }

    PI_INLINE void  dumpSurfelRadius(double r) { surfelBuffer.radius = r; }

    PI_INLINE void  dumpSurfelWeight(double w) { surfelBuffer.weight = w; }

    PI_INLINE void  dumpSurfelSamplingDensity3(double x, double y, double z)
    {
      surfelBuffer.samplingDensity[0] = x;
      surfelBuffer.samplingDensity[1] = y;
      surfelBuffer.samplingDensity[2] = z;
    }

    PI_INLINE void  dumpSurfelTextureUV(double u, double v)
    {
      surfelBuffer.textureUV[0] = u;
      surfelBuffer.textureUV[1] = v;
    }

    PI_INLINE void  dumpSurfelColorRf(colorEnum colIdx, float r)
    {
      surfelBuffer.colors[(int)colIdx][0] = (unsigned char)(255.0f*r + 0.5);
    }

    PI_INLINE void  dumpSurfelColorGf(colorEnum colIdx, float g)
    {
      surfelBuffer.colors[(int)colIdx][1] = (unsigned char)(255.0f*g + 0.5);
    }

    PI_INLINE void  dumpSurfelColorBf(colorEnum colIdx, float b)
    {
      surfelBuffer.colors[(int)colIdx][2] = (unsigned char)(255.0f*b + 0.5);
    }

    PI_INLINE void  dumpSurfelColorAf(colorEnum colIdx, float a)
    {
      surfelBuffer.colors[(int)colIdx][3] = (unsigned char)(255.0f*a + 0.5);
    }

    PI_INLINE void  dumpSurfelColorLf(colorEnum colIdx, float l)
    {
      surfelBuffer.colors[(int)colIdx][0] = (unsigned char)(255.0f*l + 0.5);
      surfelBuffer.colors[(int)colIdx][1] = surfelBuffer.colors[(int)colIdx][0];
      surfelBuffer.colors[(int)colIdx][2] = surfelBuffer.colors[(int)colIdx][1];
    }

    PI_INLINE void  dumpSurfelColorLAf(colorEnum colIdx, float l, float a)
    {
      surfelBuffer.colors[(int)colIdx][0] = (unsigned char)(255.0f*l + 0.5);
      surfelBuffer.colors[(int)colIdx][1] = surfelBuffer.colors[(int)colIdx][0];
      surfelBuffer.colors[(int)colIdx][2] = surfelBuffer.colors[(int)colIdx][1];
      surfelBuffer.colors[(int)colIdx][3] = (unsigned char)(255.0f*a + 0.5);
    }

    PI_INLINE void  dumpSurfelColorRGBf(colorEnum colIdx, float r, float g, float b)
    {
      surfelBuffer.colors[(int)colIdx][0] = (unsigned char)(255.0f*r + 0.5);
      surfelBuffer.colors[(int)colIdx][1] = (unsigned char)(255.0f*g + 0.5);
      surfelBuffer.colors[(int)colIdx][2] = (unsigned char)(255.0f*b + 0.5);
    }

    PI_INLINE void  dumpSurfelColorRGBAf(colorEnum colIdx, float r, float g, float b, float a)
    {
      surfelBuffer.colors[(int)colIdx][0] = (unsigned char)(255.0f*r + 0.5);
      surfelBuffer.colors[(int)colIdx][1] = (unsigned char)(255.0f*g + 0.5);
      surfelBuffer.colors[(int)colIdx][2] = (unsigned char)(255.0f*b + 0.5);
      surfelBuffer.colors[(int)colIdx][3] = (unsigned char)(255.0f*a + 0.5);
    }

#   define  dumpSurfelColorYf     dumpSurfelColorRf
#   define  dumpSurfelColorUf     dumpSurfelColorGf
#   define  dumpSurfelColorVf     dumpSurfelColorBf
#   define  dumpSurfelColorYUVf   dumpSurfelColorRGBf
#   define  dumpSurfelColorYUVAf  dumpSurfelColorRGBAf

    PI_INLINE void  dumpSurfelColorRuc(colorEnum colIdx, unsigned char r)
    {
      surfelBuffer.colors[(int)colIdx][0] = r;
    }

    PI_INLINE void  dumpSurfelColorGuc(colorEnum colIdx, unsigned char g)
    {
      surfelBuffer.colors[(int)colIdx][1] = g;
    }

    PI_INLINE void  dumpSurfelColorBuc(colorEnum colIdx, unsigned char b)
    {
      surfelBuffer.colors[(int)colIdx][2] = b;
    }

    PI_INLINE void  dumpSurfelColorAuc(colorEnum colIdx, unsigned char a)
    {
      surfelBuffer.colors[(int)colIdx][3] = a;
    }

    PI_INLINE void  dumpSurfelColorLuc(colorEnum colIdx, unsigned char l)
    {
      surfelBuffer.colors[(int)colIdx][0] = l;
      surfelBuffer.colors[(int)colIdx][1] = l;
      surfelBuffer.colors[(int)colIdx][2] = l;
    }

    PI_INLINE void  dumpSurfelColorLAuc(colorEnum colIdx, unsigned char l, unsigned char a)
    {
      surfelBuffer.colors[(int)colIdx][0] = l;
      surfelBuffer.colors[(int)colIdx][1] = l;
      surfelBuffer.colors[(int)colIdx][2] = l;
      surfelBuffer.colors[(int)colIdx][3] = a;
    }

    PI_INLINE void  dumpSurfelColorRGBuc(colorEnum colIdx, unsigned char r, unsigned char g, unsigned char b)
    {
      surfelBuffer.colors[(int)colIdx][0] = r;
      surfelBuffer.colors[(int)colIdx][1] = g;
      surfelBuffer.colors[(int)colIdx][2] = b;
    }

    PI_INLINE void  dumpSurfelColorRGBAuc(colorEnum colIdx, unsigned char r, unsigned char g, unsigned char b, unsigned char a)
    {
      surfelBuffer.colors[(int)colIdx][0] = r;
      surfelBuffer.colors[(int)colIdx][1] = g;
      surfelBuffer.colors[(int)colIdx][2] = b;
      surfelBuffer.colors[(int)colIdx][3] = a;
    }

#   define  dumpSurfelColorYuc     dumpSurfelColorRuc
#   define  dumpSurfelColorUuc     dumpSurfelColorGuc
#   define  dumpSurfelColorVuc     dumpSurfelColorBuc
#   define  dumpSurfelColorYUVuc   dumpSurfelColorRGBuc
#   define  dumpSurfelColorYUVAuc  dumpSurfelColorRGBAuc

    PI_INLINE void  dumpSurfelColorCoefff(colorEnum colIdx, float c)
    {
      // FIXME: Is that Tisu's semantic?
      surfelBuffer.coeffs[(int)colIdx] = (unsigned char)(255.0f*c + 0.5);
    }

    PI_INLINE void  dumpSurfelColorCoeffuc(colorEnum colIdx, unsigned char c)
    {
      surfelBuffer.coeffs[(int)colIdx] = c;
    }

    PI_INLINE void  dumpSurfelShininess(double s) { surfelBuffer.shininess = s; }

    PI_INLINE void  dumpSurfelTangentAxisOne3(double x, double y, double z)
    {
      surfelBuffer.tangentAxes[0][0] = x;
      surfelBuffer.tangentAxes[0][1] = y;
      surfelBuffer.tangentAxes[0][2] = z;
    }

    PI_INLINE void  dumpSurfelTangentAxisTwo3(double x, double y, double z)
    {
      surfelBuffer.tangentAxes[1][0] = x;
      surfelBuffer.tangentAxes[1][1] = y;
      surfelBuffer.tangentAxes[1][2] = z;
    }

    PI_INLINE void  dumpSurfelFlags(unsigned _flags) {
      surfelBuffer.flags = _flags;
    }

    // Deprecated:
    PI_INLINE void  dumpSurfelUserFlags(unsigned _flags) { dumpSurfelFlags(_flags); }

    PI_INLINE void dumpSurfelIsovalue(double v) {
      surfelBuffer.isovalue = v;
    }

    PI_INLINE void dumpSurfelOriginCell(uint32 c) {
      surfelBuffer.origin_cell = c;
    }

    PI_INLINE void  dumpSurfelCameraId(uint8 id) { surfelBuffer.cameraId = id; }
    PI_INLINE void  dumpSurfelCameraDepth(double d) { surfelBuffer.cameraDepth = d; }
    PI_INLINE void  dumpSurfelCameraXY(double x, double y) {
      surfelBuffer.cameraXY[0] = x;
      surfelBuffer.cameraXY[1] = y;
    }

    PI_INLINE void  dumpSurfelNormalInfoRayOne3(double x, double y, double z) {
      surfelBuffer.normalInfoRays[0][0] = x;
      surfelBuffer.normalInfoRays[0][1] = y;
      surfelBuffer.normalInfoRays[0][2] = z;
    }
    PI_INLINE void  dumpSurfelNormalInfoRayTwo3(double x, double y, double z) {
      surfelBuffer.normalInfoRays[1][0] = x;
      surfelBuffer.normalInfoRays[1][1] = y;
      surfelBuffer.normalInfoRays[1][2] = z;
    }

    PI_INLINE void  dumpSurfelNormalInfoThetaPhi(double phi, double theta) {
      surfelBuffer.normalInfoThetaPhi[0] = phi;
      surfelBuffer.normalInfoThetaPhi[1] = theta;
    }
  };
};

#endif

// Some Emacs-Hints -- please don't remove:
//
//  Local Variables:
//  mode:C++
//  tab-width:4
//  End:
