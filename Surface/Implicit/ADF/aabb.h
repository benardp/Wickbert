#ifndef SPACEAABB_INCLUDED // -*- C++ -*-
#define SPACEAABB_INCLUDED

#ifdef _MSC_VER
#pragma once
#pragma warning(disable : 4786)
#endif

/************************************************************************
   AABB-related routines.  (that's axis-aligned bounding box)
                               
   - Adapted from code fragments by Tomas Miller and Dave Eberly

   $Id: aabb.h,v 1.1.1.1 2004/02/04 03:19:51 wensu Exp $
 ************************************************************************/

// #include <gfx/gfx.h>
// #include <gfx/gmVector3.h>
#include "libgm/gm.h"


bool aabb_intersects_plane(const gmVector3 &length, const gmVector3 &normal, double d);
bool aabb_intersects_triangle(const gmVector3 &center, const gmVector3 &length,
                              const gmVector3 &v0, const gmVector3 &v1, const gmVector3 &v2);
bool aabb_intersects_aabb(const gmVector3 &center1, const gmVector3 &length1,
                          const gmVector3 &center2, const gmVector3 &length2);
double aabb_distance_to_point(const gmVector3 &center, const gmVector3 &length,
                              const gmVector3 &p, gmVector3 *where);

#endif
