/************************************************************************
   AABB code.  This should probably be developed into a proper object, but
   for now, it's more useful as a collection of low-level functions.
   This includes:
     - intersection tests with planes and triangles
       Adapted with very few changes from code by Tomas Möller.
     - distance to point calculation
       Adapted from OBB-point distance code by Dave Eberly.
 ************************************************************************/

#include "aabb.h"

#define X 0
#define Y 1
#define Z 2

#define FINDMINMAX(x0,x1,x2,min,max) \
  min = max = x0;   \
  if(x1<min) min=x1;\
  if(x1>max) max=x1;\
  if(x2<min) min=x2;\
  if(x2>max) max=x2;


bool aabb_intersects_plane(
    const gmVector3& length,
    const gmVector3& normal, 
    double      d )
{
    int q;
    gmVector3 vmin, vmax;
    for(q=X;q<=Z;q++)
    {
        if(normal[q]>0.0f)
        {
            vmin[q]=-length[q];
            vmax[q]=length[q];
        }
        else
        {
            vmin[q]=length[q];
            vmax[q]=-length[q];
        }
    }
    if(dot(normal,vmin)+d>0.0f) return false;
    if(dot(normal,vmax)+d>0.0f) return true;
  
    return false;
}


/*======================== X-tests ========================*/
#define AXISTEST_X01(a, b, fa, fb)			   \
	p0 = a*v0[Y] - b*v0[Z];			       	   \
	p2 = a*v2[Y] - b*v2[Z];			       	   \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
	rad = fa * length[Y] + fb * length[Z];   \
	if(min>rad || max<-rad) return false;

#define AXISTEST_X2(a, b, fa, fb)			   \
	p0 = a*v0[Y] - b*v0[Z];			           \
	p1 = a*v1[Y] - b*v1[Z];			       	   \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * length[Y] + fb * length[Z];   \
	if(min>rad || max<-rad) return false;

/*======================== Y-tests ========================*/
#define AXISTEST_Y02(a, b, fa, fb)			   \
	p0 = -a*v0[X] + b*v0[Z];		      	   \
	p2 = -a*v2[X] + b*v2[Z];	       	       	   \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
	rad = fa * length[X] + fb * length[Z];   \
	if(min>rad || max<-rad) return false;

#define AXISTEST_Y1(a, b, fa, fb)			   \
	p0 = -a*v0[X] + b*v0[Z];		      	   \
	p1 = -a*v1[X] + b*v1[Z];	     	       	   \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * length[X] + fb * length[Z];   \
	if(min>rad || max<-rad) return false;

/*======================== Z-tests ========================*/

#define AXISTEST_Z12(a, b, fa, fb)			   \
	p1 = a*v1[X] - b*v1[Y];			           \
	p2 = a*v2[X] - b*v2[Y];			       	   \
        if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
	rad = fa * length[X] + fb * length[Y];   \
	if(min>rad || max<-rad) return false;

#define AXISTEST_Z0(a, b, fa, fb)			   \
	p0 = a*v0[X] - b*v0[Y];				   \
	p1 = a*v1[X] - b*v1[Y];			           \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * length[X] + fb * length[Y];   \
	if(min>rad || max<-rad) return false;

bool aabb_intersects_triangle(
    const gmVector3& aabbcenter,
    const gmVector3& length,
    const gmVector3& vv0,
    const gmVector3& vv1, 
    const gmVector3& vv2 )
{
    // use separating axis theorem to test overlap between triangle and aabb 
    // need to test for overlap in these directions: 
    // 1) the {x,y,z}-directions 
    //    (actually, since we use the AABB of the triangle 
    //    we do not even need to test these) 
    // 2) normal of the triangle 
    // 3) crossproduct(edge from tri, {x,y,z}-directin) 
    //    this gives 3x3=9 more tests 

    gmVector3 v0, v1, v2;
    double min,max,d,p0,p1,p2,rad,fex,fey,fez;  
    gmVector3 normal, e0, e1, e2;

    // 1) first test overlap in the {x,y,z}-directions 
    //    find min, max of the triangle each direction, and test for 
    //    overlap in that direction -- this is equivalent to testing a 
    //    minimal AABB around the triangle against the AABB 
    v0 = vv0 - aabbcenter;
    v1 = vv1 - aabbcenter;
    v2 = vv2 - aabbcenter;

    FINDMINMAX(v0[X],v1[X],v2[X],min,max);
    if(min>length[X] || max<-length[X]) 
        return false;

    FINDMINMAX(v0[Y],v1[Y],v2[Y],min,max);
    if(min>length[Y] || max<-length[Y]) 
        return false;

    FINDMINMAX(v0[Z],v1[Z],v2[Z],min,max);
    if(min>length[Z] || max<-length[Z]) 
        return false;

    // 2) test if the aabb intersects the plane of the triangle 
    //    compute plane equation of triangle: normal*x+d=0 
    e0 = v1 - v0;
    e1 = v2 - v1;
    normal = cross(e0,e1);
    d=-dot(normal,v0);  // plane eq: normal.x+d=0 

    if(!aabb_intersects_plane(length, normal, d)) 
        return false;

    //    compute the last triangle edge 
    e2 = v0 - v2;

    //    3) 
    fex = fabs(e0[X]);
    fey = fabs(e0[Y]);
    fez = fabs(e0[Z]);
    AXISTEST_X01(e0[Z], e0[Y], fez, fey);
    AXISTEST_Y02(e0[Z], e0[X], fez, fex);
    AXISTEST_Z12(e0[Y], e0[X], fey, fex);

    fex = fabs(e1[X]);
    fey = fabs(e1[Y]);
    fez = fabs(e1[Z]);
    AXISTEST_X01(e1[Z], e1[Y], fez, fey);
    AXISTEST_Y02(e1[Z], e1[X], fez, fex);
    AXISTEST_Z0(e1[Y], e1[X], fey, fex);


    fex = fabs(e2[X]);
    fey = fabs(e2[Y]);
    fez = fabs(e2[Z]);
    AXISTEST_X2(e2[Z], e2[Y], fez, fey);
    AXISTEST_Y1(e2[Z], e2[X], fez, fex);
    AXISTEST_Z12(e2[Y], e2[X], fey, fex);

    return true;
}


// AABB-AABB intersection test, based on Gamasutra article by 
// Miguel Gomez; October 18, 1999

bool aabb_intersects_aabb(
    const gmVector3& center1,
    const gmVector3& length1,
    const gmVector3& center2,
    const gmVector3& length2 )
{
    gmVector3 diff = center2 - center1;
    gmVector3 ext = length1 + length2;

    return (    ( fabs(diff[0]) <= ext[0] )
             && ( fabs(diff[1]) <= ext[1] )
             && ( fabs(diff[2]) <= ext[2] ) );
}


// Point-AABB distance based on code from Dave Eberly's Magic Software:
//
// Magic Software, Inc.
// http://www.magic-software.com
// Copyright (c) 2000, 2001.  All Rights Reserved
//
// Source code from Magic Software is supplied under the terms of a license
// agreement and may not be copied or disclosed except in accordance with the
// terms of that agreement.  The various license agreements may be found at
// the Magic Software web site.  This file is subject to the license
//
// FREE SOURCE CODE
// http://www.magic-software.com/License/free.pdf

double aabb_distance_to_point(
    const gmVector3& center,
    const gmVector3& length,
    const gmVector3& p,
    gmVector3 *      where )
{
    // compute coordinates of point in aabb coordinate system
    gmVector3 diff = p - center;

    // project test point onto aabb
    double distance = 0.0;
    double delta;

    for( size_t i=0; i<3; i++ )
        if ( diff[i] < -length[i] )
        {
            delta = diff[i] + length[i];
            distance += delta*delta;
            diff[i] = -length[i];
        }
        else if ( diff[i] > length[i] )
        {
            delta = diff[i] - length[i];
            distance += delta*delta;
            diff[i] = length[i];
        }

    if ( where )
        *where = diff;

    return distance;
}
