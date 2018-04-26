#ifndef TRIANGLENEARESTPOINT_H
#define TRIANGLENEARESTPOINT_H


#include "triangle.h"
//I don't want to change this class too much
//maybe someone still uses it...
//but to make it compatible to the new vector type
//I added this conversion - Elmar
#include "libgm/gm.h"

class TriangleNearestPoint{
public:
	enum VoronoiRegion { CORNER_A_REGION = 0, CORNER_B_REGION = 1, EDGE_AB_REGION = 2, 
		    CORNER_C_REGION = 3, EDGE_AC_REGION = 4, EDGE_BC_REGION = 5, 
		    FACE_REGION = 6, INVALID_REGION = 7 };

private:
  void initFromVertices(const Vector3d &a, const Vector3d &b, const Vector3d &c);
  Vector3d corners[3];

  Vector3d unitnormal;

  Vector3d span01;
  Vector3d span02;
  Vector3d span12;
  scalar u11, u22, u12, u33, u13;

  scalar d_inv;

  Vector3d u_span01, u_span02, u_span12;

 public:
  TriangleNearestPoint(const Triangle&);
  //Lousy wrapper ... see next comment... - Elmar
  TriangleNearestPoint(const gmVector3& a, const gmVector3& b, const gmVector3& c);

  Vector3d NearestPointTo(const Vector3d& p) const;

  //changed the type of the region to Voronoi as it should be. -Elmar
  Vector3d NearestPointToVerbose(const Vector3d& p, VoronoiRegion & region) const;

  //this is just a wrapper... but I did not want to add a new class doing the same job,
  //especially because this one is really nicely implemented !!!! 
  //On the other hand I did not want to reintroduce the Vector3d type in Wickbert.
  //The function performs two unnecessary conversions and I don't like the function type 
  //but wanted to stay close to the original... anyway performance is not important for this class
  //and it is mostly a question of personal taste ;) - Elmar 
  gmVector3 NearestPointTo(const gmVector3 & pos, VoronoiRegion & region) const;


  static bool Test();
};


#endif
