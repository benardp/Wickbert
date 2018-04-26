#ifndef TRIANGLEBOXOVERLAP_H
#define TRIANGLEBOXOVERLAP_H

#include "triangle.h"

class TriangleBoxOverlap{
 private:
  Vector3d corners[3];

  Vector3d edge01;
  Vector3d edge12;
  Vector3d edge20;

  Vector3d bbox_min;
  Vector3d bbox_max;

  Vector3d normal;

///  AABoundingBox bbox;
 public:
  TriangleBoxOverlap(const Triangle&);

  bool Overlaps(const Vector3d& box_center, const Vector3d& boxhalfsize);

  static bool Test();
};


#endif
