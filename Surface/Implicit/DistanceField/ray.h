#ifndef RAY_H
#define RAY_H

#include "vector3d.h"

class Ray{
 public:
  Vector3d o;
  Vector3d d;

  Ray(){
    o = Vector3d(0,0,0);
    d = Vector3d(0,0,1);
  }

  Ray(Vector3d origin, Vector3d direction){
    o = origin;
    d = direction;
  }


  Vector3d EvalAt(scalar x) const{
    return o + x * d;
  }

  void drawOpenGL() const;

};

#endif
