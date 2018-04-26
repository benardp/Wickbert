#ifndef BOUNDINGVOLUME_H
#define BOUNDINGVOLUME_H

#include "vector3d.h"
#include "renderable.h"

class BoundingVolume : Renderable {
 public:
  virtual void Render() {}
};

class AABoundingBox : public BoundingVolume{
 private:
  bool initialized;
 public:

  Vector3d v0;
  Vector3d v1;
  
  AABoundingBox() : initialized(false) {}
  AABoundingBox(const Vector3d& a, const Vector3d& b) : initialized(true), v0(a), v1(b) {}

  void include(const Vector3d& v);
  void include(const AABoundingBox& b);

  virtual void Render();
};


class BoundingSphere : public BoundingVolume {
 private:
  bool initialized;
 public:

  Vector3d center;
  scalar radius;

  BoundingSphere() : initialized(false) {}
  BoundingSphere(const Vector3d& c, const scalar r) : initialized(true), center(c), radius(r) {}

  void include(const Vector3d& v);
  void include(const BoundingSphere& b);
 
  virtual void Render();
};

#endif

