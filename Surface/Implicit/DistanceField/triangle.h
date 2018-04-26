#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "boundingvolume.h"
#include "vector3d.h"
#include "ray.h"
#include "renderable.h"


class TriangleException {};
class DegenerateTriangleException : public TriangleException {};

class Triangle : public Renderable {
  friend std::ostream &operator<<(std::ostream&, const Triangle&);
  friend std::istream &operator>>(std::istream &input, Triangle& v);
 protected:
  Vector3d corners[3];
 public:
  Triangle(){}
  Triangle(const Vector3d& v0, const Vector3d& v1, const Vector3d& v2);
  //Not used and normals are not stored !!! - Elmar
  //Triangle(const Vector3d& v0, const Vector3d& v1, const Vector3d& v2, const Vector3d& n0, const Vector3d& n1, const Vector3d& n2);

  Vector3d Bary(scalar a, scalar b, scalar c) const;
  scalar Perimeter() const;
  scalar SurfaceArea() const;
  Vector3d RandomSample() const;

  const Vector3d& A() const {return corners[0]; }
        Vector3d& A()       {return corners[0]; }
  const Vector3d& B() const {return corners[1]; }
        Vector3d& B()       {return corners[1]; }
  const Vector3d& C() const {return corners[2]; }
        Vector3d& C()       {return corners[2]; }

  const Vector3d& operator[](int i) const{ return corners[i]; }
        Vector3d& operator[](int i)      { return corners[i]; }


	
  // The non-unitized normal of this triangle
  Vector3d Normal() const;

  bool Intersects(const Ray& ray) const;
  bool Intersection(const Ray& ray, scalar& t) const;

  AABoundingBox BoundingBox() const;

  void Translate(const Vector3d& v);

  virtual void Render();
};

inline std::ostream &operator<<(std::ostream &output, const Triangle& t){
  output << t.corners[0] << " " << t.corners[1] << " " << t.corners[2] << " ";
  //  output << t.normals[0] << " " << t.normals[1] << " " << t.normals[2];

  return output;
}

inline std::istream &operator>>(std::istream &input, Triangle& t){
  input >> t.corners[0];
  input >> t.corners[1];
  input >> t.corners[2];
  return input;
}

#endif

