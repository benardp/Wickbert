#ifndef DISTANCETRANSFORM_H
#define DISTANCETRANSFORM_H

#include "trianglemesh2.h"
#include "trianglenearestpoint.h"
#include "sparsescalarlattice.h"
#include "adf.h"

class SurfaceMesh;

class TriangleSignedDistance : private TriangleNearestPoint {
 public:
  Vector3d angle_weighted_normals[7];

  TriangleSignedDistance(const Triangle& tri) : TriangleNearestPoint(tri) {}
  TriangleSignedDistance(const gmVector3& v0, const gmVector3& v1, const gmVector3& v2) : TriangleNearestPoint(v0,v1,v2) {}

  scalar SignedDistanceTo(const Vector3d& p) const;
};

class TriangleSquaredDistance : private TriangleNearestPoint {
 public:
  TriangleSquaredDistance(const Triangle& tri) : TriangleNearestPoint(tri) {}

  scalar SquaredDistanceTo(const Vector3d& p) const;
};


class SignedDistanceTransform {
 private:
  std::vector< TriangleSignedDistance > tri_distances;
  //  std::vector< AABoundingBox > bounding_boxes;
  std::vector< Triangle > triangles;
    
 public:
  SignedDistanceTransform(const TriangleMesh2& mesh);
  //new wrapper for a surface mesh -Elmar
  SignedDistanceTransform(const SurfaceMesh& mesh);

  void Voxelize(SparseScalarLattice& lattice, scalar max_dist);
  void Voxelize(ADF& adf, SparseScalarLattice& lattice, scalar max_dist, int max_subdivisions, scalar abs_error);
  void Voxelize(ADF& adf, scalar max_dist, int max_subdivisions, scalar abs_error);
  
  double signedDistanceTo(const gmVector3 & pos) const;

  static void Test();
};


class UnsignedDistanceTransform { 
  private:
  std::vector< TriangleSquaredDistance > tri_distances;
  std::vector< Triangle > triangles;
 public:
  UnsignedDistanceTransform(const TriangleMesh2& mesh);
  void Voxelize(SparseScalarLattice& lattice, scalar max_dist);
};

#endif
