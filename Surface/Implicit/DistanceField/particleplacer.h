#ifndef PARTICLEPLACER_H
#define PARTICLEPLACER_H

#include "particle.h"
#include "trianglemesh2.h"

class ParticlePlacer {
 public:
  static void CoverMeshWithParticles(const TriangleMesh2& mesh, std::vector< Particle >& particles, scalar radius, scalar spacing);					 
  static void FillMeshWithParticles(const TriangleMesh2& mesh, std::vector< Particle >& particles, scalar radius, scalar spacing);
  static void CoverMeshSurfaceWithParticles(const TriangleMesh2& mesh, std::vector< Particle >& particles, scalar radius, scalar spacing);
};

#endif
