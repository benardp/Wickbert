/**
 * @file ParticleMaterial.h
 * @author Terry Fleury
 * This attribute stores per particle material properties
 */

#ifndef PARTICLEMATERIAL_H
#define PARTICLEMATERIAL_H

#include "Particles.h"
#include "ParticleAttribute.h"

class ParticleMaterial : public ParticleAttribute {

protected:

public:

  MAKE_PARTICLESTUFF_NAME();

  /// per particle diffuse colors
  std::vector<gmVector4> diffuseFront;
  std::vector<gmVector4> diffuseBack;

  gmVector4 all_diffuseFront;
  gmVector4 all_diffuseBack;

  gmVector4 sel_diffuseFront;
  gmVector4 sel_diffuseBack;

  /// Add per particle material properties to a system of particles.
  ParticleMaterial(Particles *ps=NULL, const std::string& name=std::string("ParticleMaterial"));

  void setParticleSystem(Particles *);

  /// Callback for particle addition.
  void particleAdded();

  /// Callback for particle removal.
  void particleRemoved(unsigned int i);

  virtual void clear();

  virtual void setMaterial(unsigned int i);
};

#endif

