/*
@file ParticlesInclude.h
@author Wen Su
This class contains all the particles header files. So it is easy to sue the particle system.
*/
#ifndef PARTICLEINCLUDE_H
#define PARTICLEINCLUDE_H

// basics
#include "Particles.h"
#include "ParticleStuff.h"
#include "ParticleAttribute.h"
#include "ParticleBehavior.h"
#include "ParticleShader.h"
#include "ParticleSystem.h"
#include "ParFileManager.h"

// attributes
#include "ParticleAge.h"
#include "AdaptiveRepulsionData.h"
#include "ParticleLocality.h"
#include "ParticleLocalityGrid.h"
#include "ParticleEigen.h"
#include "ImplicitInterrogator.h"
#include "ParticleMesh.h"
#include "ParticleNormal.h"
#include "ParticleOrientation.h"
#include "ParticleDesiredOrientation.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"
#include "ParticleMaterial.h"
#include "SurfaceInterrogator.h"
#include "ParticleBoundingBox.h"
#include "ViewDependence.h"
#include "ParticleVisibility.h"

// behaviors
#include "SurfaceAdhesion.h"
#include "SurfaceDeformation.h"
#include "ParticleFate.h"
#include "SkinMeshShape.h"
#include "ParticleCreation.h"
#include "MeshShape.h"
#include "TopologyInterrogator.h"
#include "SilhouetteAdhesion.h"
#include "SilhouetteFate.h"
#include "SingularityAdhesion.h"
#include "SingularityInterrogator.h"
#include "SingularityRepulsion.h"
#include "VelocityIntegrator.h"
#include "CriticalPointInterrogator.h"
#include "ClusterMeshInterrogator.h"
#include "KeepInBounds.h"
#include "ParticleRepulsion.h"
#include "ParticleValueDeath.h"
#include "SuggestiveContourAdhesion.h"
#include "SuggestiveContourStabilization.h"
#include "VisibleSurfaceAdhesion.h"
#include "NeighborAlignment.h"

// shaders
#include "ParticleShaderLink.h"
#include "ParticleShaderDisk.h"
#include "ElmarsShaderTest.h"
#include "ParticleShaderCylinder.h"
#include "ParticleShaderChain.h"
#include "ParticleShaderStroke.h"
#include "ParticleShaderTriangle.h"
#include "ParticleShaderGradientField.h"
#include "ParticleShaderSphere.h"
#include "ParticleShaderCluster.h"
#include "ParticleShaderRadialCurvature.h"
#include "ParticleShaderContourHider.h"
#include "ParticleShaderViewpointStability.h"
#include "ParticleShaderOriented.h"

#endif