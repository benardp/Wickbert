/*! \file VisibleSurfaceAdhesion.h
 *  \author Jared Hoberock
 *  \brief Defines the interface for a behavior to constrain particles to visible surfaces.
 */

#ifndef VISIBLE_SURFACE_ADHESION_H
#define VISIBLE_SURFACE_ADHESION_H

#include "Particles.h"
#include "ParticleBehavior.h"
#include "ImplicitInterrogator.h"
#include "ParticleOrientation.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"
#include "SilhouetteAdhesion.h"
#include "ViewDependence.h"

#include "Surface/Implicit/Implicit.h"

/*! \class VisibleSurfaceAdhesion
 *  \brief VisibleSurfaceAdhesion is a particle behavior that attempts to kick particles off portions of the surface
 *         where n.v < 0.
 */
class VisibleSurfaceAdhesion : public SilhouetteAdhesion
{
  public:
    MAKE_PARTICLESTUFF_NAME();

	  /// Creates a surface adhesion attribute for Particles p.
	  inline VisibleSurfaceAdhesion(Particles *ps=NULL, const std::string& name=std::string("VisibleSurfaceAdhesion")):SilhouetteAdhesion(ps,name)
    {
      ;
    } // end VisibleSurfaceAdhesion::VisibleSurfaceAdhesion()

	  /** Apply the visible surface adhesion constraint.
	   */
	  void applyForce();
}; // end class VisibleSurfaceAdhesion

#endif // VISIBLE_SURFACE_ADHESION_H1