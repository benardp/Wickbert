/*! \file ParticleShaderContourHider.h
 *  \author Jared Hoberock
 *  \brief Defines the interface for a particle shader class which hides contour particles failing the derivative test.
 */

#ifndef PARTICLE_SHADER_CONTOUR_HIDER_H
#define PARTICLE_SHADER_CONTOUR_HIDER_H

#include "ParticleShader.h"

class ViewDependence;
class ImplicitInterrogator;
class ParticleVisibility;

/*! \class ParticleShaderContourHider
 *  \brief This class hides particles failing the radial curvature directional derivative test.
 *  \see http://www.cs.rutgers.edu/~decarlo/pubs/sg03.pdf
 *  \see SuggestiveContourAdhesion
 */
class ParticleShaderContourHider : public ParticleShader
{
  public:

    MAKE_PARTICLESTUFF_NAME();

	  ParticleShaderContourHider(Particles *ps=NULL);

	  virtual void attachAttributes();
   
	  virtual void setVisibility(int i);

	  char *tip() {return "Hides contour particles which fail the radial curvature derivative test.";}

  protected:
    /*! A ParticleShaderContourHider must be able to interrogate the surface.
     */
  	ImplicitInterrogator *mImplicitInterrogator;

    /*! A ParticleShaderContourHider must have access to the current view.
     */
    ViewDependence *mViewDependence;

	ParticleVisibility *mVisibility;

	std::string imp_name, view_name, vis_name;

	ParticlePosition *mPosition;
}; // end class ParticleShaderContourHider

#endif // PARTICLE_SHADER_CONTOUR_HIDER_H
