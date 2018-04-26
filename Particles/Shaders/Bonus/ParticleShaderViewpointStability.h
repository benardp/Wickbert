/*! \file ParticleShaderViewpointStability.h
 *  \author Jared Hoberock
 *  \brief This file defines the interface to a shader class visualizing viewpoint stability for suggestive contours.
 */

#ifndef PARTICLE_SHADER_VIEWPOINT_STABILITY_H
#define PARTICLE_SHADER_VIEWPOINT_STABILITY_H

#include "ParticleShader.h"
class ImplicitInterrogator;
class ViewDependence;
class ParticleMaterial;

/*! \class ParticleShaderViewpointStability
 *  \brief This class visualizes the viewpoint stability of suggestive contour on regions of a surface.
 *  \see http://www.cs.rutgers.edu/~decarlo/contour.html
 */
class ParticleShaderViewpointStability : public ParticleShader
{
  public:

    MAKE_PARTICLESTUFF_NAME();

	  ParticleShaderViewpointStability(Particles *ps=NULL);
	  virtual void attachAttributes();

    /*! This method renders particle i using OpenGL.
     *  \param i The index of the particle of interest.
     */
	  virtual void drawShape(int i);

    /*! This method returns the "color" of particle i.
     *  \param i The particle to query against.
     *  \return Particle i's color.
     */
    gmVector4 getColor(int i);

    /*! This method returns mThreshold, the viewpoint stability threshold, in degrees.
     *  \return mThreshold
     */
    inline double getThreshold(void) const
    {
      return mThreshold;
    } // end ParticleShaderViewpointStability::getThreshold()

    /*! This method sets mThreshold.
     *  \param t Sets mThreshold
     */
    inline void setThreshold(const double t)
    {
      mThreshold = t;
    } // end ParticleShaderViewpointStability::setThreshold()

	  /// parameters
	  int qlen();
	  void getq(double *q);
	  void setq(double *q);
	  void qname(char **qn);

	  char *tip() {return "Sets material color according to angle between normal and direction to camera.";}
	  char *qtip(int);

  protected:
    /// ParticleShaderViewpointStability requires access to a surface
    ImplicitInterrogator *mImplicitInterrogator;

    /// ParticleShaderViewpointStability requires access to the view
    ViewDependence *mViewDependence;

	ParticleMaterial *material;

    /*! This is the threshold defined by DeCarlo in the SIGGRAPH paper.  It is
     *  specified in degrees.
     */
    double mThreshold;
}; // end class ParticleShaderViewpointStability 

#endif // PARTICLE_SHADER_VIEWPOINT_STABILITY_H
