/*! \file ParticleShaderRadialCurvature.h
 *  \author Jared Hoberock
 *  \brief Defines the interface to a class which shades particles proportional to radial curvature wrt camera direction.
 */

#ifndef PARTICLE_SHADER_RADIAL_CURVATURE_H
#define PARTICLE_SHADER_RADIAL_CURVATURE_H

#include "ParticleShader.h"
#include "ImplicitInterrogator.h"
#include "ViewDependence.h"

class ParticleShaderRadialCurvature: public ParticleShader
{

protected:

	ImplicitInterrogator *imp_int;
  ViewDependence *mViewDependence;
	double scale;

public:

	MAKE_PARTICLESTUFF_NAME();

	ParticleShaderRadialCurvature(Particles *ps=NULL);
	virtual void attachAttributes();
   
	virtual void drawShape(int i);

  /*! This method returns the "color" of particle i.
   *  \param i The particle to query against.
   *  \return Particle i's color.
   *  \author Jared Hoberock
   */
  gmVector4 getColor(int i);

	/// parameters
	int qlen();
	void getq(double *q);
	void setq(double *q);
	void qname(char **qn);

	char *tip() {return "Sets material color according to radial curvature with respect to direction to camera.";}
	char *qtip(int);

};

#endif // PARTICLE_SHADER_RADIAL_CURVATURE_H