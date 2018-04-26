/*! \file ParticleShaderViewpointStability.cpp
 *  \brief Implementation of ParticleShaderViewpointStability class. 
 */

#include "ParticleShaderViewpointStability.h"
#include "ParticleMaterial.h"
#include "ImplicitInterrogator.h"
#include "ViewDependence.h"

REGISTER_PARTICLESTUFF(ParticleShaderViewpointStability,"Shader:ParticleShaderViewpointStability");

ParticleShaderViewpointStability::ParticleShaderViewpointStability(Particles *ps)
  :ParticleShader(ps,std::string("ParticleShaderViewpointStability"))
{
	mImplicitInterrogator = 0;
  mViewDependence = 0;
  setThreshold(15);
} // end ParticleShaderViewpointStability::ParticleShaderViewpointStability()

void ParticleShaderViewpointStability::attachAttributes()
{
	ParticleShader::attachAttributes();
	attachAttribute(mImplicitInterrogator,std::string("ImplicitInterrogator"));
	attachAttribute(mViewDependence,std::string("ViewDependence"));
	attachAttribute(material,std::string("ParticleMaterial"));
} // end ParticleShaderViewpointStability::attachAttributes()

/// parameters
int ParticleShaderViewpointStability::qlen()
{
	return 1;
} // end ParticleShaderViewpointStability::qlen()

void ParticleShaderViewpointStability::getq(double *q)
{
	q[0] = getThreshold();
} // end ParticleShaderViewpointStability::getq()

void ParticleShaderViewpointStability::setq(double *q)
{
  setThreshold(q[0]);
} // end ParticleShaderViewpointStability::setq()

void ParticleShaderViewpointStability::qname(char **qn)
{
	qn[0] = "Viewpoint stability threshold.";
} // end ParticleShaderViewpointStability::qname()

char *ParticleShaderViewpointStability::qtip(int i)
{
	switch(i)
  {
		default:
		case 0:
			return "Viewpoint stability threshold.";
	} // end switch
} // end ParticleShaderViewpointStability::qtip()

gmVector4 ParticleShaderViewpointStability::getColor(int i)
{
  gmVector4 result(0,0,0,0);

	if (!mImplicitInterrogator) return result;
	if (!material) return result;

	Implicit *imp = mImplicitInterrogator->getImplicit();

	if (!imp) return result;

  // get the camera position
  gmVector3 c(0,0,0);
  if(mViewDependence)
  {
    c = mViewDependence->getCameraPosition();
  } // end if

  // get the position of this particle
  gmVector3 x = position->getPosition(i);

  // get the normal at x
  gmVector3 n = imp->normal(x);

  // get a vector pointing toward the camera
  gmVector3 v = c - x;
  v.normalize();

  // find the angle between n and v
  // note this is necessarily between 0 and 1
  double unstable = dot(v,n);
  double degrees = gmDegrees(fabs(acos(unstable)));

  // don't evaluate if degrees is > threshold
  if(degrees <= getThreshold())
  {
    // threshold is at 2 standard deviations
    double standardDeviation = getThreshold() / 2.0;
    unstable = exp(-(degrees*degrees)/(2.0*standardDeviation*standardDeviation));
  } // end if
  else
  {
    unstable = 0;
  } // end else

  // large angles are considered stable; around 0, unstable
  // interpolate between green (stable) and red (unstable)
	gmVector4 color = gmVector4(unstable, 1.0 - unstable, 0, 1.0);
  return color;
} // end ParticleShaderViewpointStability::getColor()

void ParticleShaderViewpointStability::drawShape(int i)
{
  gmVector4 color = getColor(i);

	material->setDiffuseFront(i,color);
	material->setDiffuseBack(i,0.5*color);
} // end ParticleShaderViewpointStability::drawShape()
