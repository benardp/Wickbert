/*! \file ParticleShaderRadialCurvature.cpp
 *  \author Jared Hoberock
 *  \brief Implementation of ParticleShaderRadialCurvature class.
 */

#include "ParticleShaderRadialCurvature.h"
#include "ParticleMaterial.h"

#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif // min

REGISTER_PARTICLESTUFF(ParticleShaderRadialCurvature,"Shader:ParticleShaderRadialCurvature");

ParticleShaderRadialCurvature::ParticleShaderRadialCurvature(Particles *ps)
  :ParticleShader(ps,std::string("ParticleShaderRadialCurvature"))
{
	imp_int = 0;
	material=0;
	scale = 0.00015;
  mViewDependence = 0;
} // end ParticleShaderRadialCurvature::ParticleShaderRadialCurvature()

void ParticleShaderRadialCurvature::attachAttributes()
{
	ParticleShader::attachAttributes();
	attachAttribute(imp_int,std::string("ImplicitInterrogator"));
  attachAttribute(mViewDependence,std::string("ViewDependence"));
} // end ParticleShaderRadialCurvature::attachAttributes()

/// parameters
int ParticleShaderRadialCurvature::qlen()
{
	return 1;
} // end ParticleShaderRadialCurvature::qlen()

void ParticleShaderRadialCurvature::getq(double *q)
{
	q[0] = scale;
} // end ParticleShaderRadialCurvature::getq()

void ParticleShaderRadialCurvature::setq(double *q)
{
  scale = q[0];
} // end ParticleShaderRadialCurvature::setq()

void ParticleShaderRadialCurvature::qname(char **qn)
{
	qn[0] = "Curvature Scale";
} // end ParticleShaderRadialCurvature::qname()

char *ParticleShaderRadialCurvature::qtip(int i)
{
	switch(i)
  {
		default:
		case 0:
			return "Curvature to [0,1] color scaling factor";
	} // end switch
} // end ParticleShaderRadialCurvature::qtip()

gmVector4 ParticleShaderRadialCurvature::getColor(int i)
{
  gmVector4 result(0,0,0,0);

	if (!imp_int) return result;
	if (!material) return result;

	Implicit *imp = imp_int->getImplicit();

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

  // project v onto the tangent plane to get w
  // w = v - n*(dot(n,v))
  gmVector3 w = v - dot(n,v)*n;

  // find the normal curvature in the direction of w at x
	double kr = imp->normalCurvature(x,w);

	kr *= scale;

	if (kr > 1.0) kr = 1.0;
	if (kr < -1.0) kr = -1.0;

  double convexity = 0, concavity = 0, flat = 0;
  if(kr < 0.0)
  {
    concavity = -kr;
    flat = 1.0 - concavity;
  } // end if
  else if(kr > 0.0)
  {
    convexity = kr;
    flat = 1.0 - convexity;
  } // end else if

	gmVector4 color = gmVector4(convexity, flat, concavity, 1.0);
  return color;
} // end ParticleShaderRadialCurvature::getColor()

void ParticleShaderRadialCurvature::drawShape(int i)
{
  gmVector4 color = getColor(i);

	material->setDiffuseFront(i,color);
	material->setDiffuseBack(i,0.5*color);
} // end ParticleShaderRadialCurvature::drawShape()
