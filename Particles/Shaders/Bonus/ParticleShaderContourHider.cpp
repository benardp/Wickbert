/*! \class ParticleShaderContourHider.cpp
 *  \author Jared Hoberock
 *  \brief Implementation of ParticleShaderContourHider class.
 */

#include "ParticleShaderContourHider.h"
#include "ViewDependence.h"
#include "ImplicitInterrogator.h"
#include "ParticleVisibility.h"

REGISTER_PARTICLESTUFF(ParticleShaderContourHider,"Shader:ParticleShaderContourHider");

ParticleShaderContourHider::ParticleShaderContourHider(Particles *ps)
  :ParticleShader(ps,std::string("ParticleShaderContourHider"))
{
	mImplicitInterrogator = 0;
	mViewDependence = 0;
	mVisibility = 0;

	new PSParamString(this,&imp_name,"ImplicitInterrogator",
		"impname","Implicit Interrogator","Name of the Implicit Interrogator attribute.");
	new PSParamString(this,&view_name,"ViewDependence",
		"viewname","View Dependence","Name of the View Dependence attribute.");
	new PSParamString(this,&vis_name,"ParticleVisibility",
		"visname","Visibility","Name of the Visibility attribute.");
} // end ParticleShaderContourHider::ParticleShaderContourHider()

void ParticleShaderContourHider::attachAttributes()
{
	ParticleShader::attachAttributes();
	attachAttribute(mImplicitInterrogator,imp_name);
	attachAttribute(mViewDependence,view_name);
	attachAttribute(mVisibility,vis_name);
	attachAttribute(mPosition,std::string("ParticlePosition"));
} // end ParticleShaderContourHider::attachAttributes()

/*! We need to evaluate the expression Dw Kr.  If it is not positive, then we hide the
 *  particle.  Otherwise, we set it visible.
 *  \param i The index of the particle of interest.
 */
void ParticleShaderContourHider::setVisibility(int i)
{
	Implicit *imp = mImplicitInterrogator->getImplicit();

	if (!imp) return;

	// get the camera position
	gmVector3 c(0,0,0);
	if(mViewDependence)
	{
		c = mViewDependence->getCameraPosition();
	} // end if

	// get the position of this particle
	gmVector3 x = mPosition->getPosition(i);

	// get the normal at x
	gmVector3 n = imp->normal(x);

	// get a vector pointing toward the camera
	gmVector3 v = c - x;

	// project v onto the tangent plane to get w
	// w = v - n*(dot(n,v))
	gmVector3 w = v - dot(n,v)*n;

	// find the derivative of normal curvature in direction w in the direction
	// of w
	double d = dot(imp->gradNormalCurvature(x,w), w);

	// visible only if the derivative is positive
	mVisibility->getVisibility(i) = (d > 0)?(0.0):(1.0);

	//ParticleShader::mParticleVisibility->getVisibility
} // end ParticleShaderContourHider::drawShape()
