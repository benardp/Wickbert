/*! \file SuggestiveContourStabilization.cpp

 *  \author Jared Hoberock

 *  \brief Implementation of SuggestiveContourStabilization class.

 */

#include "ImplicitInterrogator.h"

#include "ViewDependence.h"

#include "ParticlePosition.h"

#include "ParticleVelocity.h"

#include "SuggestiveContourStabilization.h"

REGISTER_PARTICLESTUFF(SuggestiveContourStabilization,"Behavior:SuggestiveContourStabilization");




SuggestiveContourStabilization::SuggestiveContourStabilization(Particles *theps,const std::string& name, double threshold, double scaleStep)
	: ParticleBehavior(theps,name)
{
		new Attached<ImplicitInterrogator>(this,&mImplicitInterrogator);
		new Attached<ViewDependence>(this,&mViewDependence);
		new Attached<ParticlePosition>(this,&position);
		new Attached<ParticleVelocity>(this,&velocity);
		setThreshold(threshold);
		setScaleStep(scaleStep);
} // end SuggestiveContourStabilization::SuggestiveContourStabilization()

// Can probably remove this now
void SuggestiveContourStabilization::attachAttributes()
{
	ParticleBehavior::attachAttributes();
} // end SuggestiveContourStabilization::attachAttributes()

int SuggestiveContourStabilization::qlen()
{
	return 2;
} // end SuggestiveContourStabilization::qlen()

void SuggestiveContourStabilization::getq(double *q)
{
	q[0] = getThreshold();
  q[1] = getScaleStep();
} // end SuggestiveContourStabilization::getq()

void SuggestiveContourStabilization::setq(double *q)
{
	setThreshold(q[0]);
  setScaleStep(q[1]);
} // end SuggestiveContourStabilization::setq()

void SuggestiveContourStabilization::qname(char **qn)
{
	qn[0] = "Stability threshold.";
  qn[1] = "Scale step size.";
} // end SuggestiveContourStabilization::qname()

/*! To stabilize suggestive contours, we check the angle between the normal and the view vector.  If they are near parallel,
 *  we must apply a force.  Near parallel is defined as when the angle between the normal and view vector is within the 
 *  threshold.  The weight of the force is gaussian weighted... in other words, exactly parallel will yield a greater force
 *  than near parallel.  The direction of the force applied is simply aimed at the silhouette.
 */
void SuggestiveContourStabilization::applyForce()
{
  // get the camera position
  const gmVector3 c = *(mViewDependence->getCameraPosition());
  gmVector3 x;
  gmVector3 v;
  gmVector3 n;
  gmVector3 step;
  gmMatrix3 hess;

  for(unsigned int i = 0; i < ps->size(); ++i)
  {
    x = position->getPosition(i);
    v = c - x;
    v.normalize();
    n = mImplicitInterrogator->normal(i);//->getImplicit()->normal(x);

    double degrees = gmDegrees(acos(dot(n,v)));
    if(degrees <= getThreshold())
    {
      double standardDeviation = getThreshold() / 2.0;
      double stepSize = exp(-(degrees*degrees)/(2.0*standardDeviation*standardDeviation));

		  hess = mImplicitInterrogator->getImplicit()->hess(x);
		  step = -((hess * v)-n);

      step.normalize();
      step *= getScaleStep()*stepSize;

      // update velocity
      velocity->v[i] += step;
    } // end if
  } // end for i
} // end SuggestiveContourStabilization::applyForce()