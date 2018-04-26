/*! \file VisibleSurfaceAdhesion.cpp

 *  \author Jared Hoberock

 *  \brief Implementation of VisibleSurfaceAdhesion class.

 */



#include "VisibleSurfaceAdhesion.h"

REGISTER_PARTICLESTUFF(VisibleSurfaceAdhesion,"Behavior:VisibleSurfaceAdhesion");



/*! VisibleSurfaceAdhesion moves a particle toward the silhouette if it is on an invisible portion

 *  of the surface.

 */

void VisibleSurfaceAdhesion::applyForce(void)

{

	gmVector3 N, V;
	gmMatrix3 Hess;
  gmVector3 c = mViewDependence->getCameraPosition();
	for(unsigned int i= 0; i < ps->size(); i++) 
	{
		N = imp_int->getImplicit()->grad(position->getPosition(i));
		N.normalize();
		V = (c - position->getPosition(i));
		V.normalize();

    if(dot(N,V) < 0.0)
    {
		  Hess = imp_int->getImplicit()->hess(position->getPosition(i));
		  velocity->v[i] -= feedback * dot(N, V) * ((Hess * V)-N);
    } // end if
	}

} // end VisibleSurfaceAdhesion::applyForce()