/**
 * Implementation of a new implicit defined as the zero set of \f$ \nabla(F) \dot V \f$, 
 * where \f$ F \f$ is an implicit and  \f$ V \f$ is vector.
 * @file GradFdotV.cpp
 * @date 19 May. 2005
 * @author Matei N. Stroila
 * @remarks
 */

#include "GradFdotV.h"
#include "Particles/Attributes/ViewDependence.h"

REGISTER_IMPLICIT(GradFdotV,"GradFdotV");

GradFdotV::GradFdotV(void)
{
	p = NULL;
	myF = NULL;
	new SurfParRefParam(this,&p,"<empty>","particles","Particles",
		"Collection of particles");
	new SurfImpRefParam(this,&myF,"<empty>","implicit","Implicit",
		"Implicit whom silhouette we want");
}

double GradFdotV::proc(const gmVector3 & x)
{
	if(!p || !myF) return 0;

	ViewDependence* view = p->getAttribute<ViewDependence>(std::string("ViewDependence"));
	gmVector3* myV = view->getCameraPosition();
	return dot(myF->grad(x), ( *myV - x));
}

gmVector3 GradFdotV::grad(const gmVector3 & x)
{
	if(!p || !myF) return gmVector3();

	ViewDependence* view = p->getAttribute<ViewDependence>(std::string("ViewDependence"));
	gmVector3* myV = view->getCameraPosition();
	return myF->hess(x) *(*myV - x) - myF->grad(x); 
}



