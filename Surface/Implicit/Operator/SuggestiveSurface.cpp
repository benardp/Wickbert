/*
 *  SuggestiveSurface.cpp
 *  gps
 *
 *  Created by Matei Stroila on 6/20/05.
 *  
 *
 */

#include "Surface/Implicit/Implicit.h"
#include "Particles/Attributes/ViewDependence.h"
#include "libgm/gm.h"
#include "SuggestiveSurface.h"

REGISTER_IMPLICIT(SuggestiveSurface,"SuggestiveSurface");

SuggestiveSurface::SuggestiveSurface(Implicit *F,  gmVector3 *cameraP)
{
	p = NULL;
	myF = F;
	myV = cameraP;
}

SuggestiveSurface::SuggestiveSurface(void)
{
	p = NULL;
	myF = NULL;
	new SurfParRefParam(this,&p,"<empty>","particles","Particles",
						"Collection of particles");
	new SurfImpRefParam(this,&myF,"<empty>","implicit","Implicit",
						"Implicit whom silhouette we want");
}

double SuggestiveSurface::proc(const gmVector3 & x)
{
	if(!(p && myF)) return 0.0; 
	view = p->getAttribute<ViewDependence>(std::string("ViewDependence"));
	myV = view->getCameraPosition();
	gmVector3 n(myF->grad(x));
	n.normalize();
	gmVector3 w = (*myV - x) - dot(n,(*myV - x)) * n;
	return myF->normalCurvature(x, w);
}

//TO DO: implement this
Intervald SuggestiveSurface::proc(const Box<double>& x)
{
	return Intervald();
}

gmVector3 SuggestiveSurface::grad(const gmVector3 & x)
{
	if(!(p && myF)) return gmVector3(); 
	view = p->getAttribute<ViewDependence>(std::string("ViewDependence"));
	myV = view->getCameraPosition();
	gmVector3 n(myF->grad(x));
	n.normalize();
	gmVector3 w = (*myV - x) - dot(n,(*myV - x)) * n;
	return myF->gradNormalCurvature(x, w);
}

 void SuggestiveSurface::_setq(double* q) 
 {
	 (*myV)[0] = q[0];
	 (*myV)[1] = q[1];
	 (*myV)[2] = q[2];
 }
 
 
 