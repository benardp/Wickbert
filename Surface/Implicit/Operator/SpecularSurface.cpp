/**
 * Implementation of a specular implicit 
 * @file SpecularSurface.cpp
 * @date 16 March. 2006
 * @author Matei N. Stroila
 * @remarks
 */

#include "SpecularSurface.h"
#include "Particles/Attributes/ViewDependence.h"
#include "Particles/Attributes/LightPosition.h"

REGISTER_IMPLICIT(SpecularSurface,"SpecularSurface");

SpecularSurface::SpecularSurface(void)
{
	p = NULL;
	myF = NULL;
	new SurfParRefParam(this,&p,"<empty>","particles","Particles",
		"Collection of particles");
	new SurfImpRefParam(this,&myF,"<empty>","implicit","Implicit",
		"Implicit whom specular surface we want");
	new SurfParamDouble(this,&shine,0.9,"shine","size of specular",
					  "The arccos of this number determines the specular shadow curve size");
}

double SpecularSurface::proc(const gmVector3 & x)
{
	if(!p || !myF) return 0; 
	
	ViewDependence* view = p->getAttribute<ViewDependence>(std::string("ViewDependence"));
	LightPosition* light = p->getAttribute<LightPosition>(std::string("LightPosition"));
	gmVector3* myV = view->getCameraPosition();

	gmVector4 lightPos4;
	light->getLightPosition(lightPos4);
	
	gmVector3 lightPos3(lightPos4[0],lightPos4[1],lightPos4[2]);
		
	gmVector3 gradient = myF->grad(x);
	double gradLength = gradient.length();
		
	gmVector3 ci = *myV - x;
	ci.normalize();
	
	gmVector3 li = lightPos3 - x;
	li.normalize();

	gmVector3 si = li + ci;
	 
	double siLength = si.length(); 
				 
	return dot(myF->grad(x),si) - shine * gradLength * siLength;
}

/*
gmVector3 SpecularSurface::grad(const gmVector3 & x)
{
	if(!p || !myF) return  gmVector3();
	ViewDependence* view = p->getAttribute<ViewDependence>(std::string("ViewDependence"));
	LightPosition* light = p->getAttribute<LightPosition>(std::string("LightPosition"));
	gmVector3* myV = view->getCameraPosition();

	gmVector4 lightPos4;
	light->getLightPosition(lightPos4);
	bool directionalLight = light->isDirectional; //is a directional light?
	
	gmVector3 lightPos3;
	if(directionalLight)
		for(int j = 0; j <3 ; j++) lightPos3[j] = lightPos4[j];
	else
		for(int j = 0; j <3 ; j++) lightPos3[j] = lightPos4[j]/lightPos4[3];
		
	gmVector3 gradient = myF->grad(x);
	double gradLength = gradient.length();
		
	gmVector3 ci = *myV - x;
	ci.normalize();
	
	gmVector3 li = lightPos3 - x;
	li.normalize();	

	gmVector3 si = li + ci;
	 
	double siLength = si.length(); 
		
	gmMatrix3 Hess = myF->hess(x);
													  
	gmVector3 shadowNormal;	
		shadowNormal = (Hess * si) - 2 * gradient -  shine *(siLength * (1/gradLength) * (Hess * gradient)
														- 2 * gradLength * (1/siLength) * si);	
			 
	return shadowNormal;
}
*/
