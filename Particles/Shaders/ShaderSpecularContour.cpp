/**
 * Implementation of a specular contour shader
 * @file ShaderSpecularContour.cpp
 * @date 04/07/2006
 * @author Matei N. Stroila
 * @remarks
 */

#include "ShaderSpecularContour.h"
#include "LightPosition.h"
#include "Attributes/ParticlePosition.h"
#include "Attributes/ParticleVelocity.h"
#include "Attributes/ImplicitInterrogator.h"
#include "Attributes/ViewDependence.h"
#include "Attributes/ParticleVisibility.h"
#include "Attributes/ParticleLocality.h"
#include "Attributes/AdaptiveRepulsionData.h"
#include "Attributes/ParticleVector.h"
#include "Attributes/ParticleBoundingBox.h"


REGISTER_PARTICLESTUFF(ShaderSpecularContour,"Shader:ShaderSpecularContour");

ShaderSpecularContour::ShaderSpecularContour(Particles *ps)
:ShaderContour(ps,std::string("ShaderSpecularContour"))
{
	
	new PSParamDouble(this,&shine,0.9,"shine","size of specular",
					  "The arccos of this number determines the specular shadow curve size");
	new Attached<LightPosition>(this,&light);
	brushColorf = gmVector3(1.0,1.0,1.0);
}

void ShaderSpecularContour::findTangent(const unsigned int i)
{
	gmVector3 myV = *(view->getCameraPosition());

	gmVector3 ni, posi; //particle normal and position

	posi = position->getPosition(i);
	ni = impInt->grad(i);

	gmVector4 lightPos4;
	light->getLightPosition(lightPos4);
	gmVector3 lightPos3(lightPos4[0],lightPos4[1],lightPos4[2]);	
		
	gmVector3 s = myV.normalize() + lightPos3.normalize();
	
	double niLength = ni.length();
		
	double sLength = s.length(); 
		
	gmMatrix3 Hess = impInt->hess(i);
													  
	gmVector3 specularNormal;	

	specularNormal = (Hess * s) - shine *( sLength * (1/niLength) *( Hess * ni));	

	gmVector3 tangent = -1 * cross(ni,specularNormal);

	if(tangent.lengthSquared() > gmEPSILON) tangent.normalize();
	
	tangentAttr->setVector(i,tangent);
	
}

ShaderSpecularContour::~ShaderSpecularContour(){};
